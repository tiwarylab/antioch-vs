#!/bin/bash
#
# setup_MD_from_PDB.sh — build an AMBER MD system (topology + coordinates) from a receptor PDB.
# Used for both the AF2RAVE SPIB MD and the CTMD MD.
#
# Common prep (every receptor, soluble or membrane):
#   clean_pdb.py -> PROPKA protonation -> Reduce/CPPTRAJ prepareforleap
#   -> disulfide detection + required-disulfide eligibility gate (skip ineligible cleanly)
# Then the box is built, differing only in how waters/ions are placed:
#   soluble  (default)                      -> tleap setup.leap (solvateoct TIP3P + neutralizing ions)
#   membrane (/mnt/receptor_files/MEMBRANE) -> pack_membrane (POPC bilayer + water + 0.15 M salt)
#                                              -> tleap membrane_setup.leap
# Both build the topology in tleap with the CPPTRAJ-detected disulfide bonds inserted.
#
# The membrane mode is selected by the /mnt/receptor_files/MEMBRANE marker (set by `init -m`).
#
# Target-specific prep parameters are environment variables with sensible defaults (nothing hardcoded to
# one target); see ctmd_module/md_submodule/README.md:
#   DOCKFLOW_PH (7.0), DOCKFLOW_PKA_CUTOFF (2.0), DOCKFLOW_PROTECT_RES (""), DOCKFLOW_REQUIRED_DISULFIDES
#   ("" => no gate), DOCKFLOW_DISULFIDE_SG_CUTOFF (2.3), DOCKFLOW_REDUCE_EXE (reduce), DOCKFLOW_CPPTRAJ_EXE
#   (cpptraj). The image must provide propka3 + AmberTools (reduce/cpptraj/tleap) and, for membranes,
#   packmol-memgen.
#
# Force field: ff14SB + GAFF2 + TIP3P (+ Lipid21 for membranes).

pdbfile="${1:-}"
if test -z "$pdbfile"
then
    echo "Please choose a PDB file to start from"
    exit 1
fi

if test -f "$pdbfile"
then
    echo "PDB file $pdbfile found"
else
    echo "PDB file $pdbfile NOT found"
    exit 1
fi

output_folder="MD_setup"
if test -n "${2:-}"
then
    output_folder="$2"
fi
echo "Results will be saved to $output_folder"

# ─── Paths ────────────────────────────────────────────────────────────────────
# DOCKFLOW_ROOT defaults to the in-SIF location but can be overridden for a bare-metal checkout.
DOCKFLOW_ROOT="${DOCKFLOW_ROOT:-/usr/software/dockflow/}"
RECEPTOR_PREP="$DOCKFLOW_ROOT/docking_module/receptor_prep"

# Bare-metal runs point DOCKFLOW_ENV_SETUP at a script that sources AmberTools (tleap/cpptraj/reduce),
# propka3, and packmol-memgen onto PATH; unset inside the SIF, where the tools are already on PATH.
# Sourced before `set -u` so an unbound var in the env script does not abort the run.
if [ -n "${DOCKFLOW_ENV_SETUP:-}" ]; then source "$DOCKFLOW_ENV_SETUP"; fi

# ─── Prep config (env with defaults; nothing hardcoded to one target) ───────────
PH="${DOCKFLOW_PH:-7.0}"
PKA_CUTOFF="${DOCKFLOW_PKA_CUTOFF:-2.0}"
SG_CUTOFF="${DOCKFLOW_DISULFIDE_SG_CUTOFF:-2.3}"
PROTECT_RES="${DOCKFLOW_PROTECT_RES:-}"
REQUIRED_DISULFIDES="${DOCKFLOW_REQUIRED_DISULFIDES:-}"
REDUCE_EXE="${DOCKFLOW_REDUCE_EXE:-reduce}"
CPPTRAJ_EXE="${DOCKFLOW_CPPTRAJ_EXE:-cpptraj}"
# Skip CPPTRAJ prepareforleap's sugar pass by default (receptors here carry no modelled glycans, and the
# sugar pass can fail on atoms whose element it cannot resolve). Set DOCKFLOW_KEEP_SUGARS=1 to re-enable it.
NOSUGAR_ARG="--no-sugars"; [ -n "${DOCKFLOW_KEEP_SUGARS:-}" ] && NOSUGAR_ARG=""
SALT_CONC="${DOCKFLOW_SALT_CONC:-0.15}"   # bulk salt molarity (M); shared by both paths
export DOCKFLOW_SALT_CONC="$SALT_CONC"    # so pack_membrane uses the same value

if ! mkdir -p "$output_folder"; then
    echo "Could not create folder $output_folder !!!"
    exit 2
fi

set -euo pipefail

pdbfile_abs="$(realpath "$pdbfile")"
cd "$output_folder"

BOX_DIMS_DATA=""        # only set in the membrane branch; initialised so the soluble path is safe under set -u
N_PAIRS="0"             # soluble bulk-salt ion pairs (computed from the water count below)
DISULFIDE_LEAP=""       # path to the CPPTRAJ-detected disulfide LEaP fragment, if any
# Membrane mode: the SIF orchestration sets it via the /mnt/receptor_files/MEMBRANE marker; a bare-metal
# run (no writable /mnt) sets DOCKFLOW_MEMBRANE=1 in the environment instead.
MEMBRANE="${DOCKFLOW_MEMBRANE:-0}"
if test -f "/mnt/receptor_files/MEMBRANE"; then MEMBRANE=1; fi

# ─── Common prep: clean -> PROPKA -> Reduce/CPPTRAJ -> disulfide gate (soluble AND membrane) ─────────────
"$DOCKFLOW_ROOT/ctmd_module/clean_pdb.py" "$pdbfile_abs" clean.pdb

# Membrane: orient the protein so its transmembrane span aligns with the bilayer normal BEFORE protonation.
# MEMEMBED scores only standard residue names, so it must run on the clean PDB (not the Amber-renamed
# CYX/HIE prep output, on which it silently returns zero energy and writes nothing). The oriented
# coordinates carry through PROPKA/Reduce/prepareforleap unchanged, and pack_membrane is then run with
# --preoriented so packmol does not re-run (and fail) MEMEMBED on the prepared protein.
prep_input="clean.pdb"
if [ "$MEMBRANE" -eq 1 ]; then
    MEMEMBED="$(python3 -c 'import os,packmol_memgen;print(os.path.join(os.path.dirname(packmol_memgen.__file__),"lib","memembed","bin","memembed"))' 2>/dev/null || true)"
    if [ -n "$MEMEMBED" ] && [ -x "$MEMEMBED" ]; then
        "$MEMEMBED" -o clean_oriented.pdb clean.pdb > memembed.log 2>&1 || true
    fi
    if [ -s clean_oriented.pdb ]; then
        prep_input="clean_oriented.pdb"
        echo "Membrane: oriented clean protein with MEMEMBED -> clean_oriented.pdb"
    else
        echo "[!] MEMEMBED orientation failed (see memembed.log) — continuing with unoriented clean.pdb"
    fi
fi

# Optional residue-range truncation (generic): excise disordered loop(s) that inflate the solvation/membrane
# box and cap the new chain ends with ACE/NME (built from the flanking residues' own backbone, no linker
# modelling). DOCKFLOW_TRUNCATE_RANGES is a space/comma list of START-END ranges in the input PDB numbering;
# DOCKFLOW_TRUNCATE_CHAIN optionally restricts to one chain. Applied AFTER MEMEMBED orientation so the
# orientation still sees the full, standard-named protein; the prep/pack/tleap steps below then run on the
# smaller capped structure. Disulfide detection runs downstream on this truncated PDB, so its LEaP bond
# indices stay correct.
if [ -n "${DOCKFLOW_TRUNCATE_RANGES:-}" ]; then
    trunc_out="${prep_input%.pdb}_trunc.pdb"
    chain_arg=""; [ -n "${DOCKFLOW_TRUNCATE_CHAIN:-}" ] && chain_arg="--chain ${DOCKFLOW_TRUNCATE_CHAIN}"
    python3 "$DOCKFLOW_ROOT/protein_prep_module/truncate_and_cap.py" \
        "$(realpath "$prep_input")" "$trunc_out" --remove ${DOCKFLOW_TRUNCATE_RANGES} $chain_arg
    prep_input="$trunc_out"
    echo "Truncated ranges [${DOCKFLOW_TRUNCATE_RANGES}] with ACE/NME caps -> $prep_input"
fi

protect_args=""; for r in $PROTECT_RES; do protect_args="$protect_args --protect $r"; done
req_args=""; for d in $REQUIRED_DISULFIDES; do req_args="$req_args --required-disulfide $d"; done

# Run the two prep stages without set -e so an ineligible model (a normal skip) and a real failure can be
# told apart from the manifest/outputs rather than aborting the run.
set +e
python3 "$RECEPTOR_PREP/prepare_amber_protonation.py" "$(realpath "$prep_input")" \
    --out-dir propka_prep --ph "$PH" --pka-cutoff "$PKA_CUTOFF" --cysteine-mode keep $protect_args
states_pdb="$(find propka_prep -name '*_amber_states.pdb' 2>/dev/null | head -1 || true)"
if [ -n "$states_pdb" ]; then
    python3 "$RECEPTOR_PREP/prepare_for_tleap.py" "$states_pdb" \
        --out-dir tleap_prep --reduce --reduce-exe "$REDUCE_EXE" \
        --detect-disulfides --disulfide-sg-cutoff "$SG_CUTOFF" $req_args $NOSUGAR_ARG \
        --leap-unit rec --cpptraj-exe "$CPPTRAJ_EXE"
fi
set -e

manifest="$(find tleap_prep -name '*_prepareforleap_manifest.json' 2>/dev/null | head -1 || true)"
prepped_pdb="$(find tleap_prep -name '*_prepareforleap.pdb' 2>/dev/null | head -1 || true)"
DISULFIDE_LEAP="$(find tleap_prep -name '*_tleap_modifications.leap' 2>/dev/null | head -1 || true)"
if [ -z "$manifest" ] || [ -z "$prepped_pdb" ]; then
    echo "Receptor prep failed (no prepareforleap output). See propka_prep/ and tleap_prep/ logs."
    exit 3
fi

# Eligibility gate: a model missing a required disulfide is skipped cleanly (not a batch failure). With no
# required disulfides configured, every model is eligible (manifest -> True).
eligible="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("eligible_for_membrane_preparation"))' "$manifest")"
if [ "$eligible" != "True" ]; then
    echo "Model is INELIGIBLE (missing required disulfide) — skipping $pdbfile_abs"
    touch SKIPPED_INELIGIBLE
    exit 0
fi

cp "$prepped_pdb" protein.pdb

# ─── Box build: membrane adds the bilayer; soluble lets tleap add the water box below ─────────────────────
if [ "$MEMBRANE" -eq 1 ]
then
    echo "# Membrane setup for $pdbfile_abs"
    # pack_membrane neutralizes and adds DOCKFLOW_SALT_CONC (0.15 M) bulk salt via packmol-memgen --saltcon
    "$DOCKFLOW_ROOT/ctmd_module/md_submodule/pack_membrane" protein.pdb
    mv protein.pdb protein_vacuum.pdb
    mv protein_packed.pdb protein.pdb
    # packmol-memgen writes no CRYST1 record without --parametrize; the periodic box edges are the
    # x_len/y_len/z_len it reports under "Box information" in packmol-memgen.log.
    BOX_DIMS_DATA="$(awk '/x_len/{x=$NF} /y_len/{y=$NF} /z_len/{z=$NF} END{if(x!=""&&y!=""&&z!="") printf "%.3f %.3f %.3f", x, y, z}' packmol-memgen.log 2>/dev/null || true)"
    echo "Box dims: $BOX_DIMS_DATA"
else
    # Soluble: tleap addions is count-based, so to hit DOCKFLOW_SALT_CONC we solvate once to count waters,
    # then compute the NaCl pair count (n_pairs = conc/55.5 * n_water) injected into setup.leap below.
    printf 'source leaprc.protein.ff14SB\nsource leaprc.water.tip3p\nrec = loadpdb protein.pdb\ns = combine {rec}\nsolvateoct s TIP3PBOX 12.0\nsavepdb s watercount.pdb\nquit\n' > count_waters.leap
    tleap -f count_waters.leap > count_waters.log 2>&1 || true
    nwat="$(awk '$1=="ATOM" && $3=="O" && $4=="WAT"{n++} END{print n+0}' watercount.pdb 2>/dev/null || echo 0)"
    N_PAIRS="$(python3 -c "import sys; print(int(round(float('$SALT_CONC')/55.5*int(sys.argv[1]))))" "$nwat" 2>/dev/null || echo 0)"
    echo "Soluble: $nwat waters -> $N_PAIRS NaCl pair(s) for ${SALT_CONC} M bulk salt"
fi

# ─── Write the tleap setup script ────────────────────────────────────────────────
if test -f "setup.leap"
then
    echo "setup.leap found. Using this script"
else
    if [ "$MEMBRANE" -eq 1 ]; then
        template="$DOCKFLOW_ROOT/ctmd_module/md_submodule/membrane_setup.leap"
    elif test -f "/mnt/MD_setup.leap"; then
        template="/mnt/MD_setup.leap"
    else
        template="$DOCKFLOW_ROOT/ctmd_module/md_submodule/setup.leap"
    fi
    sed "s/PDBFILE/protein.pdb/" "$template" > setup.leap

    # insert the CPPTRAJ-detected disulfide bonds at the DISULFIDE_BONDS marker (both soluble and membrane)
    if [ -n "$DISULFIDE_LEAP" ] && [ -f "$DISULFIDE_LEAP" ]; then
        awk -v f="$DISULFIDE_LEAP" '/DISULFIDE_BONDS/{while((getline l < f) > 0) print l; next} {print}' \
            setup.leap > setup.leap.tmp && mv setup.leap.tmp setup.leap
    else
        grep -v "DISULFIDE_BONDS" setup.leap > setup.leap.tmp && mv setup.leap.tmp setup.leap
    fi

    # membrane: substitute the periodic box edges; soluble: substitute the bulk-salt ion-pair count
    if [ "$MEMBRANE" -eq 1 ]; then
        sed "s/BOX_DIMS/$BOX_DIMS_DATA/g" setup.leap > setup.leap.tmp && mv setup.leap.tmp setup.leap
    else
        sed "s/SALT_PAIRS/$N_PAIRS/g" setup.leap > setup.leap.tmp && mv setup.leap.tmp setup.leap
    fi
    echo "[NOTE]: Setup LEAP file written to $output_folder/setup.leap. Edit it and re-run to change the setup."
fi

# tleap may fail (e.g. bad PDB) — continue anyway so we can inspect the outputs/logs
tleap -f setup.leap > leap_output.log 2>&1 || true

if test -f "system.inpcrd" && test -f "system.prmtop"
then
    echo "Setup completed successfully."
    cp -r "$DOCKFLOW_ROOT/ctmd_module/md_submodule/amber_MD_template" .
    # start_MD.sh is the in-SIF launcher (singularity exec of the AmberTools image). It is only written
    # when the SIF's system.conf is present; a bare-metal run uses its own MD runner and doesn't need it,
    # so a missing system.conf is NOT a failure (the topology already built successfully above).
    ambertools_img=""
    if [ -f /usr/software/common/system.conf ]; then
        ambertools_img=$(grep AMBERTOOLS /usr/software/common/system.conf | cut -d"=" -f2 || true)
    fi
    if [ -n "$ambertools_img" ]; then
        # single printf so the shebang is not clobbered (the old `echo >` sequence overwrote line 1)
        printf '#!/bin/bash\nmodule load singularity\nsingularity exec --nv --bind .:/scratch %s /scratch/amber_MD_template/singularity_passer.sh ./ 0\ntouch MD_COMPLETE\n' \
            "$ambertools_img" > start_MD.sh
        chmod +x start_MD.sh
    else
        echo "[note] No SIF system.conf found — skipping start_MD.sh (bare-metal run uses its own MD runner)."
    fi
    exit 0
fi

echo "Setup failed! See $output_folder/leap.log or $output_folder/leap_output.log"
exit 3
