#!/bin/bash
# ============================================================================
# prepare_ctmd.sh — Automate CTMD AMBER MD setup from raw inputs
#
# Usage:
#   cd <ctmd_project>
#   bash /path/to/antioch-vs/ctmd_module/prepare_ctmd.sh <receptor_pdb> [add_ions]
#
# Arguments:
#   receptor_pdb — Path to any prepared receptor PDB file.
#                  If it is not already ./rec_prepped.pdb, it will be
#                  copied to that name in the current directory.
#   add_ions     — Optional. "True" to add Na+/Cl- ions after solvation
#                  Default: False (no ions added beyond what setup.leap specifies)
#
# What this script does:
#   1. Computes net charge for each .mol2 file by summing atomic charges
#   2. Generates hit_data_mol2.spc (ligand list + charges)
#   3. Copies CTMD pipeline scripts (setup.leap, _setup.sh, etc.)
#   4. Patches scripts for mawk compatibility and missing 'mul' binary
#   5. Computes protein net charge via tleap + parmed
#   6. Edits setup.leap to neutralize the system with correct ion count
#   7. Creates per-ligand folders with mol2 copies and Qnet files
#   8. Runs antechamber → parmchk2 → tleap for each ligand (up to 12 parallel)
#
# Output:
#   Each ligand folder (named after its input MOL2 file) will contain:
#     - system.pdb        (solvated complex coordinates)
#     - system.prmtop     (AMBER topology)
#     - system.inpcrd     (AMBER coordinates)
#     - protlig_vac.*     (gas-phase files for MMGBSA)
#
# Requirements:
#   - micromamba chem_env activated (antechamber, parmchk2, tleap, sqm, parmed)
#   - python3 available
#   - CTMD module scripts in antioch-vs/ctmd_module/
# ============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Ensure chem_env is active (required for RDKit + AMBER tools)
# This handles cases where the script is called from a subshell that lost env.
# ---------------------------------------------------------------------------
if ! command -v antechamber &>/dev/null; then
    source /scratch/software/micromamba_root/sourceme.sh 2>/dev/null || {
        export MAMBA_ROOT_PREFIX=/scratch/software/micromamba_root
        export CONDA_PKGS_DIRS=/scratch/software/micromamba_cache/pkgs
        export XDG_CACHE_HOME=/scratch/software/micromamba_cache/xdg
        eval "$(/scratch/software/micromamba_install/bin/micromamba shell hook -s bash)"
    }
    micromamba activate chem_env 2>/dev/null
fi

# ---------------------------------------------------------------------------
# Pre-flight: Verify required AMBER tools are available
# ---------------------------------------------------------------------------
echo "--- Checking dependencies ---"

MISSING_TOOLS=()
for tool in antechamber parmchk2 tleap sqm python3; do
    if ! command -v "$tool" &>/dev/null; then
        MISSING_TOOLS+=("$tool")
    fi
done

# Check parmed Python module separately
if ! python3 -c "import parmed" &>/dev/null; then
    MISSING_TOOLS+=("parmed (Python module)")
fi

if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
    echo ""
    echo "ERROR: Required tools not found in PATH:"
    for t in "${MISSING_TOOLS[@]}"; do
        echo "  - $t"
    done
    echo ""
    echo "Make sure the AMBER chem_env is activated before running this script:"
    echo "  source /scratch/software/micromamba_root/sourceme.sh"
    echo ""
    exit 1
fi

echo "  antechamber : $(command -v antechamber)"
echo "  parmchk2    : $(command -v parmchk2)"
echo "  tleap       : $(command -v tleap)"
echo "  sqm         : $(command -v sqm)"
echo "  python3     : $(command -v python3)"
echo "  parmed      : $(python3 -c 'import parmed; print(parmed.__file__)')"
echo ""

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
if [ $# -lt 1 ]; then
    echo "Usage: bash $0 <receptor_pdb> [add_ions]"
    echo ""
    echo "  receptor_pdb — Path to prepared receptor PDB file"
    echo "  add_ions     — 'True' to add ions (default: False)"
    exit 1
fi

RECEPTOR_PDB="$1"
ADD_IONS="${2:-False}"
CTMD_SRC="$(cd "$(dirname "$0")" && pwd)"
INPUT_DIR="$(pwd)"
LOG_FILE="$INPUT_DIR/ctmd_setup.log"

# ---------------------------------------------------------------------------
# Resolve receptor PDB → always work with ./rec_prepped.pdb
# ---------------------------------------------------------------------------
echo "=== CTMD Preparation Script ==="
echo "Receptor input    : $RECEPTOR_PDB"
echo "Add ions          : $ADD_IONS"
echo "CTMD source       : $CTMD_SRC"
echo "Working directory : $INPUT_DIR"
echo ""

if [ ! -f "$RECEPTOR_PDB" ]; then
    echo "ERROR: PDB file not found: $RECEPTOR_PDB"
    exit 1
fi

# If the input is not already rec_prepped.pdb in this directory, copy it.
TARGET="$INPUT_DIR/rec_prepped.pdb"
if [ "$(realpath "$RECEPTOR_PDB")" != "$(realpath "$TARGET")" ]; then
    echo "Copying $RECEPTOR_PDB → $TARGET"
    cp "$RECEPTOR_PDB" "$TARGET"
fi

shopt -s nullglob
MOL2_FILES=("$INPUT_DIR"/molecules/*.mol2)
shopt -u nullglob
MOL_COUNT=${#MOL2_FILES[@]}
if [ "$MOL_COUNT" -eq 0 ]; then
    echo "ERROR: No .mol2 files found in $INPUT_DIR/molecules/"
    exit 1
fi

echo "Found $MOL_COUNT ligand mol2 files."
echo ""

# ---------------------------------------------------------------------------
# Verify CTMD source has required files
# ---------------------------------------------------------------------------
if [ ! -d "$CTMD_SRC" ]; then
    echo "ERROR: CTMD source directory not found: $CTMD_SRC"
    exit 1
fi

MISSING_CTMDS=()
for f in setup.leap _setup.sh 1_folder_structure.sh 2_setup.sh; do
    if [ ! -f "$CTMD_SRC/$f" ]; then
        MISSING_CTMDS+=("$f")
    fi
done

if [ ${#MISSING_CTMDS[@]} -gt 0 ]; then
    echo "ERROR: Missing CTMD source files in $CTMD_SRC:"
    for f in "${MISSING_CTMDS[@]}"; do
        echo "  - $f"
    done
    exit 1
fi

# ---------------------------------------------------------------------------
# Step 1: Compute net charges and generate hit_data_mol2.spc
# ---------------------------------------------------------------------------
echo "--- Step 1: Generating hit_data_mol2.spc ---"

> "$INPUT_DIR/hit_data_mol2.spc"

for molfile in "$INPUT_DIR"/molecules/*.mol2; do
    name=$(basename "$molfile")

    # Compute net charge from MOL2 file using a three-tier approach:
    # 1. RDKit formal charge perception (most reliable — reads atom types correctly)
    # 2. UNITY formal charges (written by obabel from SDF CHG records)
    # 3. Gasteiger partial-charge sum (fallback, often rounds to 0)
    #
    # Using the correct charge is critical: antechamber SQM crashes with
    # "odd number of electrons" when -nc is wrong.

    # Tier 1: RDKit formal charge
    rd_charge=$(python3 -c "
import sys
try:
    from rdkit import Chem
    mol = Chem.MolFromMol2File(sys.argv[1], removeHs=False)
    if mol is None:
        mol = Chem.MolFromMol2File(sys.argv[1], removeHs=True)
    if mol:
        print(sum(a.GetFormalCharge() for a in mol.GetAtoms()))
    else:
        print('RDFAIL')
except Exception:
    print('RDFAIL')
" "$molfile" 2>/dev/null)

    if [ "$rd_charge" != "RDFAIL" ] && [ -n "$rd_charge" ]; then
        net_charge=$rd_charge
    else
        # Tier 2: UNITY formal charges (from SDF CHG records via obabel)
        unity_charge=$(awk '
            /@<TRIPOS>UNITY_ATOM_ATTR/ { flag=1; next }
            /^@<TRIPOS>/              { flag=0; next }
            flag && /^charge/         { sum += $2 }
            END                       { print sum+0 }
        ' "$molfile")

        if [ "$unity_charge" != "0" ]; then
            net_charge=$unity_charge
        else
            # Tier 3: Fallback — sum Gasteiger partial charges (column 10) and round
            net_charge=$(awk '
                /^@<TRIPOS>ATOM/ { flag=1; next }
                /^@<TRIPOS>/ && flag { flag=0 }
                flag && NF>=10     { sum += $10 }
                END { printf "%d", (sum >= 0) ? int(sum + 0.5) : -int(-sum + 0.5) }
            ' "$molfile")
        fi
    fi

    echo "$name True $net_charge" >> "$INPUT_DIR/hit_data_mol2.spc"
    echo "  $name → net charge = $net_charge"
done

echo "  → hit_data_mol2.spc written."
echo ""

# ---------------------------------------------------------------------------
# Step 2: Copy CTMD scripts and patch for compatibility
# ---------------------------------------------------------------------------
echo "--- Step 2: Copying and patching CTMD scripts ---"

cp "$CTMD_SRC/setup.leap"       "$INPUT_DIR/setup.leap"
cp "$CTMD_SRC/_setup.sh"        "$INPUT_DIR/_setup.sh"
cp "$CTMD_SRC/1_folder_structure.sh" "$INPUT_DIR/1_folder_structure.sh"
cp "$CTMD_SRC/2_setup.sh"       "$INPUT_DIR/2_setup.sh"

# Patch 1: mawk does not support -e flag → remove it
for script in "$INPUT_DIR/_setup.sh" "$INPUT_DIR/1_folder_structure.sh" "$INPUT_DIR/2_setup.sh"; do
    sed -i "s|awk -e '|awk '|g" "$script"
done

# Patch 2: 'mul' binary not available → use python3 for ion count math
sed -i 's|^nat=`mul .*|nat=$(python3 -c "print(int(float(\x27$vol\x27) * 9.033e-05))") #9.033e-05 atoms/A^3 is target conc|' "$INPUT_DIR/_setup.sh"

# Patch 3: Speed up bcc charge equalization with -pl 10
sed -i 's|-c bcc |-c bcc -pl 10 |g' "$INPUT_DIR/_setup.sh"

echo "  Scripts copied and patched."
echo ""

# ---------------------------------------------------------------------------
# Step 2b: Compute protein net charge and patch setup.leap for neutralization
# ---------------------------------------------------------------------------
echo "--- Step 2b: Computing protein charge ---"

# Use tleap to generate a temporary prmtop, then parmed (Python) to read charge.
# This requires chem_env (tleap + parmed). Both are available there.
TMP_LEAP=$(mktemp /tmp/charge_test_XXXXXX.leap)
TMP_PRMTOP=$(mktemp /tmp/charge_test_XXXXXX.prmtop)
TMP_INPCRD=$(mktemp /tmp/charge_test_XXXXXX.inpcrd)

cat > "$TMP_LEAP" << LEAPEOF
source leaprc.protein.ff14SB
rec = loadpdb $INPUT_DIR/rec_prepped.pdb
saveamberparm rec $TMP_PRMTOP $TMP_INPCRD
quit
LEAPEOF

tleap -s -f "$TMP_LEAP" > /dev/null 2>&1
if [ $? -ne 0 ] || [ ! -f "$TMP_PRMTOP" ]; then
    echo ""
    echo "ERROR: tleap failed to generate prmtop for charge calculation."
    echo "  Leap input : $TMP_LEAP"
    echo "  Expected   : $TMP_PRMTOP"
    echo ""
    echo "Common causes:"
    echo "  - rec_prepped.pdb has bad residue names or missing atoms"
    echo "  - Force field (leaprc.protein.ff14SB) not found in chem_env"
    echo "  - tleap crashed (check for segfaults with valgrind if needed)"
    echo ""
    rm -f "$TMP_LEAP" "$TMP_PRMTOP" "$TMP_INPCRD"
    exit 1
fi

# Read net charge via parmed (sum of per-atom charges), round to nearest int
PROT_CHARGE=$(python3 -c "
import parmed as pmd, sys
try:
    s = pmd.load_file('$TMP_PRMTOP')
    total = sum(a.charge for a in s.atoms)
    print(round(total))
except Exception as e:
    print(f'ERROR: {e}', file=sys.stderr)
    sys.exit(1)
") || {
    echo ""
    echo "ERROR: Failed to read protein charge from prmtop via parmed."
    echo "  File: $TMP_PRMTOP"
    echo ""
    rm -f "$TMP_LEAP" "$TMP_PRMTOP" "$TMP_INPCRD"
    exit 1
}

rm -f "$TMP_LEAP" "$TMP_PRMTOP" "$TMP_INPCRD"

echo "  Protein net charge: $PROT_CHARGE"

# Determine which ion and how many to neutralize
if [ "$PROT_CHARGE" -gt 0 ]; then
    ION_TYPE="Cl-"
    ION_COUNT=$PROT_CHARGE
elif [ "$PROT_CHARGE" -lt 0 ]; then
    ION_TYPE="Na+"
    ION_COUNT=$(( -PROT_CHARGE ))
else
    ION_TYPE=""
    ION_COUNT=0
fi

if [ "$ION_COUNT" -gt 0 ]; then
    echo "  Neutralizing with $ION_COUNT x $ION_TYPE"

    # Patch setup.leap: uncomment the matching addions line and set count,
    # comment out the other ion line.
    if [ "$PROT_CHARGE" -gt 0 ]; then
        # Positive protein → needs Cl- counterions
        sed -i "s|^ # addions solvcomplex Cl- .*| addions solvcomplex Cl- $ION_COUNT|" "$INPUT_DIR/setup.leap"
        sed -i "s|^ # addions solvcomplex Na+ .*| # addions solvcomplex Na+ 0|" "$INPUT_DIR/setup.leap"
    else
        # Negative protein → needs Na+ counterions
        sed -i "s|^ # addions solvcomplex Cl- .*| # addions solvcomplex Cl- 0|" "$INPUT_DIR/setup.leap"
        sed -i "s|^ # addions solvcomplex Na+ .*| addions solvcomplex Na+ $ION_COUNT|" "$INPUT_DIR/setup.leap"
    fi
else
    echo "  Protein is neutral — no ions needed."
fi

echo ""

# ---------------------------------------------------------------------------
# Step 3: Create directory structure
# ---------------------------------------------------------------------------
echo "--- Step 3: Creating directory structure ---"

# extracted_poses/ — where mol2 files live for the pipeline
mkdir -p "$INPUT_DIR/extracted_poses"
cp "$INPUT_DIR"/molecules/*.mol2 "$INPUT_DIR/extracted_poses/"
cp "$INPUT_DIR/hit_data_mol2.spc" "$INPUT_DIR/extracted_poses/"

# MD_components/ — where the protein lives
mkdir -p "$INPUT_DIR/MD_components"
cp "$INPUT_DIR/rec_prepped.pdb" "$INPUT_DIR/MD_components/rec_prepped.pdb"

echo "  extracted_poses/ : $(ls "$INPUT_DIR/extracted_poses/"*.mol2 | wc -l) mol2 files"
echo "  MD_components/   : rec_prepped.pdb"
echo ""

# ---------------------------------------------------------------------------
# Step 4: Run 1_folder_structure.sh (creates per-ligand folders)
# ---------------------------------------------------------------------------
echo "--- Step 4: Creating per-ligand folders ---"

cd "$INPUT_DIR"
if ! echo "y" | bash 1_folder_structure.sh; then
    echo ""
    echo "ERROR: 1_folder_structure.sh failed."
    echo "  Check that hit_data_mol2.spc and extracted_poses/ are correct."
    exit 1
fi

FOLDER_COUNT=$MOL_COUNT
echo "  Created $FOLDER_COUNT ligand folders."
echo ""

# ---------------------------------------------------------------------------
# Step 5: Run antechamber → parmchk2 → tleap per ligand (parallel, 1 CPU each)
# ---------------------------------------------------------------------------
echo "--- Step 5: Running AMBER setup for all ligands ---"
echo "  Max parallel jobs : 12"
echo "  CPUs per job      : 1"
echo "  Log               : $LOG_FILE"
echo ""

cd "$INPUT_DIR"

# Read ligand list from hit_data_mol2.spc (column 1 = mol2 filename)
LIGANDS=()
while IFS=' ' read -r fname _rest; do
    LIGANDS+=("$fname")
done < hit_data_mol2.spc

# Function to process a single ligand inside its folder with 1 CPU limit
process_ligand() {
    local molfile="$1"
    local folname="${molfile%.mol2}"
    local logfile="$INPUT_DIR/${folname}_setup.log"

    echo "[$(date '+%H:%M:%S')] Starting $molfile (PID $$)" >> "$LOG_FILE"

    if [ ! -d "$folname" ]; then
        echo "[$(date '+%H:%M:%S')] ERROR: folder $folname not found" >> "$logfile"
        return 1
    fi

    cd "$folname"

    # AMBER setup is a single-process job. Scheduler CPU affinity, when present,
    # is preserved instead of pinning every parallel job to CPU 0.
    if OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 bash ../_setup.sh ligand_raw.mol2 "$(cat Qnet)" "$ADD_IONS" >> "$logfile" 2>&1; then
        echo "[$(date '+%H:%M:%S')] OK: $molfile" >> "$LOG_FILE"
    else
        # BCC/SQM failed — retry with Gasteiger charges (no SQM needed)
        echo "[$(date '+%H:%M:%S')] BCC failed, retrying with Gasteiger charges..." >> "$logfile"
        # Patch _setup.sh to use -c gas instead of -c bcc
        sed -i 's|-c bcc -pl [0-9]* |-c gas |g' ../_setup.sh
        if OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 bash ../_setup.sh ligand_raw.mol2 "$(cat Qnet)" "$ADD_IONS" >> "$logfile" 2>&1; then
            echo "[$(date '+%H:%M:%S')] OK (Gasteiger fallback): $molfile" >> "$LOG_FILE"
        else
            echo "[$(date '+%H:%M:%S')] FAIL: $molfile (exit $? — see $logfile)" >> "$LOG_FILE"
            return 1
        fi
    fi
}
export -f process_ligand
export INPUT_DIR ADD_IONS LOG_FILE

# Run up to 12 ligands in parallel using background jobs + wait
RUNNING=0
MAX_JOBS=12
FAIL_COUNT=0

for molfile in "${LIGANDS[@]}"; do
    process_ligand "$molfile" &
    RUNNING=$((RUNNING + 1))

    # If we hit the limit, wait for any one job to finish before launching more
    if [ "$RUNNING" -ge "$MAX_JOBS" ]; then
        wait -n
        RUNNING=$((RUNNING - 1))
    fi
done

# Wait for all remaining jobs
wait

echo ""
echo "  All ligand jobs finished. See $LOG_FILE for details."
echo ""

# ---------------------------------------------------------------------------
# Step 6: Summary
# ---------------------------------------------------------------------------
echo ""
echo "=== Setup Complete ==="
echo ""

SUCCESS=0
FAIL=0

for molfile in "${LIGANDS[@]}"; do
    name="${molfile%.mol2}"
    d="$INPUT_DIR/$name"
    if [ -f "$d/system.pdb" ] && [ -f "$d/system.prmtop" ] && [ -f "$d/system.inpcrd" ]; then
        echo "  $name — system.pdb, system.prmtop, system.inpcrd ready"
        SUCCESS=$((SUCCESS + 1))
    else
        echo "  $name — incomplete (check log for errors)"
        FAIL=$((FAIL + 1))
    fi
done

echo ""
echo "  Successful : $SUCCESS / $MOL_COUNT"
echo "  Failed     : $FAIL / $MOL_COUNT"
echo ""

if [ "$FAIL" -gt 0 ]; then
    echo "  Check $LOG_FILE for details on failed ligands."
fi

echo ""
echo "=== Next Steps ==="
echo "Each ligand folder now contains:"
echo "  system.pdb      — Solvated complex (for visualization)"
echo "  system.prmtop   — AMBER topology file"
echo "  system.inpcrd   — AMBER coordinate file"
echo "  protlig_vac.*   — Gas-phase files (for MMGBSA)"
echo ""
echo "To launch MD, use pmemd or sander:"
echo "  pmemd -O -i mdin_file -o mdout -p system.prmtop -c system.inpcrd -r restart.rst7 -x trajectory.mdcrd"
