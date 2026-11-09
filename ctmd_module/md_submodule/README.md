# MD setup (`md_submodule`)

Builds an AMBER MD system (topology + coordinates) from a receptor PDB, for both the AF2RAVE SPIB MD
and the CTMD MD. Entry point: **`setup_MD_from_PDB.sh <receptor.pdb> [output_folder]`**.

```
input PDB
  └─ COMMON PREP (every receptor):  clean_pdb.py → PROPKA → Reduce/CPPTRAJ prepareforleap
                                    → disulfide detection + required-disulfide gate (skip ineligible)
     ├─ soluble   → tleap setup.leap          (solvateoct TIP3P + neutralize + 0.15 M salt)  + disulfides
     └─ membrane  → pack_membrane (POPC + water + neutralize + 0.15 M salt) → tleap membrane_setup.leap + disulfides
```

The prep is **common to soluble and membrane** — both receptors are PROPKA-protonated, disulfide-checked
(and gated), solvated and neutralized, and get a disulfide-aware topology; they differ only in *how* the
box/waters are built (octahedral water box vs POPC bilayer). Force field everywhere: **ff14SB + GAFF2 +
TIP3P** (+ **Lipid21** for membranes).

## Soluble vs membrane

The mode is chosen by the marker file **`/mnt/receptor_files/MEMBRANE`** (created by `init -m`). No marker
→ the original soluble path (unchanged). Marker present → the membrane path below.

## Common prep (both modes)

Self-contained — `setup_MD_from_PDB.sh` runs the prep and reads the gate itself, so the only orchestration
requirement for membranes is the `MEMBRANE` marker. Per receptor, regardless of mode:

1. **`clean_pdb.py`** — strip hetero/waters/ions, resolve altlocs, fix formatting.
2. **PROPKA protonation** — `docking_module/receptor_prep/prepare_amber_protonation.py` (`--cysteine-mode keep`).
3. **Reduce + CPPTRAJ `prepareforleap`** — `prepare_for_tleap.py` (`--leap-unit rec`, `--detect-disulfides`),
   producing the prepared PDB, the `*_prepareforleap_manifest.json` (eligibility), and the
   `*_tleap_modifications.leap` disulfide `bond` lines.
4. **Eligibility gate** — if `eligible_for_membrane_preparation` is false (a required disulfide did not
   form), the receptor is **skipped cleanly** (`touch SKIPPED_INELIGIBLE`, `exit 0`) — not a batch failure.
   With no required disulfides configured, every model is eligible (detected bonds are still applied).

## Box build (where the two modes diverge)

- **soluble** — `tleap setup.leap`: loads the prepared protein, inserts the detected disulfide bonds,
  `solvateoct TIP3PBOX`, neutralizes (`addions … 0`), **and adds bulk salt to `DOCKFLOW_SALT_CONC`**. Since
  tleap's `addions` is count-based, the driver first solvates once to count waters, computes the NaCl pair
  count (`n_pairs = conc/55.5 × n_waters`), and injects it.
- **membrane** — `pack_membrane` packs a POPC bilayer + water, neutralizes, and adds `DOCKFLOW_SALT_CONC`
  bulk salt (packmol-memgen `--saltcon`; packing only, **no `--parametrize`**); then `tleap
  membrane_setup.leap` builds the topology from the packed PDB, inserts the disulfide bonds, and sets the
  box from the packmol-memgen `CRYST1` record.

Both paths therefore end up **neutralized + `DOCKFLOW_SALT_CONC` (default 0.15 M) bulk salt**.

Both build the topology in tleap (not packmol-memgen's `--parametrize`), which is what lets the CPPTRAJ-detected
disulfides be applied.

## Configuration (environment variables, with defaults)

Target-specific prep is driven by env vars so nothing is hardcoded to one protein. The orchestration layer
(or the user) exports them:

| Variable | Default | Meaning |
|---|---|---|
| `DOCKFLOW_PH` | `7.0` | PROPKA target pH |
| `DOCKFLOW_PKA_CUTOFF` | `2.0` | PROPKA pKa cutoff |
| `DOCKFLOW_PROTECT_RES` | *(empty)* | space-separated residues protected from PROPKA, e.g. `A:48 A:82 A:99` |
| `DOCKFLOW_REQUIRED_DISULFIDES` | *(empty)* | space-separated obligatory SS pairs, e.g. `A:75-A:160 A:153-A:159`. Empty ⇒ **no gate** (every model eligible; detected bonds still applied). |
| `DOCKFLOW_DISULFIDE_SG_CUTOFF` | `2.3` | SG–SG distance cutoff (Å) for disulfide detection |
| `DOCKFLOW_SALT_CONC` | `0.15` | bulk salt molarity (M); shared by both paths (soluble `addions` count calc + membrane `pack_membrane --saltcon`) |
| `DOCKFLOW_REDUCE_EXE` | `reduce` | Reduce executable |
| `DOCKFLOW_CPPTRAJ_EXE` | `cpptraj` | CPPTRAJ executable |

### ADRB2 (β2-adrenergic receptor) values

```bash
export DOCKFLOW_REQUIRED_DISULFIDES="A:75-A:160 A:153-A:159"   # the two native extracellular disulfides
export DOCKFLOW_PROTECT_RES="A:48 A:82 A:99"                    # functionally important residues
export DOCKFLOW_PH=7.0
```
With the two required pairs, the 20 ADRB2 AF2 cluster centers split **16 eligible / 4 skipped**.

## Orchestration / SIF wiring (for the image rebuild)

`setup_MD_from_PDB.sh` reads the `DOCKFLOW_*` variables **from its environment at runtime**. It is invoked
by `proceed`'s `cluster_and_run_MD` loop inside the image — and `proceed` / `.init.py` live in the image's
`/bin/ext`, **not in this repo**. So the exports must be wired in there when the image is (re)built.

> ⚠️ **If the vars are unset, the membrane branch silently runs with empty defaults — no disulfide gate and
> no protected residues — producing a topology without the native disulfides.** For any disulfide-bearing
> target (e.g. ADRB2) they MUST be set.

**Suggested wiring (one place):** have `.init.py` write an env file when `-m/--membrane` is given, and have
`proceed` source it before the MD loop.

In `.init.py`, where it already does `touch /mnt/receptor_files/MEMBRANE`, also write the params it was
given (as CLI flags or from a config):
```python
with open("/mnt/receptor_files/membrane_env.sh", "w") as f:
    f.write(f'export DOCKFLOW_REQUIRED_DISULFIDES="{args.required_disulfides}"\n')
    f.write(f'export DOCKFLOW_PROTECT_RES="{args.protect_res}"\n')
    f.write(f'export DOCKFLOW_PH="{args.ph}"\n')
```
In `proceed`'s `cluster_and_run_MD`, before the `setup_MD_from_PDB.sh` loop:
```bash
if test -f /mnt/receptor_files/membrane_env.sh; then
    source /mnt/receptor_files/membrane_env.sh
fi
```
Nothing else needs to change — `setup_MD_from_PDB.sh` picks the values up from the environment.

**Verify after rebuild** (ADRB2): run with `-m` and the ADRB2 values; confirm an eligible receptor's
`MD_setup/setup.leap` contains the `bond rec.<i>.SG rec.<j>.SG` lines, and that the 4 ineligible centers
each produce a `SKIPPED_INELIGIBLE` marker (16 eligible / 4 skipped).

## Image requirements

The membrane image must provide **propka3, reduce, cpptraj, tleap, and packmol-memgen** (AmberTools).
`pack_membrane` sources `/usr/software/sourceme.sh` and activates `/usr/software/my_env`.

## Files

| File | Role |
|---|---|
| `setup_MD_from_PDB.sh` | entry point; soluble + membrane branches, eligibility gate, bug fixes |
| `setup.leap` | soluble tleap template (TIP3P octahedral box + neutralize + bulk salt); `PDBFILE` / `DISULFIDE_BONDS` / `SALT_PAIRS` substitution points |
| `membrane_setup.leap` | membrane tleap template (ff14SB + Lipid21 + TIP3P); `PDBFILE` / `DISULFIDE_BONDS` / `BOX_DIMS` substitution points |
| `pack_membrane` | packmol-memgen bilayer packing (pack-only) |
| `extract_mdout.sh`, `strip.in`, `framedump.in` | post-MD trajectory extraction (COLVAR prep) |
| `amber_MD_template/` | equilibration/production `.in` files + launcher |

> **Note (membrane post-processing):** `extract_mdout.sh`/`strip.in` currently strip water; for membrane
> trajectories they may also need to strip lipids (POPC) before COLVAR generation. Verify before the SPIB step.
