# Receptor Preparation

This folder contains utilities for preparing receptor and ligand starting files
for the docking workflows.

## Fetch and Prepare a PDB Complex

`fetch_prepare_receptor.py` can download a structure by PDB ID from RCSB, clean
and protonate the receptor with PDBFixer/OpenMM, extract a ligand by residue
name, and write ligand PDB/SDF files for downstream docking preparation.
It can also write workflow-specific ligand files for DOCK6 and receptor/ligand
files for Vina when requested.
It also writes a preparation report so users can inspect what PDBFixer changed.
For ligand SDF output, the script uses the RCSB Chemical Component Dictionary
ideal SDF by default to restore ligand bond orders that are not reliably stored
in PDB files.

### Dependencies

Install these in the environment where you run the script:

```bash
micromamba install -n docking_pipeline -c conda-forge pdbfixer openmm biopython rdkit
```

The existing Vina environment already includes `biopython` and `rdkit`, but may
need `pdbfixer` and `openmm` added.
The Vina workflow also expects Meeko's `mk_prepare_receptor.py` and
`mk_prepare_ligand.py` commands to be available on `PATH`.

## Optional AlphaFold Amber Relaxation

`relax_receptor_amber.py` performs the restrained Amber minimization used by
AlphaFold's relaxation workflow. It is useful for resolving local clashes in
predicted receptor structures before preparing them for docking. It restrains
heavy atoms with a harmonic potential; it is not molecular dynamics or a
replacement for structure refinement.

The script requires an AlphaFold or local-ColabFold Python environment that
provides `alphafold.relax` and OpenMM. It is not part of the standard
`docking_pipeline` environment. For CPU jobs, the script automatically sets
`JAX_PLATFORMS=cpu` to avoid CUDA probing on nodes without a GPU.

```bash
COLABFOLD_PYTHON=/path/to/localcolabfold/.pixi/envs/default/bin/python

"$COLABFOLD_PYTHON" docking_module/receptor_prep/relax_receptor_amber.py \
  --input-pdb predicted_receptor.pdb \
  --output-dir relaxed_receptors \
  --max-iterations 2000 \
  --tolerance 2.39 \
  --stiffness 10.0 \
  --max-outer-iterations 3
```

For a stronger local relaxation, increase `--max-iterations`, decrease
`--tolerance`, or decrease `--stiffness`. Inspect the output structure before
using it for docking, especially around a repaired binding site. Each relaxed
PDB receives an accompanying `.amber_relax.json` report.

### Example: Fetch by PDB ID

```bash
python docking_module/receptor_prep/fetch_prepare_receptor.py \
  --pdb-id 1xgj \
  --ligand-resname LIG \
  --output-dir receptor_prep_out
```

### Example: Use a Local Complex PDB

```bash
python docking_module/receptor_prep/fetch_prepare_receptor.py \
  --input-pdb /absolute/path/to/complex.pdb \
  --receptor-chains A B \
  --ligand-resname LIG \
  --ligand-chain A \
  --ligand-resseq 501 \
  --output-dir receptor_prep_out
```

### Outputs

For a prefix such as `1xgj`, the script writes:

- `1xgj_receptor_prepared.pdb`
- `1xgj_complex_prepared.pdb` when `--write-prepared-complex` is used
- `1xgj_LIG_reference_ligand.mol2` when `--docking-software dock6` is used
- `1xgj_LIG_ligand.mol2` when requested, copied from `--prepared-ligand-mol2`, or needed by a workflow
- `1xgj_receptor_prepared.pdbqt` when `--docking-software vina` is used
- `ligands/1xgj_LIG_ligand.pdbqt` when `--docking-software vina` is used
- `1xgj_prep_report.txt`

By default, intermediate files are temporary and are removed after the run:
downloaded complex PDB, all-chain receptor PDB, extracted ligand PDB, ligand
SDF, and fetched CCD template SDF. Use `--keep-all-files` to keep these
debugging/intermediate files in the output folder.

Use the prepared receptor as the receptor input for downstream docking modules.
Use `--write-prepared-complex` if you also want a prepared receptor plus the
original ligand pose for visual inspection.
Use the ligand MOL2/PDBQT workflow outputs as the ligand starting point for
downstream docking.

### Docking Software Mode

Use `--docking-software` to keep the outputs focused on the next docking tool.
The default is `generic`, which writes the receptor and report unless MOL2
output or retained intermediates are explicitly requested.

For DOCK6, the script writes a no-hydrogen reference ligand MOL2 from the
crystal ligand. Use this reference ligand to define the docking site. If you
also provide a carefully prepared ligand MOL2, the script copies it into the
output folder for the actual docking ligand:

```bash
python docking_module/receptor_prep/fetch_prepare_receptor.py \
  --pdb-id 1xgj \
  --receptor-chains B \
  --ligand-resname HTC \
  --ligand-chain B \
  --output-dir receptor_prep_out \
  --docking-software dock6 \
  --prepared-ligand-mol2 /absolute/path/to/HTC_prepared.mol2
```

For Vina, the script converts the prepared receptor to PDBQT with
`mk_prepare_receptor.py`, writes or copies the target ligand MOL2, then converts
the ligand to PDBQT with `mk_prepare_ligand.py`. The preparation report also
records the center of mass of the reference ligand, which can be used as the
Vina box center:

```bash
python docking_module/receptor_prep/fetch_prepare_receptor.py \
  --pdb-id 1iep \
  --receptor-chains A \
  --ligand-resname STI \
  --ligand-chain A \
  --output-dir receptor_prep_out \
  --docking-software vina
```

This public 1IEP example intentionally does not use `--add-missing-residues`.
In this structure, missing-residue rebuilding can produce a receptor PDB that
Meeko/RDKit rejects during receptor PDBQT conversion. For a lightweight Vina
workflow test, use the deposited chain as-is.

The Vina docking config can then use the generated receptor PDBQT directly so
the docking pipeline bypasses receptor conversion:

```yaml
paths:
  receptor_pdbqt: "receptor_prep_out/1iep_receptor_prepared.pdbqt"
  receptor_pdb: ""
  ligand_dir: "receptor_prep_out/ligands"
```

Use the `ligands/` subfolder for `ligand_dir`; the top-level output directory
also contains the receptor PDBQT, which should not be included in Vina batch
ligand input.

If Meeko commands are not on your `PATH`, pass them explicitly:

```bash
--mk-prepare-receptor-command /absolute/path/to/mk_prepare_receptor.py
--mk-prepare-ligand-command /absolute/path/to/mk_prepare_ligand.py
```

If `--vina-receptor-source crystal` is used on a receptor with incomplete
crystal side chains, Meeko may refuse residues that do not match residue
templates. You can allow Meeko to remove those residues with:

```bash
--meeko-allow-bad-res
```

Use this only after checking that the reported incomplete residues are not in
the binding site you care about. For lower-level Meeko options, pass a quoted
string with `--mk-prepare-receptor-args`, for example
`--mk-prepare-receptor-args "--default_altloc A"`.

Use `--receptor-chains` for multichain structures when you only want to save
specific protein chains in the final receptor:

```bash
python docking_module/receptor_prep/fetch_prepare_receptor.py \
  --pdb-id 1xgj \
  --receptor-chains A \
  --ligand-resname LIG \
  --output-dir receptor_prep_out
```

The script still prepares the full receptor first, then saves the selected
chain(s) at the end. This preserves missing-residue records from the original
PDB so PDBFixer can model missing residues before final chain filtering. The
ligand is also extracted from the original complex PDB, so ligand chain
selection remains independent from receptor chain selection.

When `--add-missing-residues` is used, ANTIOCH-VS models internal missing-residue
gaps but skips missing residues at chain termini by default. Terminal tails are
often unresolved because they are flexible and can be poorly constrained for
docking prep. Use `--add-terminal-missing-residues` together with
`--add-missing-residues` only when you intentionally want PDBFixer to rebuild
N- or C-terminal missing segments.

For structures with engineered insertions, fusion proteins, or other unwanted
protein blocks, remove residue ranges after receptor preparation with
`--remove-receptor-residues`. Use `START-END` to remove that residue-number
range from all final receptor chains, or `CHAIN:START-END` to target a specific
chain:

```bash
python docking_module/receptor_prep/fetch_prepare_receptor.py \
  --pdb-id 3ny8 \
  --receptor-chains A \
  --ligand-chain A \
  --ligand-resname JRZ \
  --docking-software dock6 \
  --add-missing-residues \
  --remove-receptor-residues A:340-1161
```

The preparation report records each requested range and how many protein
residues were removed. After removing the ranges, the script reruns PDBFixer on
the trimmed receptor to add missing atoms/hydrogens and re-evaluate terminal
atoms at the new chain breaks, but it does not remodel the removed ranges as
missing residues. Ligand extraction still uses the original complex, so do not
remove the ligand residue with this option.

### Ligand Bond Orders and Optimization

PDB files usually do not preserve ligand bond orders. If RDKit builds an SDF
directly from a ligand PDB, double or aromatic bonds can become single bonds and
extra hydrogens can be added. To avoid that, the script fetches:

```text
https://files.rcsb.org/ligands/download/<RESNAME>_ideal.sdf
```

and transfers those template bond orders onto the extracted ligand coordinates
before writing the final SDF.

For custom ligands, provide your own template:

```bash
python docking_module/receptor_prep/fetch_prepare_receptor.py \
  --input-pdb /absolute/path/to/complex.pdb \
  --ligand-resname LIG \
  --ccd-template-sdf /absolute/path/to/LIG_template.sdf \
  --output-dir receptor_prep_out
```

If you intentionally want RDKit to infer bond orders from the PDB geometry, use
`--no-template-bond-orders`.

By default, the script adds hydrogens and runs a short UFF cleanup while keeping
the ligand heavy atoms fixed. This relaxes added hydrogens without moving the
crystallographic ligand pose away from the binding site. Use
`--no-optimize-ligand` to skip this cleanup.

Use `--write-ligand-mol2` in generic mode if you also want MOL2 output. The
MOL2 file is written from the same corrected ligand molecule used for SDF, so
it keeps the template bond orders and final coordinates. Its charge column uses
RDKit Gasteiger partial charges by default:

```bash
python docking_module/receptor_prep/fetch_prepare_receptor.py \
  --pdb-id 1xgj \
  --ligand-resname HTC \
  --ligand-chain B \
  --output-dir receptor_prep_out \
  --write-ligand-mol2
```

Use `--ligand-mol2-charge-model formal` only if you want formal charges instead
of Gasteiger partial charges in the MOL2 charge column.

For reproducible DOCK6 jobs, use a ligand MOL2 that was prepared with the
desired protonation state and charges. Pass it with `--prepared-ligand-mol2` to
copy it into the receptor-prep output without changing its chemistry:

```bash
python docking_module/receptor_prep/fetch_prepare_receptor.py \
  --pdb-id 1xgj \
  --ligand-resname HTC \
  --ligand-chain B \
  --output-dir receptor_prep_out \
  --docking-software dock6 \
  --prepared-ligand-mol2 /absolute/path/to/HTC_charged.mol2
```

The preparation report records the ligand formal charge, atom counts, template
SDF, and MOL2 charge sum. This receptor-prep tool is intended for quick tests.
For careful docking jobs, inspect the ligand and prepare the desired
protonation/charge state separately. If you have a guide SDF with the desired
ligand chemistry, pass it with `--ccd-template-sdf` so the script uses that
template instead of the default CCD SDF.

### Preparation Report

The report lists:

- missing residues detected by PDBFixer
- whether missing residues were modeled or left out
- missing-residue segments near the selected ligand, based on nearby flanking residues
- nonstandard residues replaced by PDBFixer
- missing heavy atoms and terminal atoms added
- residue-numbering gaps that remain in the final prepared receptor
- residue ranges removed from the receptor and whether post-removal cleanup ran
- whether restrained minimization was requested and applied
- ligand formal charge, atom counts, template source, and optional MOL2 charge sum
- reference ligand center of mass for choosing a Vina box center

By default, the script adds missing atoms and hydrogens but does not model
missing residues or loops. Use `--add-missing-residues` only when you want
PDBFixer to attempt to build internal missing residues. Terminal missing
residue segments are still skipped unless you also pass
`--add-terminal-missing-residues`.

If many missing residues are close to the ligand, the report prints a warning.
That is a sign that the binding site may contain unresolved loops, and PDBFixer
may not be the best tool for modeling them. Consider a dedicated loop/modeling
workflow and manual inspection before docking.

The final residue-numbering continuity check is a quick sanity check for gaps
that remain after preparation. Some numbering gaps are expected for terminal
truncations or author numbering schemes, but internal gaps near the ligand
should be inspected before docking.

### Optional Restrained Receptor Minimization

Use `--minimize-receptor` if you want a quick local relaxation after adding
atoms and hydrogens:

```bash
python docking_module/receptor_prep/fetch_prepare_receptor.py \
  --pdb-id 1xgj \
  --ligand-resname LIG \
  --output-dir receptor_prep_out \
  --minimize-receptor
```

The minimization uses OpenMM with positional restraints on non-hydrogen atoms.
This is meant to reduce local clashes without letting the receptor drift far
from the crystal structure. Keep it short and inspect the output if the binding
site is sensitive.

### Useful Options

- `--keep-water`: keep crystallographic waters during receptor preparation
- `--receptor-chains A B`: save only selected receptor chains after full receptor preparation
- `--remove-receptor-residues A:100-200`: remove protein residue ranges after receptor preparation
- `--docking-software generic|dock6|vina`: choose workflow-specific docking outputs
- `--write-prepared-complex`: also write a prepared receptor plus original ligand-pose PDB
- `--keep-all-files`: keep downloaded/intermediate PDB, SDF, and CCD files
- `--add-missing-residues`: let PDBFixer add internal missing residues
- `--add-terminal-missing-residues`: also let PDBFixer add N-/C-terminal missing residues
- `--ph 7.4`: set the pH used for hydrogen addition
- `--no-optimize-ligand`: skip RDKit hydrogen/local UFF cleanup for SDF output
- `--write-ligand-mol2`: also write the ligand as MOL2
- `--prepared-ligand-mol2`: copy a user-prepared MOL2 instead of generating MOL2 with RDKit
- `--ligand-mol2-charge-model gasteiger`: MOL2 charge model, either `gasteiger` or `formal`
- `--ccd-template-sdf`: provide a custom ligand SDF template for bond orders
- `--no-template-bond-orders`: skip CCD/template bond-order assignment
- `--mk-prepare-ligand-command`: Vina ligand PDBQT conversion command
- `--mk-prepare-receptor-command`: Vina receptor PDBQT conversion command
- `--mk-prepare-receptor-args`: extra quoted arguments for Vina receptor PDBQT conversion
- `--meeko-allow-bad-res`: pass `--allow_bad_res` to Meeko for incomplete crystal residues
- `--report-file`: write the preparation report to a custom path
- `--missing-residue-near-ligand-cutoff`: distance cutoff for pocket-proximity warnings
- `--missing-residue-warning-threshold`: number of nearby missing residues needed to trigger a warning
- `--minimize-receptor`: run short restrained OpenMM minimization
- `--minimize-max-iterations`: cap minimization iterations
- `--minimize-restraint-k`: set heavy-atom restraint strength in kcal/mol/A^2
