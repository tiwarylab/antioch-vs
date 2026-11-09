# DUD-E Dock6 Examples

This folder contains two small Dock6 redocking examples based on DUD-E kinase
targets:

- `cdk2/`: CDK2, PDB `1h00`, ligand residue `FAP`.
- `abl1/`: ABL1 kinase, PDB `2hzi`, ligand residue `JIN`.

These are intended as ANTIOCH-VS workflow tests, not full DUD-E enrichment
benchmarks. They show how to fetch a receptor/ligand from a PDB ID, prepare the
receptor for Dock6 with Chimera, and run Dock6 redocking with example settings.

The `prepared/` folders include compact prepared files so the Dock6 step can be
tested without repeating receptor preparation:

- `*_receptor_prepared.pdb`: receptor from `fetch_prepare_receptor.py`.
- `*_reference_ligand.mol2`: ligand without hydrogens used to select spheres.
- `*_ligand.mol2`: ligand used for redocking and RMSD reference.
- `*_surf_d5.dms`: DMS surface from Chimera.
- `*_wH_wcH.mol2`: charged receptor Mol2 from Chimera.
- `*_flex_scored.mol2`: example Dock6 redocking output.

## Recreate The Prepared Inputs

From the ANTIOCH-VS repository root:

```bash
python docking_module/receptor_prep/fetch_prepare_receptor.py \
  --pdb-id 1h00 \
  --ligand-resname FAP \
  --docking-software dock6 \
  --output-dir test_files/docking_module_test/dock6_dude_examples/cdk2/prepared
```

```bash
python docking_module/receptor_prep/fetch_prepare_receptor.py \
  --pdb-id 2hzi \
  --ligand-resname JIN \
  --add-missing-residues \
  --docking-software dock6 \
  --receptor-chains A \
  --ligand-chain A \
  --output-dir test_files/docking_module_test/dock6_dude_examples/abl1/prepared
```

Then run Chimera receptor preparation for the target:

```bash
python docking_module/dock6_pipeline/chimera_pipeline.py \
  --config test_files/docking_module_test/dock6_dude_examples/cdk2/example_chimera_config.yaml
```

or:

```bash
python docking_module/dock6_pipeline/chimera_pipeline.py \
  --config test_files/docking_module_test/dock6_dude_examples/abl1/example_chimera_config.yaml
```

Update each Chimera config for your local Chimera execution mode. On an HPC
cluster using a Singularity or Apptainer image, you may need to load the module
first:

```bash
module load singularity
```

## Run Dock6

Run from the ANTIOCH-VS repository root:

```bash
python docking_module/dock6_pipeline/dock6_pipeline.py \
  --config test_files/docking_module_test/dock6_dude_examples/cdk2/example_dock6_config.yaml \
  --overwrite
```

```bash
python docking_module/dock6_pipeline/dock6_pipeline.py \
  --config test_files/docking_module_test/dock6_dude_examples/abl1/example_dock6_config.yaml \
  --overwrite
```

The configs use relative paths and write outputs under each target's
`dock6_run/` directory.

Set `environment.dock6_dir` for a native Dock6 install. To use a Singularity
image, uncomment `environment.sif_image`, set `dock6_dir` to the Dock6
installation path inside the image, and add any required host bind paths under
`environment.singularity_binds`.

## CDK2 vs ABL1

ABL1 is the easier redocking example here. In the included example output, the
top scored ABL1 pose has a Dock6 heavy-atom RMSD of `0.368 A` relative to the
input crystal ligand.

CDK2 is more sensitive to ligand sampling. The CDK2 config intentionally uses
expanded sampling:

```yaml
parameters:
  max_orientations: 5000
  extra_flex_args:
    num_preclustered_conformers: 50
    num_final_scored_poses: 50
    write_orientations: yes
    simplex_final_min: yes
    simplex_max_cycles: 100
    simplex_anchor_max_iterations: 1000
    simplex_grow_max_iterations: 1000
```

In the included CDK2 output, the best-scored pose has `HA_RMSDs 3.139 A`, while
another saved pose reaches `2.413 A`. This mirrors the DUD-E paper's point that
some redocking cases depend strongly on ligand conformer/orientation sampling:
the correct chemistry alone does not guarantee the correct pose is sampled and
ranked first.

Dock6 reports RMSD in the scored Mol2 headers when `calculate_rmsd yes` and
`use_rmsd_reference_mol yes` are enabled. In this pipeline, the RMSD reference is
the `ligand_file`, symlinked locally as `ligands.mol2`.

Quickly inspect RMSDs and scores:

```bash
grep "HA_RMSDs\|Grid_Score" cdk2/prepared/1h00_receptor_prepared_wH_wcH_flex_scored.mol2
grep "HA_RMSDs\|Grid_Score" abl1/prepared/2hzi_receptor_prepared_wH_wcH_flex_scored.mol2
```
