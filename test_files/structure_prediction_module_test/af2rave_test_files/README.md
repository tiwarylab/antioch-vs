# AF2RAVE Test Files

This example runs AF2RAVE feature selection, filtering, and clustering on PDB
structures produced by a separate prediction runner.

The main module files live in:

```text
structure_prediction_module/af2rave_pipeline/
```

## Setup

Install AF2RAVE outside ANTIOCH-VS. A separate environment is recommended.

```bash
micromamba create -n af2rave_env python=3.10
micromamba activate af2rave_env
pip install git+https://github.com/tiwarylab/af2rave.git
```

Update the predicted-structure directory and output directory:

```yaml
paths:
  structures_dir: "/path/to/predicted_pdbs"
  output_dir: "/path/to/af2rave_test_output"
```

## Run AF2RAVE Analysis

```bash
python structure_prediction_module/af2rave_pipeline/af2rave_pipeline.py \
  --config test_files/structure_prediction_module_test/af2rave_test_files/example_af2rave_config.yaml
```

## Combine Prediction Folders

AF2RAVE can combine predictions from several independent folders before it
starts filtering and clustering. Use a YAML list:

```yaml
paths:
  structures_dir: null
  structure_dirs:
    - "/path/to/localcolabfold_output/runs/rmsa_8_16"
    - "/path/to/localcolabfold_output/runs/rmsa_16_32"
```

or a text file containing one directory path per line:

```yaml
paths:
  structure_dirs_file: "/path/to/prediction_directories.txt"
```

All PDBs are copied to `output_dir/structures/`. If two inputs have the same
filename, ANTIOCH-VS prefixes the collected name and records its origin in
`structure_collection_manifest.csv`.

## Existing PDB Structures

Use the feature-selection config for ColabFold, native AlphaFold 2, normalized
AlphaFold 3, or other existing PDB structures.

Copy or edit:

```text
example_af2rave_feature_selection_config.yaml
```

Important fields:

```yaml
paths:
  structures_dir: "/path/to/existing_pdb_structures"
  structure_dirs: []
  structure_dirs_file: null
  reference_structure: null
  reference_selection: "best_plddt"
  output_dir: "/path/to/af2rave_test_output"

filtering:
  rmsd_selection: "name CA"
  rmsd_cutoff: 6.0
  steric_clash_cutoff: 1.0
  nonbonded_cutoff: null
  peptide_bond_mean: null
  peptide_bond_std: null

analysis:
  write_full_feature_stats: false
```

Set a cutoff to turn on a filter. Leave it as `null` to skip that filter.
For the DDR1 SI-style test, use only the 6 A RMSD filter.
With `write_full_feature_stats: false`, expensive nonbonded and peptide-bond
statistics are skipped unless their filters are enabled.

`reference_structure` controls the AF2RAVE reference used for RMSD alignment,
RMSD filtering, and regular-space clustering. Set it to a specific PDB if you
already know the reference structure:

```yaml
paths:
  reference_structure: "/path/to/pred_1.pdb"
```

If `reference_structure` is `null`, `reference_selection: "best_plddt"` uses
the PDB with the highest residue-averaged pLDDT from the PDB B-factor column.
Each residue contributes once to this average. If pLDDT values are not
available, the pipeline falls back to AF2RAVE's default behavior: the first
loaded PDB is used as the reference. To force that original behavior, use:

```yaml
paths:
  reference_selection: "first"
```

Then run:

```bash
python structure_prediction_module/af2rave_pipeline/af2rave_pipeline.py \
  --config test_files/structure_prediction_module_test/af2rave_test_files/example_af2rave_feature_selection_config.yaml
```

## Expected Outputs

The output folder contains collected structures, filtering statistics, selected
features, cluster-center structures, PCA plots, and simulation boxes when
`simulation_box.enabled` is true.

Important files:

- `feature_stats.csv`
- `feature_stats.png`
- `plddt_scores.csv`
- `filtered_structures.txt`
- `selected_features.csv`
- `feature_cv_histogram.png`
- `cluster_centers.txt`
- `pca_clusters.png`
- `simulation_boxes/boxed_*.pdb`
- `simulation_boxes/index.npy`
- `af2rave_pipeline.log` with the config snapshot and progress messages
