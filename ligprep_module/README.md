# ligprep_module

`ligprep_module` prepares ligand libraries for downstream docking workflows by converting:

1. **SMILES → UMD → UMF** (internal Universal Molecular Format)
2. **UMF → docking-ready formats** (for example `pdbqt` or `mol2`)

The main entrypoint is:

- `ligprep_pipeline.py`

## What this module does

- Accepts either a SMILES file (`source.type: "smiles"`) or an existing UMF file (`source.type: "umf"`).
- When the source is SMILES, runs:
  - `smiles2umd` to generate UMD
  - `umfgen` to generate UMF
- Runs `umfdump` to export UMF into the target format (`pdbqt`, `mol2`, `sdf`, `smi`, or `umd`).
- Supports:
  - single-folder output or multi-block output (`num_blocks`)
  - one-file-per-molecule or combined file output (`split_files`)

## Required tools

This module orchestrates external UMDFP tools:

- `smiles2umd`
- `umfgen`
- `umfdump`
- `python3` (for `smiles2umd`)

Binary paths are configured in YAML (`umdfp` section). If not absolute, the pipeline tries:

1. `../../UMDFP/bin/<tool>`
2. the system `PATH`

## Configuration

Use YAML configs with `ligprep_pipeline.py`.  
Template and examples:

- `config_template.yaml`
- `config_case1.yaml` — PDBQT, single output folder, split files
- `config_case2.yaml` — PDBQT, 4 blocks, split files
- `config_case3.yaml` — Mol2, 4 blocks, combined file per block

Example input data:

- `input_data/example.smi`

## Usage

From `ligprep_module/`:

```bash
python ligprep_pipeline.py config_case1.yaml
```

Or with your own config:

```bash
python ligprep_pipeline.py /absolute/path/to/config.yaml
```

## Outputs for docking pipeline

Set `processing.target_format` based on the downstream docking workflow:

- `pdbqt` for AutoDock Vina pipelines
- `mol2` for workflows expecting Mol2 ligand libraries

Final prepared ligands are written to `processing.output_dir` (or block-suffixed folders when `num_blocks > 1`).
