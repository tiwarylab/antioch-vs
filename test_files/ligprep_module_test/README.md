# Ligprep Module Test Files

This directory contains small example inputs and YAML configs for smoke-testing
the ligand preparation pipeline (`ligprep_module/ligprep_pipeline.py`).

## Contents

- `input_data/`
  - `example.smi` — Small SMILES file (~20 compounds) for testing SMILES → PDBQT conversion
  - `example.umf` / `example.umfp` — Pre-built UMF file for testing UMF → PDBQT conversion (skips the SMILES step)
- `config_template.yaml` — Full template with all configurable options documented
- `config_case1.yaml` — Minimal config: SMILES → PDBQT, single output folder, one file per molecule

## Setting Up the Python Environment

The pipeline itself depends on a few Python packages:

| Package | Purpose |
|---|---|
| `pyyaml` | YAML config parsing |
| `rdkit` | Used by `smiles2umd` for SMILES → UMD conversion |
| `openff-toolkit` | Partial charge calculation (AM1-BCC) during smiles2umd |

You can create a dedicated micromamba environment with these dependencies:

```bash
micromamba create -n ligprep_env \
  -c conda-forge \
  python=3.10 rdkit openff-toolkit pyyaml -y
```

Then activate it before running the pipeline:

```bash
micromamba activate ligprep_env
```

> **Note:** The `chem_env` environment used in development includes additional
> packages for docking and MD workflows. The command above installs only what
> is needed for ligand preparation. Adjust the Python version as desired.

## External Tools Required (UMDFP)

The pipeline orchestrates three external UMDFP binaries. These must be on `PATH`
or referenced via absolute paths in the config.

First, clone the UMDFP repository alongside antioch-vs:

```bash
git clone https://github.com/venkataadury/UMDFP.git
```

Place it in the same parent directory as antioch-vs so that the sibling path
`../../UMDFP/bin/<tool>` resolves correctly (i.e., both repos share a common
parent):

```
<parent>/
├── antioch-vs/
│   └── ligprep_module/
└── UMDFP/
    └── bin/
        ├── smiles2umd
        ├── umfgen
        └── umfdump
```

| Binary | Purpose | Type |
|---|---|---|
| `smiles2umd` | SMILES → UMD conversion | Python script |
| `umfgen` | UMD → UMF conversion | C++ binary |
| `umfdump` | UMF → pdbqt/mol2/sdf/smi/umd export | C++ binary |

By default, the pipeline resolves these via:

1. Absolute path (if configured in YAML under `umdfp:`)
2. Relative to `<dockflow_root>/UMDFP/bin/<tool>`
3. System `PATH` lookup

## Smoke Test (SMILES → PDBQT)

1. Ensure the config paths are correct for your setup (`config_case1.yaml`
   references `input_data/example.smi` relative to the working directory).
2. Run from the `test_files/ligprep_module_test/` directory:

```bash
python ligprep_pipeline.py config_case1.yaml
```

This converts the example SMILES file into PDBQT files (one per molecule) in
the `case1_pdbqt_single/` output directory.

Unlike some other pipelines in antioch-vs, this pipeline does not accept an
`--overwrite` flag. Instead, it automatically cleans the output directory at
the start of each run by removing and recreating it. This means re-running with
the same config will always produce a fresh set of output files, overwriting
any previous results in that directory.

## Smoke Test (UMF → PDBQT, skipping SMILES conversion)

To test only the UMF dump step (skipping smiles2umd/umfgen), create a config
with `source.type: "umf"` pointing to `input_data/example.umf`:

```yaml
source:
  type: "umf"
  path: "input_data/example.umf"

processing:
  target_format: "pdbqt"
  output_dir: "case_umf_pdbqt"
  num_blocks: 1
  split_files: true
  scratch_dir: "./scratch/"
```

Then run:

```bash
python ligprep_pipeline.py <that_config>.yaml
```

## Configuration Options

Copy `config_template.yaml` and edit the following sections as needed:

- **`umdfp:`** — Paths to UMDFP binaries (leave defaults if on PATH)
- **`source:`** — Input type (`smiles` or `umf`) and file path
- **`smiles_conversion:`** — CPU count, batch size, charge method (SMILES only)
- **`processing:`** — Target format (`pdbqt`, `mol2`, `sdf`, `smi`, `umd`),
  output directory, block count, split vs. combined files, scratch directory
- **`logging:`** — Log file path and verbosity

## Output Formats for Downstream Docking

Set `processing.target_format` based on your docking workflow:

- `pdbqt` → AutoDock Vina pipelines
- `mol2` → DOCK6 or workflows expecting Mol2 ligand libraries

## Syntax Checks

Validate the pipeline and config files without running the full conversion:

```bash
python3 -c 'import ast; ast.parse(open("ligprep_pipeline.py").read()); print("Python syntax OK")'
ruby -e 'require "yaml"; YAML.load_file("config_case1.yaml"); puts "YAML syntax OK"'
```
