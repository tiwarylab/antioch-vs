# Vina Pipeline Test Files

This folder contains public receptor and ligand files for smoke-testing the
Vina CPU and Vina-GPU pipelines.

The files come from the official AutoDock Vina basic docking tutorial for
imatinib docked to c-Abl from PDB entry 1IEP:

- `receptor_dir/1iep_receptor.pdbqt`
- `receptor_dir/1iep_receptorH.pdb`
- `ligand_dir/1iep_ligand.pdbqt`
- `source_files/1iep_ligand.sdf`

The tutorial uses this search box:

```text
center = 15.190 53.903 16.917
size   = 20 20 20
```

Because the prepared receptor PDBQT is provided, these smoke-test configs skip
receptor preparation and do not require an interface-residue CSV.
The ligand directory intentionally contains only `.pdbqt` files so batch
docking tools do not try to scan source-format files.

## Setup

Create and activate the docking environment from the repository root:

```bash
micromamba env create -f docking_module/vina_pipeline/vina_env.yml
micromamba activate docking_pipeline
```

If you use `mamba` or `conda` instead of `micromamba`, adjust the commands and
the `environment.mamba_exe` value in the example config files.

## Test the Vina CPU Pipeline

The CPU example is configured for native `vina` execution through the selected
Mamba environment:

```yaml
environment:
  mamba_env: "docking_pipeline"
  mamba_exe: "micromamba"
  vina_exe: "vina"
```

Run from the repository root:

```bash
python docking_module/vina_pipeline/vina_cpu_pipeline.py \
  --config test_files/docking_module_test/vina_test_files/example_vina_cpu_config.yaml
```

The pipeline writes the generated Vina config and docked output under
`test_files/docking_module_test/vina_test_files/docking_out_cpu`.

## Test the Vina-GPU Pipeline

Before running the GPU example, update `environment.vina_gpu_bin` in
`example_vina_gpu_config.yaml` to the absolute path of your Vina-GPU executable
or Singularity image.

### Native Vina-GPU Executable

Use this mode when Vina-GPU is installed as a direct executable on the host:

```yaml
environment:
  mamba_env: "docking_pipeline"
  mamba_exe: "micromamba"
  vina_gpu_bin: "/absolute/path/to/AutoDock-Vina-GPU-2-1"
  vina_gpu_cmd: "AutoDock-Vina-GPU-2-1"
  singularity_binds: []
```

When `vina_gpu_bin` does not end in `.sif`, the pipeline treats it as a native
binary and runs it directly from its containing directory.

### Singularity Container

Use this mode when Vina-GPU is packaged as a Singularity image:

```yaml
environment:
  mamba_env: "docking_pipeline"
  mamba_exe: "micromamba"
  vina_gpu_bin: "/absolute/path/to/Vina.sif"
  vina_gpu_cmd: "run_gpu"
  singularity_binds:
    - "/absolute/path/to/input-data"
    - "/absolute/path/to/OpenCL/vendors"
```

When `vina_gpu_bin` ends in `.sif`, the pipeline runs:

```bash
singularity exec --nv -B <binds> /absolute/path/to/Vina.sif run_gpu --config <generated-vina-config>
```

Include every host directory that the container needs to read or write in
`singularity_binds`. At minimum, bind the directories that contain the receptor,
ligands, output folder, and any OpenCL vendor files needed by your system.

Then run from the repository root:

```bash
python docking_module/vina_pipeline/vina_gpu_pipeline.py \
  --config test_files/docking_module_test/vina_test_files/example_vina_gpu_config.yaml
```

The pipeline writes the generated Vina-GPU config and docked output under
`test_files/docking_module_test/vina_test_files/docking_out_gpu`.

## Optional Vina-GPU Ligand Compatibility Split

Use the helper below before a Vina-GPU batch run when ligands were made with a
recent `mk_prepare_ligand.py`. It creates GPU-safe and CPU-safe input
directories, a `vina_unsupported/` directory for types neither built-in Vina
profile accepts, and per-file compatibility manifests without rewriting PDBQTs:

```bash
python docking_module/vina_pipeline/classify_vina_ligands.py \
  --input test_files/docking_module_test/vina_test_files/ligand_dir \
  --output-dir /tmp/vina_ligands_by_engine \
  --cpu-ligand-set cpu_only
```

Use `/tmp/vina_ligands_by_engine/vina_gpu_compatible` in the Vina-GPU YAML or
`/tmp/vina_ligands_by_engine/vina_cpu_compatible` in the Vina CPU YAML. The
`cpu_only` setting keeps those engine input directories disjoint; omit it when
you instead want an all-compatible CPU run.

## Optional Interface-Residue Extraction

The public 1IEP smoke test uses an explicit box, so interface-residue extraction
is not needed. For a different complex where you want automatic box placement,
first generate a CSV with:

```bash
python docking_module/vina_pipeline/get_interface_residues.py \
  --input /path/to/complex.pdb \
  --output /path/to/interface_residues.csv \
  --receptor <receptor-chain-id> \
  --ligand <ligand-chain-id>
```

Then set `box.center: "auto"`, `box.size: "auto"`, `paths.receptor_pdb`, and
`paths.target_residues_file` in the Vina config.

## Quick Syntax Checks

These checks validate the Python entry points and YAML config syntax without
running docking:

```bash
python3 -c 'import ast, pathlib; files=["docking_module/vina_pipeline/classify_vina_ligands.py","docking_module/vina_pipeline/get_interface_residues.py","docking_module/vina_pipeline/vina_cpu_pipeline.py","docking_module/vina_pipeline/vina_gpu_pipeline.py"]; [ast.parse(pathlib.Path(p).read_text(), filename=p) for p in files]; print("Python syntax OK")'
ruby -e 'require "yaml"; ARGV.each { |p| YAML.load_file(p); puts "YAML syntax OK: #{p}" }' test_files/docking_module_test/vina_test_files/example_vina_cpu_config.yaml test_files/docking_module_test/vina_test_files/example_vina_gpu_config.yaml
python3 -m unittest docking_module.vina_pipeline.test_classify_vina_ligands
```
