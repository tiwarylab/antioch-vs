# Docking Module Guide

This guide covers the working parts of the repository under `docking_module/`.
It does not cover `ctmd_module/`, which is currently a placeholder.

There are four runnable paths:

1. `receptor_prep/` for fetching PDB complexes and preparing receptor/ligand starting files
2. `dock6_pipeline/` for receptor prep with UCSF Chimera plus DOCK6 docking
3. `vina_pipeline/vina_cpu_pipeline.py` for AutoDock Vina CPU batch docking
4. `vina_pipeline/vina_gpu_pipeline.py` for AutoDock Vina GPU batch docking

## What You Need

The scripts in this repo orchestrate external tools. They do not install those tools for you.

Required for the DOCK6 path:

- UCSF Chimera available locally, through micromamba, or inside a Singularity `.sif`
- A DOCK6 installation with `bin/` and `parameters/`
- A receptor PDB file
- A reference ligand in `.mol2` format for defining the docking site, or a
  ANTIOCH-VS rank-1 pocket-definition YAML from `pocket_detection_module/`
- A ligand library in `.mol2` format

Required for the Vina paths:

- `micromamba` or another executable name you will place in the YAML config
- A Conda or Mamba environment with the packages from `vina_pipeline/vina_env.yml`
- AutoDock Vina CPU or Vina-GPU as a native executable, or as a Singularity `.sif`
- A receptor `.pdb` or prebuilt receptor `.pdbqt`
- A directory of ligand `.pdbqt` files

If you run containerized workflows on an HPC cluster, you may need to load
Singularity or Apptainer first, for example:

```bash
module load singularity
```

Python packages used by the Vina scripts:

- `numpy`
- `biopython`
- `pyyaml`
- `mdanalysis`
- `meeko`

## Repo Layout

- `receptor_prep/fetch_prepare_receptor.py` fetches a PDB complex, prepares the receptor, and extracts a ligand
- `dock6_pipeline/chimera_pipeline.py` prepares a receptor for DOCK6
- `dock6_pipeline/dock6_pipeline.py` runs DOCK6 sphere generation, box generation, grid generation, and flexible docking
- `vina_pipeline/get_interface_residues.py` extracts interface residues from a complex PDB into CSV
- `vina_pipeline/vina_cpu_pipeline.py` runs CPU batch docking
- `vina_pipeline/vina_gpu_pipeline.py` runs GPU batch docking

## Optional: Fetch and Prepare a Receptor/Ligand Pair

Use `receptor_prep/fetch_prepare_receptor.py` when you want to start from a PDB
ID or a local receptor-ligand complex PDB:

```bash
python docking_module/receptor_prep/fetch_prepare_receptor.py \
  --pdb-id 1xgj \
  --ligand-resname LIG \
  --output-dir receptor_prep_out
```

This writes a prepared receptor PDB, extracted ligand PDB, and ligand SDF. See
`docking_module/receptor_prep/README.md` for ligand chain/residue selectors and
dependency notes.

## Expected Inputs

### DOCK6 inputs

- `receptor.pdb`
  - Raw receptor structure used by the Chimera prep pipeline
- `crystal_ligand.mol2`
  - Reference ligand used to define the binding site for sphere selection
- `ligands_to_dock.mol2`
  - Single ligand or multi-molecule ligand library in Mol2 format

### Vina inputs

- `receptor.pdb` or `receptor.pdbqt`
  - Use `.pdb` if you want the pipeline to prepare receptor PDBQT
  - Use `.pdbqt` if you want to skip receptor preparation
- `ligands/`
  - Directory containing one or more ligand `.pdbqt` files
- `interface_residues.csv`
  - Needed only if `box.center` or `box.size` is set to `"auto"`
- `docking_boxes/<receptor>_best_pocket.yaml`
  - Optional pocket-detection output. Set `paths.pocket_definition` to this
    file to override `box.center` and `box.size` in either Vina runner.

### Interface residue CSV format

The Vina pipelines expect a CSV with this header:

```csv
Chain,ResidueID,ResidueName
A,45,TYR
A,46,GLY
C,102,LEU
```

You can generate that CSV with `get_interface_residues.py`.

## Environment Setup

### Vina environment

From the repo root:

```bash
micromamba create -f docking_module/vina_pipeline/vina_env.yml
```

If your executable is not named `micromamba`, update the `mamba_exe` field in the YAML config files.

### Chimera notes

The Chimera prep pipeline supports three execution modes:

- `local`: runs `chimera_cmd --nogui <generated_com_file>`
- `micromamba`: runs `micromamba run -n <mamba_env> <chimera_cmd> --nogui <generated_com_file>`
- `singularity`: runs `singularity exec ... <sif_image> chimera --nogui <generated_com_file>`

Use `dock6_pipeline/chimera_config.yaml` for the most complete configuration
template.

## Workflow A: DOCK6

This path has two stages:

1. Prepare the receptor with Chimera
2. Run DOCK6 with the prepared receptor files

### Step A1: Prepare the receptor with Chimera

From the repo root:

```bash
python docking_module/dock6_pipeline/chimera_pipeline.py \
  --config docking_module/dock6_pipeline/chimera_config.yaml
```

Minimal local-mode command:

```bash
python docking_module/dock6_pipeline/chimera_pipeline.py \
  --mode local \
  -i /absolute/path/to/receptor.pdb \
  -o /absolute/path/to/chimera_output \
  -d 5
```

Singularity-mode config example:

```yaml
input_pdb: "/absolute/path/to/receptor.pdb"
out_dir: "/absolute/path/to/chimera_output"
density: 5
execution_mode: "singularity"
sif_image: "/absolute/path/to/chimera_frozen.sif"
sif_binds:
  - "/absolute/path/to/project"
chimera_cmd: "chimera"
```

What this writes:

- `receptorname_wH.pdb` if hydrogens were added
- `receptorname_wH_wcH.mol2`
- `receptorname_surf_d5.dms`
- `receptorname_chimera_pipeline.log`

The log includes a config snapshot for the run, followed by command output.

Use the generated `.mol2` and `.dms` files in the DOCK6 config.

### Step A2: Fill in the DOCK6 YAML config

Start from `dock6_pipeline/dock6_config.yaml` and replace the placeholders:

```yaml
paths:
  receptor_dms: "/absolute/path/to/chimera_output/receptorname_surf_d5.dms"
  receptor_mol2: "/absolute/path/to/chimera_output/receptorname_wH_wcH.mol2"
  ref_ligand_mol2: "/absolute/path/to/crystal_ligand.mol2"
  pocket_definition: null
  ligand_file: "/absolute/path/to/ligands_to_dock.mol2"
  output_dir: "/absolute/path/to/dock6_output"

parameters:
  sphere_radius: 15.0
  box_margin: 8.0
  max_orientations: 1000

environment:
  dock6_dir: "/absolute/path/to/dock6"
```

Important:

- `dock6_dir` must be the top-level DOCK6 install directory
- The script expects `dock6_dir/bin/` and `dock6_dir/parameters/`
- `output_dir` must either be empty or you must run with `--overwrite`

When screening from a predicted pocket rather than a crystallographic ligand,
set `paths.pocket_definition` to one of the
`docking_boxes/<receptor>_best_pocket.yaml` files written by the pocket
detection module and leave `ref_ligand_mol2` blank. ANTIOCH-VS writes a one-atom
MOL2 pseudo-reference ligand at the predicted center, then runs the normal
DOCK6 `sphgen`, `sphere_selector`, and `showbox` stages. The default
`sphere_radius: 15.0` is used unless overridden in the docking YAML.

If your DOCK6 installation has an MPI-enabled executable, keep the setup steps
serial and launch only the final flexible docking step through MPI:

```yaml
environment:
  dock6_dir: "/absolute/path/to/dock6"
  dock6_executable: "/absolute/path/to/dock6/bin/dock6.mpi"
  dock6_launcher: "srun --mpi=pmix -n ${SLURM_NTASKS}"
```

For Slurm, request MPI tasks rather than only CPUs per task:

```bash
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
```

`dock6_launcher` is prepended only to the `dock6 -i flex.in -o flex.out`
command. The earlier `sphgen`, `sphere_selector`, `showbox`, and `grid` stages
still run once.

To run DOCK6 from a Singularity image, set `dock6_dir` to the installation
directory inside the image and provide the image path:

```yaml
environment:
  dock6_dir: "/opt/dock6"
  sif_image: "/absolute/path/to/dock6.sif"
  singularity_binds:
    - "/absolute/path/to/project"
```

When `sif_image` is set, every pipeline stage runs as `singularity exec ...`.
The internal commands default to their standard names: `sphgen`,
`sphere_selector`, `showbox`, `grid`, and `dock6`. For example, flexible
docking runs as:

```bash
singularity exec -B <binds> /absolute/path/to/dock6.sif dock6 -i flex.in -o flex.out
```

If an image exposes different names or paths, override only those entries:

```yaml
environment:
  dock6_dir: "/opt/dock6"
  sif_image: "/absolute/path/to/dock6.sif"
  container_commands:
    sphgen: "/opt/dock6/bin/sphgen"
    dock6: "/opt/dock6/bin/dock6"
```

`dock6_dir` must point to the installation directory inside the image so the
pipeline can reference its `parameters/` files.
The pipeline automatically binds the input and output directories; add any
other required host paths to `singularity_binds`. On an HPC system, you may
first need to run `module load singularity` or its site-specific equivalent.

### Step A3: Run DOCK6

From the repo root:

```bash
python docking_module/dock6_pipeline/dock6_pipeline.py \
  --config docking_module/dock6_pipeline/dock6_config.yaml
```

If reusing an existing output directory:

```bash
python docking_module/dock6_pipeline/dock6_pipeline.py \
  --config docking_module/dock6_pipeline/dock6_config.yaml \
  --overwrite
```

What this creates under your output directory:

- `01_spheres/`
- `02_gridbox/`
- `03_flex_dock/`
- `*_dock6_pipeline.log`

The log includes a config snapshot for the run, progress messages, and external
DOCK6 command output.

Main output to look for:

- `03_flex_dock/<prefix>_flex_scored.mol2`

### Restricting Orientations to Part of the Site (Critical Spheres)

Sphere selection defines a *region*; it cannot say that a pose must engage a particular spot. DOCK6's
critical-point filter can: spheres are tagged with a cluster number, and any orientation that does not
use a sphere from **every** tagged cluster is discarded. Typical use is to require a kinase ligand to
reach the hinge, or a GPCR ligand to reach the conserved aspartate, instead of drifting into a
neighbouring sub-pocket.

Add `parameters.critical_spheres`. One mapping is one cluster; a list gives several (at most 10):

```yaml
parameters:
  critical_radius: 3.0            # Å from the named atoms to a sphere CENTRE (default 3.0)
  critical_spheres:
    selection: "resid 930 932 and name N O"
```

Spheres are pseudo-atoms filling the pocket, so their centres sit roughly a hydrogen bond away from
the atoms you name: a sphere 2.9 Å from a backbone carbonyl is where a ligand's donor would sit.
Raising `critical_radius` much above 3 Å starts tagging spheres on the far side of those atoms, which
no ligand can reach.

A cluster is defined by exactly one of:

| Key | Meaning |
| --- | --- |
| `selection` | An [MDAnalysis selection](https://userguide.mdanalysis.org/stable/selections.html) against the receptor MOL2, e.g. `resid 930 932 and name N O`, `resid 113 and name OD1 OD2`, `resid 930 and backbone`. The selection language is the atom filter — `backbone`, `not backbone`, `name ...` all work. Hydrogens are always excluded. |
| `selection` + `reference` | The same, but the residue numbers are written in `reference`'s numbering (a PDB or MOL2) and mapped onto this receptor by sequence alignment. |
| `points` | A MOL2/PDB whose heavy atoms mark the site — a ligand fragment, a structural water, an ion. |
| `center` | Explicit `[x, y, z]`. |

Use `reference` when one residue list has to serve a whole ensemble. Receptors in one benchmark
rarely share a numbering scheme (a crystal may carry UniProt numbers while a predicted model starts
at 1), so a literal residue list silently points at the wrong residues:

```yaml
  critical_spheres:
    selection: "resid 930 932 and name N O"     # JAK2 hinge, in 3LPB numbering
    reference: /path/to/3LPB_prepared.mol2
```

Only residues whose type matches across the alignment are mapped, and an unmappable residue is an
error rather than a silent omission.

Notes:

- `critical_points yes` is set automatically when `critical_spheres` is configured. An explicit
  `extra_flex_args.critical_points` still wins, so spheres can be tagged without filtering on them.
- With several clusters a pose must satisfy **all** of them, which is a much harsher filter than one
  cluster. A sphere belongs to at most one cluster; where clusters overlap, the lower number wins.
- A cluster that tags no spheres is a fatal error: DOCK6 would otherwise reject every orientation and
  return no poses with nothing naming the cause. Tagging more than half the spheres only logs a
  warning, since the filter is then barely restricting anything.
- The filter constrains the **orientation match**, which in flexible docking places the anchor
  fragment. It requires the anchor near the tagged spheres, not a particular hydrogen bond. Expect
  some ligands to stop docking, and check coverage alongside enrichment.
- Counts, the radius and the resolved selections are recorded under `configuration.critical_spheres`
  in `dockflow_dock6_summary.json`.
- The selection is resolved when the pipeline is constructed, before sphgen runs, so a typo, a
  missing `reference` or a missing dependency fails immediately rather than after a surface has
  been built.

Needs MDAnalysis and Biopython, both already used elsewhere in this module.

### Reproducibility: What a Rerun Has to Match

DOCK6 is reproducible for a given receptor, parameter set and execution environment: two separate
runs of one config here produced byte-identical scored poses, because `simplex_random_seed` is a
fixed seed and orientation generation, pruning and clustering are deterministic.

Whether it is reproducible *across* execution environments has not been tested here. The rank count
set by `dock6_launcher` and the CPU a job lands on are both plausible ways for results to shift —
a different rank count changes how work is distributed, and a different CPU changes floating-point
rounding on near-tied poses — but neither has been measured, and a difference between two of this
project's runs that at first looked environmental turned out to be a parameter difference
(`simplex_final_min`) instead. Treat the environment as something to hold fixed and record, not as
a known source of variation.

None of this is visible in `parameters`, so every run records it under `configuration.execution` in
`dockflow_dock6_summary.json`:

```json
"execution": {
  "dock6_launcher": "srun --mpi=pmix -n ${SLURM_NTASKS}",
  "dock6_launcher_expanded": "srun --mpi=pmix -n 32",
  "mpi_ranks": 32,
  "slurm_ntasks": "32",
  "slurm_job_id": "236558",
  "hostname": "...",
  "cpu": "..."
}
```

`mpi_ranks` is read from the launcher after environment expansion; it is 1 when no launcher is set,
and null when the launcher cannot be parsed (rather than guessing). When comparing receptors or
conditions, keep `dock6_launcher` identical across every run in the comparison.

### DOCK6 Name Length Limits

DOCK6 reads file names and MOL2 title lines into fixed-width buffers. Past
roughly 80 characters, `sphgen` silently truncates the name it was handed and
`grid` segfaults on an over-long receptor title. Neither failure reports its own
cause, so the pipeline keeps what it generates short:

- Intermediates use fixed names — `01_spheres/receptor.sph`,
  `02_gridbox/receptor.box.pdb`, `02_gridbox/receptor_grid*`. Each run has its
  own working directory, so these do not need to be unique, and no receptor file
  name can lengthen them. Earlier runs named these after the receptor prefix.
- The log and the final `03_flex_dock/<prefix>_flex_scored.mol2` keep the
  descriptive prefix, so the documented output contract is unchanged.
- If the receptor MOL2 title exceeds the limit, a copy with the title
  `receptor` is staged for grid generation and the run records a warning in
  `dockflow_dock6_summary.json`. The receptor file itself is never modified, and
  a title within the limit is symlinked as before. Ligand MOL2 titles are never
  rewritten — those are the ligand IDs the extractor reads back as `Ligand`.

The threshold is `parameters.dock6_name_limit` (default 80).

### DOCK6 Sampling Presets

DOCK6 sampling has a direct speed-versus-thoroughness tradeoff. Recent DOCK6
benchmarking with hierarchical database search found that increasing
orientation sampling and keeping more candidate poses can improve enrichment
or pose recovery, but increases runtime. Our pipeline currently uses the
classic DOCK6 flexible anchor-and-grow workflow, so these are practical
presets rather than exact HDB parameter translations.

Fast HTVS:

```yaml
parameters:
  max_orientations: 1200
  extra_flex_args:
    num_preclustered_conformers: 1
    num_final_scored_poses: 1
    write_orientations: no
    simplex_final_min: no
```

Balanced HTVS:

```yaml
parameters:
  max_orientations: 3600
  extra_flex_args:
    num_preclustered_conformers: 10
    num_final_scored_poses: 1
    write_orientations: no
    simplex_final_min: yes
    simplex_max_cycles: 100
```

Thorough redocking/debugging:

```yaml
parameters:
  max_orientations: 5000
  extra_flex_args:
    num_preclustered_conformers: 50
    num_final_scored_poses: 10
    write_orientations: yes
    simplex_final_min: yes
    simplex_max_cycles: 200
    simplex_anchor_max_iterations: 1000
    simplex_grow_max_iterations: 1000
```

For large screens, start with the fast or balanced preset. Use the thorough
preset for crystal-ligand redocking, difficult systems, or when you want to
inspect whether failures are sampling failures or scoring failures. Ligand
preparation remains critical: protonation, charges, and input conformer quality
can matter as much as the DOCK6 sampling parameters.

## Workflow B: AutoDock Vina CPU

### Step B1: Optional interface residue extraction

If you want automatic box calculation, generate interface residues from a receptor-ligand complex PDB:

```bash
micromamba run -n docking_pipeline python docking_module/vina_pipeline/get_interface_residues.py \
  -i /absolute/path/to/complex.pdb \
  -o /absolute/path/to/interface_residues.csv \
  -c 5.0 \
  -r A \
  -l C
```

Adjust chain IDs to match your structure.

### Step B2: Fill in the Vina CPU config

Start from `vina_pipeline/vina_cpu_config.yaml`:

```yaml
paths:
  receptor_pdbqt: ""
  receptor_pdb: "/absolute/path/to/raw_receptor.pdb"
  ligand_dir: "/absolute/path/to/ligands"
  output_dir: "/absolute/path/to/vina_cpu_output"
  target_residues_file: "/absolute/path/to/interface_residues.csv"

box:
  center: "auto"
  size: "auto"
  padding: [0, 0, 15.0]
  target_chain: "A"

docking:
  cpu: 32
  exhaustiveness: 40
  num_modes: 10
  min_rmsd: 3
  energy_range: 3

environment:
  mamba_env: "docking_pipeline"
  mamba_exe: "micromamba"
  vina_exe: "/absolute/path/to/vina"
  vina_cmd: "vina"
  singularity_binds: []
```

Notes:

- If `receptor_pdbqt` is set, the script skips receptor preparation
- If `box.center` and `box.size` are explicit values, `target_residues_file` is not needed
- `ligand_dir` should contain `.pdbqt` ligands
- If `vina_exe` ends in `.sif`, the CPU pipeline uses Singularity and runs `vina_cmd` inside the container
- If `vina_exe` does not end in `.sif`, the CPU pipeline runs it through the configured Mamba environment

### Step B3: Run Vina CPU

From the repo root:

```bash
python docking_module/vina_pipeline/vina_cpu_pipeline.py \
  --config docking_module/vina_pipeline/vina_cpu_config.yaml
```

What this writes:

- `<prefix>_vina_config.txt`
- `batch_docking.log`
- Docked ligand outputs in the configured output directory

`batch_docking.log` includes a config snapshot for the run, progress messages,
and Vina command output.

## Workflow C: AutoDock Vina GPU

### Step C1: Optional interface residue extraction

Use the same `get_interface_residues.py` step shown above if you want auto box calculation.

### Step C1a: Optional ligand compatibility split

Recent `mk_prepare_ligand.py` versions can emit macrocycle closure pseudo atom
types such as `CG0`, `CG1`, `G0`, and `G1`. Current AutoDock Vina CPU accepts
these PDBQTs, but Vina-GPU 2.1 rejects them. Before a GPU screen, separate an
existing prepared-ligand file or directory without changing its chemistry:

```bash
python docking_module/vina_pipeline/classify_vina_ligands.py \
  --input /absolute/path/to/prepared_ligands \
  --output-dir /absolute/path/to/prepared_ligands_by_engine \
  --cpu-ligand-set cpu_only
```

This writes four directories and two provenance files:

- `vina_gpu_compatible/`: PDBQTs whose atom types match the built-in legacy
  Vina-GPU 2.1 profile. Set this as the GPU YAML's `paths.ligand_dir`.
- `vina_cpu_compatible/`: PDBQTs selected for the built-in current AutoDock
  Vina CPU profile. With the default `--cpu-ligand-set all_compatible`, this
  includes GPU-compatible and CPU-only inputs. With `cpu_only`, it contains
  only inputs rejected by GPU but accepted by CPU, so the engine input sets do
  not overlap. Set this as the CPU YAML's `paths.ligand_dir`.
- `vina_unsupported/`: syntactically readable PDBQTs whose atom types are not
  accepted by either built-in profile, such as `B`. Re-prepare these ligands or
  use an engine with explicit support for their atom types.
- `invalid_pdbqt/`: files without readable PDBQT atom records; do not dock
  these until they are prepared again.
- `vina_ligand_compatibility.csv` and
  `vina_ligand_compatibility_summary.json`: per-file atom-type decisions and
  run provenance.

The default is `--materialize copy`, making both output directories usable when
only the classified directory is bind-mounted into a Vina-GPU container. Use
`--materialize hardlink` to save space on the same filesystem, or `symlink`
only when the source directory is also visible inside the container.

Use `--cpu-ligand-set cpu_only` when GPU is the primary screen and CPU should
run only the ligands GPU cannot parse. The default, `all_compatible`, retains
the complete CPU-compatible set for an independent CPU screen.

The built-in profiles target Vina-GPU 2.1 and current AutoDock Vina CPU. For a
different build, pass its exact comma-separated parser allowlist through
`--gpu-allowed-types` or `--cpu-allowed-types`; the helper records both in the
summary JSON.

### Step C2: Fill in the Vina GPU config

Start from `vina_pipeline/vina_gpu_config.yaml`:

```yaml
paths:
  receptor_pdbqt: ""
  receptor_pdb: "/absolute/path/to/raw_receptor.pdb"
  ligand_dir: "/absolute/path/to/ligands"
  output_dir: "/absolute/path/to/vina_gpu_output"
  target_residues_file: "/absolute/path/to/interface_residues.csv"

box:
  center: "auto"
  size: "auto"
  padding: [0, 0, 15.0]
  target_chain: "C"

docking:
  thread: 8000
  search_depth: 40
  num_modes: 10

environment:
  mamba_env: "docking_pipeline"
  mamba_exe: "micromamba"
  vina_gpu_bin: "/absolute/path/to/AutoDock-Vina-GPU-2-1"
  vina_gpu_cmd: "AutoDock-Vina-GPU-2-1"
  singularity_binds: []
```

Notes:

- `vina_gpu_bin` should point to the executable itself, not just its directory
- If `vina_gpu_bin` ends in `.sif`, the GPU pipeline uses `singularity exec --nv` and runs `vina_gpu_cmd` inside the container
- If `vina_gpu_bin` does not end in `.sif`, the script treats it as a native executable and changes into its directory before launch
- `ligand_dir` should contain `.pdbqt` ligands

Singularity GPU example:

```yaml
environment:
  mamba_env: "docking_pipeline"
  mamba_exe: "micromamba"
  vina_gpu_bin: "/absolute/path/to/Vina.sif"
  vina_gpu_cmd: "run_gpu"
  singularity_binds:
    - "/absolute/path/to/project"
    - "/absolute/path/to/OpenCL/vendors"
```

### Step C3: Run Vina GPU

From the repo root:

```bash
python docking_module/vina_pipeline/vina_gpu_pipeline.py \
  --config docking_module/vina_pipeline/vina_gpu_config.yaml
```

What this writes:

- `<prefix>_vina_config.txt`
- `docking.log`
- Docked ligand outputs in the configured output directory

`docking.log` includes a config snapshot for the run and pipeline progress
messages.

## Post-Docking Result Extraction

After a Vina or DOCK6 run, use the shared post-docking analysis entry point:

```bash
python docking_module/post_docking_analysis/extract_docking_results.py \
  --software dock6 \
  -i /absolute/path/to/dock6_output \
  --out_dir /absolute/path/to/dock6_summary \
  --out_parquet dock6_summary.parquet \
  --out_csv top_hits.csv \
  --top_n 20 \
  -c 8
```

Set `--software vina` for Vina CPU/Vina-GPU outputs or `--software dock6` for
DOCK6 outputs. The master script forwards the remaining flags to the matching
extractor.

Create the optional extraction/export environment without modifying your main
docking environment:

```bash
micromamba create -f docking_module/post_docking_analysis/environment.yaml
micromamba activate dockflow-post-docking
```

This environment includes `pandas`, `pyarrow`, `tqdm`, Open Babel, and Meeko.
Parquet writing requires `pyarrow` or `fastparquet`. DOCK6 SDF export requires
Open Babel's `obabel`; Vina SDF export requires Meeko's `mk_export.py`.

### Vina Extraction

The Vina extractor scans for `*out.pdbqt` files by default. For each file, it
reads the first `REMARK VINA RESULT:` line, records the top Vina score, stores
the absolute file path, and sorts ligands from best to worst score. More
negative scores are ranked first.

```bash
python docking_module/post_docking_analysis/extract_docking_results.py \
  --software vina \
  -i /absolute/path/to/vina_output \
  --out_dir /absolute/path/to/vina_summary \
  --out_parquet vina_summary.parquet \
  --out_csv top_hits.csv \
  --export_sdf \
  --sdf_n 20 \
  -c 8
```

### DOCK6 Extraction

The DOCK6 extractor scans for `*_flex_scored*.mol2` files by default. Each
scored Mol2 block is parsed as one pose. The metadata written by DOCK6 before
each pose, including `Grid_Score`, `Grid_vdw_energy`, `Grid_es_energy`,
`HA_RMSDs`, ligand properties, source file, and Mol2 line range, is stored in
the output table. More negative `Grid_Score` values are ranked first.

```bash
python docking_module/post_docking_analysis/extract_docking_results.py \
  --software dock6 \
  -i /absolute/path/to/dock6_output \
  --out_dir /absolute/path/to/dock6_summary \
  --out_parquet dock6_summary.parquet \
  --out_csv top_ligands.csv \
  --best_per_ligand \
  --export_mol2 \
  --export_sdf \
  --sdf_n 20 \
  -c 8
```

For DOCK6 multi-Mol2 inputs, multiple scored poses or conformers can share the
same ligand name. The master parquet keeps every scored pose. Add
`--best_per_ligand` when you want the terminal table, CSV, and exported
structures to include only the best-scoring pose per ligand.

For large runs organized into subfolders, process each subfolder as a batch and
then merge the results:

```bash
python docking_module/post_docking_analysis/extract_docking_results.py \
  --software dock6 \
  -i /absolute/path/to/dock6_output \
  --out_dir /absolute/path/to/dock6_summary \
  --batch_by_subdirs \
  --batch_discovery recursive_leaf \
  --out_csv top_hits.csv \
  --csv_n 1000 \
  -c 16
```

## Common Failure Points

- Missing external binaries
  - `chimera`, `vina`, `mk_prepare_receptor.py`, `obabel`, or the Vina-GPU executable are not on the system or not referenced correctly
- Wrong DOCK6 root path
  - `dock6_dir` must contain both `bin/` and `parameters/`
- Wrong ligand formats
  - DOCK6 expects Mol2 input
  - Vina expects ligand `.pdbqt` files in a directory
- Singularity not available in the shell
  - On HPC systems, run `module load singularity` or the site-specific equivalent before using `.sif` images
- Missing container bind paths
  - Add receptor, ligand, output, and support directories to `sif_binds` or `singularity_binds`
- Auto box with no residue CSV
  - If Vina uses `"auto"` box settings, the residue CSV must exist and match chain IDs in the receptor PDB
- Non-empty DOCK6 output folder
  - Use `--overwrite` if you want the script to clear previous stage folders

## Quick Command Summary

From the repo root:

```bash
# 1. Prepare a DOCK6 receptor
python docking_module/dock6_pipeline/chimera_pipeline.py \
  --config docking_module/dock6_pipeline/chimera_config.yaml

# 2. Run DOCK6
python docking_module/dock6_pipeline/dock6_pipeline.py \
  --config docking_module/dock6_pipeline/dock6_config.yaml

# 3. Generate interface residues for Vina auto-boxing
micromamba run -n docking_pipeline python docking_module/vina_pipeline/get_interface_residues.py \
  -i /absolute/path/to/complex.pdb \
  -o /absolute/path/to/interface_residues.csv \
  -c 5.0 \
  -r A \
  -l C

# 4. Run Vina CPU
python docking_module/vina_pipeline/vina_cpu_pipeline.py \
  --config docking_module/vina_pipeline/vina_cpu_config.yaml

# 5. Run Vina GPU
python docking_module/vina_pipeline/vina_gpu_pipeline.py \
  --config docking_module/vina_pipeline/vina_gpu_config.yaml
```
