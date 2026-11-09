# Structure Prediction Module

This module separates structure prediction from AF2RAVE ensemble analysis.
Prediction runners generate structures in backend-specific output folders; the
AF2RAVE runner then accepts one or more PDB folders and analyzes the combined
ensemble.

## Prediction Runners

Each runner uses its own YAML and supports a `runs:` list. Every entry writes to
`paths.output_dir/runs/<name>/`, keeping different reduced-MSA settings, seeds,
or sampling settings separate until you intentionally combine them for AF2RAVE.

```bash
python structure_prediction_module/prediction/run_localcolabfold.py \
  --config structure_prediction_module/prediction/localcolabfold_config.yaml

python structure_prediction_module/prediction/run_alphafold2.py \
  --config structure_prediction_module/prediction/alphafold2_config.yaml

python structure_prediction_module/prediction/run_alphafold3.py \
  --config structure_prediction_module/prediction/alphafold3_config.yaml
```

- `run_localcolabfold.py` uses a LocalColabFold CLI executable and can generate
  one deep MSA before launching multiple reduced-MSA prediction runs.
- `run_alphafold2.py` launches native AlphaFold 2 with Apptainer, a user-supplied
  container, and user-supplied database paths.
- `run_alphafold3.py` launches AlphaFold 3 with Apptainer from its JSON input and
  converts output mmCIF structures to PDB with Biopython.

The runners do not install prediction software, models, or databases. Set those
paths in a run-specific YAML rather than editing a ANTIOCH-VS template.

All prediction runners write `<runner>.log` in `paths.output_dir` by default.
Set `execution.log_file` to choose another filename. Relative paths are placed
under `paths.output_dir`; absolute paths are written at the specified location.

```yaml
execution:
  overwrite: false
  log_file: "logs/localcolabfold_rmsa.log"
```

For LocalColabFold, `execution.overwrite: true` also passes
`--overwrite-existing-results` to ColabFold. Keep it `false` by default: it
prevents replacement of existing run outputs and deep A3M files.

## LocalColabFold Inputs and Model Parameters

The default `localcolabfold.mode: predict` accepts a FASTA file or inline
sequence, generates one deep A3M, then reuses that alignment for every
configured reduced-MSA prediction run. It never performs one MSA query per
`max_msa` setting. To predict from a precomputed `.a3m`, use
`mode: prediction_only` instead. ANTIOCH-VS copies supplied A3Ms, retaining their
filenames, into `<output_dir>/localcolabfold_a3m_input/` and launches
LocalColabFold from that A3M-only directory. It refuses to run if that
controlled workspace contains a FASTA file. The first sequence in the `.a3m`
must be the query sequence.

```yaml
paths:
  sequence_file: "/project/msas/jak2.a3m"
  output_dir: "/project/predictions/localcolabfold"

localcolabfold:
  mode: "prediction_only"
```

`msa_mode` is omitted for `.a3m` prediction input. The per-run `max_msa` value
can still be used to produce several reduced-MSA predictions from the same deep
alignment.

### MSA-only runs

To generate an alignment without predicting structures, supply an inline
sequence or FASTA input and set `localcolabfold.mode: msa_only`. ANTIOCH-VS passes
`--msa-only`, writes one deep A3M directly in `output_dir`, and omits
prediction-only options including `max_msa`, seeds, recycles, model type,
dropout, and `--save-all`. MSA-only mode requires exactly one `runs` entry;
the entry is only a label for the run summary.

Run only one `localcolabfold.mode` per YAML. Treat MSA generation and
prediction as separate configurations so a mixed run list cannot accidentally
query an MSA server while another entry is intended to reuse an A3M.

```yaml
paths:
  sequence_file: "/project/sequences/jak2.fasta"
  output_dir: "/project/localcolabfold_msas"

localcolabfold:
  mode: "msa_only"
  msa_mode: "mmseqs2_uniref_env"

runs:
  - name: "deep_msa"
```

Do not use `.a3m` input with `msa_only`: that file is already an alignment and
is meant for the prediction workflow. `max_msa` does not make a reduced A3M;
it is applied later, while ColabFold constructs prediction features.

### Continue From MSA-only Output

To predict structures from a previous MSA-only run, set
`localcolabfold.mode: prediction_only`. In this mode, `paths.output_dir` is the
prediction destination, while `paths.sequence_file` must point to one deep
`.a3m` file, produced by MSA-only mode or supplied by the user. ANTIOCH-VS copies
that A3M, retaining its filename, into `output_dir/localcolabfold_a3m_input/`,
runs LocalColabFold from that A3M-only workspace, and reuses it for every YAML
`runs[].name`. Each reduced-MSA prediction is written to
`output_dir/runs/<name>/`. It does not pass `--msa-mode`, so no new MSA
retrieval is requested.

Each prediction directory includes `dockflow_a3m_provenance.json`, which
records `input_mode: precomputed_a3m`, `msa_server_queried: false`, the staged
A3M path, and its SHA-256 digest. This is the ANTIOCH-VS provenance record for
the MSA source. ColabFold's own `config.json` may retain its internal default
`msa_mode` string even when it consumes a supplied A3M.

### End-to-end Prediction

Use `localcolabfold.mode: predict` when starting from a sequence or FASTA and
want both stages in one invocation. ANTIOCH-VS first executes the same deep-MSA
phase as `msa_only`, then automatically stages the generated A3M and performs
the same prediction phase as `prediction_only` for every `runs` entry. This is
the usual mode for an AF2RAVE reduced-MSA ensemble.

```yaml
paths:
  sequence_file: "/project/sequences/jak2.fasta"
  output_dir: "/project/localcolabfold_predictions"

localcolabfold:
  mode: "predict"
  msa_mode: "mmseqs2_uniref_env"
  num_seeds: 32
  num_recycles: 1
  model_type: "alphafold2_ptm"

runs:
  - name: "rmsa_8_16"
    settings:
      max_msa: "8:16"
  - name: "rmsa_16_32"
    settings:
      max_msa: "16:32"
```

ANTIOCH-VS refuses to run if the controlled A3M workspace contains a FASTA. Any
FASTA found directly in a pre-existing prediction run directory is renamed to
`<name>.fasta.disabled` (with a numeric suffix if needed), and a warning is
written to the LocalColabFold log. This prevents a stale FASTA from being used
as a new sequence query.

```yaml
paths:
  sequence: null
  sequence_file: "/project/localcolabfold_msas/jak2.a3m"
  output_dir: "/project/localcolabfold_msas"

localcolabfold:
  mode: "prediction_only"
  num_seeds: 32
  num_recycles: 1
  use_dropout: true
  save_all: true
  model_type: "alphafold2_ptm"
  extra_args:
    - "--data"
    - "/project/colabfold_cache"

runs:
  - name: "msa_8_16"
    settings:
      max_msa: "8:16"
  - name: "msa_16_32"
    settings:
      max_msa: "16:32"
```

To reuse an existing LocalColabFold model-parameter cache, pass the cache root
with `--data` via `localcolabfold.extra_args`:

```yaml
localcolabfold:
  extra_args:
    - "--data"
    - "/project/colabfold_cache"
```

The path supplied to `--data` is the parent directory that contains `params/`,
not the `params/` directory itself. For example, with model files at
`/project/colabfold_cache/params/params_model_1_ptm.npz`, use
`/project/colabfold_cache`. Supplying the nested `params/` path makes
LocalColabFold create and populate an unintended `params/params/` directory.

## AF2RAVE Analysis

Run AF2RAVE after prediction:

```bash
python structure_prediction_module/af2rave_pipeline/af2rave_pipeline.py \
  --config structure_prediction_module/af2rave_pipeline/af2rave_config.yaml
```

The AF2RAVE runner requires AF2RAVE, NumPy, Matplotlib, and Biopython in the
active Python environment. A separate environment is recommended:

```bash
micromamba create -n af2rave_env python=3.10
micromamba activate af2rave_env
pip install git+https://github.com/tiwarylab/af2rave.git biopython
```

## AF2RAVE Inputs

Provide one of these inputs:

- `paths.structures_dir`: one PDB directory.
- `paths.structure_dirs`: a YAML list of PDB directories.
- `paths.structure_dirs_file`: a text file with one PDB directory path per line.

AF2RAVE copies all discovered PDBs to `paths.output_dir/structures/` before
analysis, assigning short deterministic collection IDs so same-named PDBs from
recursive subdirectories or multiple source folders cannot overwrite one
another. `structure_collection_manifest.csv` records each source and collected
filename.

### Input Order and Clustering Reproducibility

The order of directories in `paths.structure_dirs` or
`paths.structure_dirs_file` is preserved while structures are collected. Keep
that order fixed when reproducing an analysis. Although a fixed
`clustering.randomseed` makes a run reproducible for a particular input order,
it does not make clustering invariant to a reordered structure list. K-means++
can choose different initial structures after a row-order change and converge
to a different local solution. AF2RAVE regular-space clustering is likewise
order-sensitive because it selects centers greedily. PDBs within each input
directory are sorted consistently by filename.

For a published or benchmarked ensemble, retain the original directory-list
file with the analysis outputs and reuse it unchanged for later comparisons.

```yaml
paths:
  structure_dirs:
    - "/project/predictions/localcolabfold/runs/rmsa_8_16"
    - "/project/predictions/localcolabfold/runs/rmsa_16_32"
  output_dir: "/project/af2rave_output"
```

AF2RAVE needs PDB input. The AlphaFold 3 runner writes normalized PDBs under
each run's `structures/` directory for this handoff.

## AF2RAVE Options

RMSD and peptide-bond filters are opt-in: set a cutoff to enable one, or set
it to `null` to leave it off. Steric-clash filtering defaults to AF2RAVE's
`1.0` A minimum nonbonded heavy-atom distance; set
`steric_clash_cutoff: null` to disable it. The older `nonbonded_cutoff` key is
still accepted as an alias when `steric_clash_cutoff` is absent.
Expensive nonbonded and peptide-bond statistics are calculated only when their
filters are enabled. Set `analysis.write_full_feature_stats: true` to calculate
those diagnostics even when they are not used for filtering.

The documented `filtering` selection and cutoff keys are authoritative. Do not
repeat `selection`, `rmsd_cutoff`, `min_non_bonded_cutoff`, `mean_cutoff`, or
`std_cutoff` inside the corresponding `*_filter_kwargs` dictionary; ANTIOCH-VS
rejects such conflicting settings so the applied filters and reported
diagnostics remain consistent.

Set `paths.reference_structure` to a specific PDB when you know the AF2RAVE
reference for RMSD alignment, RMSD filtering, and regular-space clustering. If
it is `null`, `paths.reference_selection: "best_plddt"` chooses the structure
with the highest residue-averaged pLDDT from PDB B-factor values. Set
`paths.reference_selection: "first"` to use the first structure instead.

`best_plddt` is appropriate only when every input PDB stores pLDDT in its
B-factor column. ANTIOCH-VS requires an explicit common scale through
`paths.plddt_bfactor_scale`: use `"0_100"` for standard AlphaFold/ColabFold
PDBs or `"0_1"` for unit-interval values, which are normalized before ranking.
An invalid scale setting always stops the run. When `best_plddt` is used,
values outside the configured range, including mixed scales, also stop the run.
With an explicit reference, invalid B-factor values are reported per structure
and recorded as blank pLDDT scores while valid structures remain in the table.
Real experimental B-factors can overlap the pLDDT numeric range and cannot be
identified reliably from PDB columns alone; for experimental or mixed inputs,
provide `paths.reference_structure` explicitly. The complete sorted ranking is
written to `plddt_scores.csv`.

By default pLDDT is averaged across all protein residues. To choose a reference
based on a folded core instead, set `paths.plddt_residue_ranges`, for example:

```yaml
paths:
  plddt_residue_ranges: ["41-187"]
```

For chain-specific structures, use entries such as `"A:41-187"`.

Set `clustering.method: "regular_space"` for AF2RAVE's regular-space clustering,
or `clustering.method: "kmeans"` for sklearn k-means. K-means representatives
are real PDBs closest to their centroids, so downstream workflows still use
physical structures.

When `clustering.copy_cluster_centers: true`, ANTIOCH-VS exports every selected
center in the coordinate frame of the AF2RAVE pipeline reference by default:
`paths.reference_structure`, the highest-pLDDT fallback, or the first loaded
structure. This affects only the copied PDBs in `cluster_center_structures/`;
it does not alter AF2RAVE feature calculation or either clustering method.
Alignment uses common C-alpha atoms and applies the resulting rigid
transformation to every atom in each exported PDB. Set
`clustering.align_cluster_centers: false` to retain original coordinates, or
set `clustering.center_alignment_reference: first_center` to use the first
ranked cluster center instead. The `cluster_center_structures.csv` manifest
records the reference structure, aligned C-alpha atom count, and
`alignment_common_ca_rmsd_A` for each center.

Aligned center exports are intended for predicted protein structures. Biopython
`PDBIO` preserves atom coordinates and B-factors but does not preserve PDB
records such as `CRYST1`, `CONECT`, or `SSBOND`; use a dedicated complex
preparation workflow rather than these exports when retaining a holo ligand or
crystallographic metadata is required.

`pca_clusters.png` is calculated from the same top-ranked feature subset used
for clustering (`features.n_features`), so the displayed separation and center
positions represent the clustering space rather than all available features.

For regular-space clustering, `max_centers` is the AF2RAVE search capacity and
defaults to 100. It is not a requested center count: AF2RAVE raises an error
when the capacity is reached. To request a fixed number of centers, set
`clustering.target_centers` and keep `max_centers` larger than that target.
ANTIOCH-VS starts at `dist_cutoff` (or its `min_dist` alias), increases or
decreases the distance to bracket the requested count, then performs a
deterministic binary search. After finding the requested count, it scans the
upper bracket at an effective resolution of the larger of
`target_search_min_dist_tolerance` and
`target_search_relative_tolerance * min_dist`, then selects the largest sampled
`min_dist` that still yields the requested number of centers. This favors the
most separated representatives at the requested count. The default 0.05 A
absolute floor and 1% relative tolerance prevent excessive full clustering
passes when feature-space distances are large. The scan, rather than a second
binary search, is necessary because greedy regular-space center counts can make
small non-monotonic jumps. It writes every attempt to
`regular_space_target_search.csv`; `target_search_max_evaluations` defaults to
1,000 and caps the scan. The run stops with a clear error if no tested distance
can produce exactly the requested count.

```yaml
clustering:
  method: "regular_space"
  dist_cutoff: 3.5
  max_centers: 100
  target_centers: 8
```

## Main Outputs

- `structures/`: collected PDB structures used by AF2RAVE.
- `structure_collection_manifest.csv`: source directory and collected name for
  every input PDB.
- `plddt_scores.csv`: residue-averaged pLDDT for each collected structure.
- `feature_stats.csv` and `feature_stats.png`: filtering diagnostics.
- `filtered_structures.txt`, `selected_features.csv`, and
  `feature_cv_histogram.png`.
- `cluster_centers.txt`, `cluster_assignments.csv` for k-means, and
  `pca_clusters.png`.
- `regular_space_target_search.csv` when `clustering.target_centers` is set.
- `simulation_boxes/` when `simulation_box.enabled` is true.
- `af2rave_pipeline.log`, `af2rave_report.txt`, and a JSON run summary.

Successful JSON summaries list only artifacts enabled by the active
configuration and present on disk; failed summaries report only the output
directory and current log, avoiding stale outputs from an earlier run. A feature
cache loaded through `reuse_feature_cache` is recorded as an input; only a
cache written during the current run appears in `outputs`. The
`reproducibility` section records UTC time, Python/AF2RAVE/NumPy/scikit-learn/
MDTraj/Biopython/matplotlib versions, and a SHA-256 fingerprint of the
source-relative collected-structure order. Regular-space clustering is
deterministic but order-sensitive, so its seed is recorded as `not_applicable`;
k-means records its effective random seed instead.

When AF2RAVE was installed without package metadata, its version field falls
back to `af2rave.__version__` or a SHA-256 hash of its source module. This
source-content identity is also used for feature-cache validation.

Feature-cache reuse is guarded by the ordered filtered PDB list, feature
selection, `rank_feature_kwargs`, and AF2RAVE version. A legacy cache without
this metadata, or a cache generated with a different feature definition or
AF2RAVE version, is rejected; regenerate it with `reuse_feature_cache: false`.

Runnable AF2RAVE examples live in:

```text
test_files/structure_prediction_module_test/af2rave_test_files/
```

## Unbiased MD

After clustering, build an AMBER unbiased-MD system for each cluster center:

```bash
python structure_prediction_module/af2rave_pipeline/run_unbiased_md.py \
  --config structure_prediction_module/af2rave_pipeline/unbiased_md_config.yaml
```

`run_unbiased_md.py` wraps `ctmd_module/md_submodule/setup_MD_from_PDB.sh` (the shared
soluble/membrane system builder) in a YAML-driven runner. It is build-only: it builds
`system.prmtop`/`system.inpcrd` per center and writes `submit_md.sh` for the production runs;
it does not run MD itself.

A single `membrane.enabled` flag selects the box: `false` = octahedral TIP3P (soluble),
`true` = POPC bilayer (packmol-memgen + Lipid21). Target-specific prep is passed through as the
`DOCKFLOW_*` environment variables documented in `ctmd_module/md_submodule/README.md` (pH, bulk
salt, protected residues, required disulfides with the eligibility gate, and optional loop
truncation with ACE/NME caps). The soluble OpenMM alternative remains `af2rave_md_pipeline.py`.

```yaml
paths:
  cluster_centers_dir: "/project/af2rave_output/cluster_center_structures"
  output_dir: "/project/unbiased_md_output"
membrane:
  enabled: true
prep:
  required_disulfides: ["A:75-A:160", "A:153-A:159"]
  truncate_ranges: ["200-231"]
```

Each center is classified built / skipped_ineligible / failed in `unbiased_md_manifest.csv`,
alongside `unbiased_md_summary.json` and `unbiased_md_pipeline.log`. A model missing a required
disulfide is skipped cleanly, not failed. The emitted `submit_md.sh` runs production for every
built system; set `launcher.md_run_command` to your per-system MD runner (01_Min on CPU,
02_Min2..prod on GPU) and `launcher.gpus_per_job` to pair systems onto GPUs. It submits only when
`launcher.submit: true`.

## SPIB Stage

Turn the unbiased-MD trajectories into SPIB reaction coordinates and metastable states:

```bash
python structure_prediction_module/af2rave_pipeline/run_spib_stage.py \
  --config structure_prediction_module/af2rave_pipeline/spib_stage_config.yaml
```

`run_spib_stage.py` orchestrates the whole chain: stage each trajectory (autoimage + strip),
optionally relabel, select the top-N features, generate COLVAR files, run the SPIB time-lag
sweep, and run the analysis. For membrane systems set `membrane.enabled: true` to strip the
lipids (POPC) as well as water. If the MD used truncation, set `features.drop_ranges` to drop the
cut residues from the feature set (equivalent to re-ranking on the truncated system, since the
coefficient of variation is per-pair independent) and `relabel.enabled: true` to restore the
original numbering on each `only_protein.pdb` (tleap renumbers the topology, so the feature names
would not otherwise map).

```yaml
paths:
  md_output_dir: "/project/unbiased_md_output"
  features_csv: "/project/af2rave_output/selected_features.csv"
membrane:
  enabled: true
features:
  n_features: 200
  drop_ranges: ["200-231"]
relabel:
  enabled: true        # split/offset auto-derived from a single drop_range
```

The COLVAR step needs MDAnalysis + PLUMED; SPIB and the analysis need AF2RAVE (and, for the
analysis, mdtraj + scikit-learn + matplotlib). When those live in different environments, set the
`environment.colvar_python` / `environment.spib_python` / `environment.cpptraj` command prefixes
(for example `micromamba run -n metad_env python`). Outputs land under `paths.output_dir`
(default: `md_output_dir`): per-folder `COLVAR` + `spib_cv_dt*.dat`, `spib_results/*.pkl`,
`graph_outputs/` (time-lag, FES, and state plots), and `protein_outputs_MD/` representative-frame
PDBs, with `spib_stage_summary.json` and `spib_stage.log`. Set `paths.af2_structures_dir` to the
AF2 run folders to add the static-structure overlay on the free-energy surface; it is skipped
when unset. The analysis logic is `spib_analysis.py`, an environment-variable-driven adaptation
of `prediction/SPIBAnalysis.py`.
