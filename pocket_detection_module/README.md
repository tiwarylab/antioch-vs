# Pocket Detection Module

This module adapts external pocket-prediction methods to ANTIOCH-VS docking
inputs. It supports Graph Attention Site Prediction (GrASP) and P2Rank through
one configuration file: `pocket_detection_config.yaml`.
GrASP predicts a binding-site probability for each protein atom and its
`color_pdb.py` step stores those probabilities in the PDB B-factor column.
P2Rank directly predicts and ranks pockets from protein structures.

ANTIOCH-VS does not vendor or install GrASP. Install its separately maintained,
legacy runtime from the [GrASP repository](https://github.com/tiwarylab/GrASP).
Its documented production workflow parses PDB/MOL2 inputs, runs pretrained
inference, and colors PDB files with its predictions.

## Environments

Three separate environments are involved. Only the first is ANTIOCH-VS's.

| Environment | Contains | Declared by |
|-------------|----------|-------------|
| `pocket_detection_env` | the ANTIOCH-VS runners (stdlib + PyYAML; the GrASP pocket postprocessor also needs numpy, scikit-learn, scipy) | `pocket_detection_env.yml` |
| P2Rank's | a JDK, to run the `prank` launcher | you |
| GrASP's | Python 3.7.10, PyTorch, PyTorch Geometric | you, following GrASP's README |

Run the ANTIOCH-VS runners from `pocket_detection_env`, never from GrASP's. They
invoke GrASP's interpreter by absolute path, so they never need to be inside it:

```bash
micromamba env create -f pocket_detection_module/pocket_detection_env.yml
micromamba activate pocket_detection_env
```

### Versions in use

These are the versions this module has actually been run against. They are
recorded so a result can be traced, not because other versions are forbidden.

| Component | Version |
|-----------|---------|
| P2Rank | 2.5.1 (MIT) |
| JDK for P2Rank | OpenJDK 17.0.18 - upstream requires Java 17 or later |
| GrASP | commit `0a01251`, 2026-03-02 (MIT) - the repository has no tags |
| GrASP Python | 3.7.10 |
| GrASP PyTorch | 1.13.1 |
| GrASP PyTorch Geometric | 2.3.1 |

Note the JDK sits exactly at P2Rank's minimum. There is no margin for a
downgrade. GrASP's own README pins NetworkX and OpenBabel but not PyTorch
Geometric, so the two PyTorch versions above are more specific than upstream
documents and cannot be recovered from it.

### P2Rank

P2Rank needs no installation - unpack a release and point at its launcher.
Download the 2.5.1 archive from the
[releases page](https://github.com/rdk/p2rank/releases) and unpack it, then
supply a JDK. A dedicated environment keeps it off the system Java:

```bash
micromamba create -n p2rank openjdk=17
micromamba activate p2rank
```

Set `p2rank.executable` in `pocket_detection_config.yaml` to the unpacked
`prank` launcher, for example `<...>/P2Rank-2.5.1/prank`. Confirm it runs:

```bash
"<...>/P2Rank-2.5.1/prank" -v
```

### GrASP

GrASP is a separately maintained legacy runtime. Clone it and build its
environment following its own README - the commands below are reproduced from
[the GrASP repository](https://github.com/tiwarylab/GrASP) and are theirs, not
ours:

```bash
git clone https://github.com/tiwarylab/GrASP
mamba create -n grasp python==3.7.10
mamba install conda-forge::cython openbabel=2.4.1 rdkit mdtraj mdanalysis
pip install networkx==2.5
```

PyTorch and PyTorch Geometric are additionally required and are not pinned by
those instructions; see the versions table above for what this module was run
against. GrASP also compiles Cython extensions on first use through
`pyximport`, so the first run is slower than later ones.

Then set three keys in `pocket_detection_config.yaml`:

- `grasp.root` - the clone
- `grasp.python` - the Python inside GrASP's environment
- `grasp.color_python` - the Python used for the `color_pdb.py` step, usually
  the same interpreter

### Keeping them working is yours

ANTIOCH-VS does not vendor, install, update, or test P2Rank or GrASP. It calls
them by path and reads what they produce. If either stops working - a JDK
upgrade past what P2Rank accepts, a GrASP dependency that no longer resolves on
Python 3.7.10, a CUDA or PyTorch mismatch - that is a problem in your
environment, not in this module, and ANTIOCH-VS cannot detect or repair it.

What ANTIOCH-VS does guarantee is that failures are recorded, but read the exit
codes carefully - they are not uniform:

- **P2Rank, some receptors failed:** exit code is **0**. The failures are
  isolated deliberately so one bad structure does not discard a whole batch. A
  `WARNING` is printed, and each failure is a `status=FAILED` row in
  `p2rank_manifest.csv` with an error message.
- **P2Rank, every receptor failed:** raises, so the exit code is non-zero.
- **GrASP:** a failing GrASP subprocess propagates immediately, so the exit code
  is non-zero. Colored PDBs that never appear are recorded as
  `color_status=MISSING` in the manifest.

The consequence for anything scripted around this module: **checking the exit
code alone will accept a partially failed P2Rank run.** Check
`n_failed_receptors` in the run JSON, or the manifest rows, before treating the
output as complete.

## Unified Runner

Set `method` to either `"grasp"` or `"p2rank"` in
`pocket_detection_config.yaml`, populate only the matching method section, and
run:

```bash
python pocket_detection_module/run_pocket_detection.py \
  --config pocket_detection_module/pocket_detection_config.yaml
```

Both methods use `paths.input_structures` and `paths.output_dir`. GrASP accepts
one PDB or a directory of PDBs; P2Rank also accepts mmCIF and BinaryCIF input.
The remaining sections are method-specific, so unused settings can stay as
their template values.

For every receptor with at least one detected pocket, both methods write a
rank-1 docking-box definition under `docking_boxes/<receptor>_best_pocket.yaml`
and index them in `docking_boxes/best_pockets.csv`. Each definition contains
the predicted `box.center` and `box.size` as numeric coordinates and as the
space-separated strings expected by ANTIOCH-VS's Vina YAML `box` section. The
same center and size are also provided under `dock6.box` for generating a
DOCK6 grid box.

## GrASP

Set `paths.input_structures` to an AF2RAVE output such as:

```text
/path/to/af2rave_output/cluster_center_structures
```

Set `grasp.root` to an existing GrASP clone, `grasp.python` to the Python in its
inference environment, and `grasp.color_python` to the Python environment that
can run GrASP's `color_pdb.py`.

The runner creates one private workspace per input PDB under
`grasp_workspaces/`, links the external GrASP source and pretrained models into
each workspace, and stages the receptor into GrASP's required production input
path. It then runs `parse_files.py production`, `infer_test_set.py`, and
`color_pdb.py` separately for each receptor. The colored PDBs are copied to
`grasp_colored_pdbs/` and immediately exported as all ranked pockets under
`pockets/`.

GrASP names its colored per-atom-probability artifact `<receptor>_probs.pdb`.
ANTIOCH-VS preserves that suffix only for the colored-PDB artifact; ranked pocket
tables, pocket directories, and `docking_boxes/<receptor>_best_pocket.yaml`
strip it so their receptor IDs match the original input names and P2Rank
exports.

An individual receptor that fails GrASP does not discard the successful
receptors. `grasp_input_manifest.csv` records a `FOUND` or `MISSING` status and
failure reason for every staged input; `grasp_run.json` reports the successful
and failed receptor counts. The pocket postprocessor similarly writes
`pockets/grasp_postprocess_manifest.csv` if a colored PDB cannot be processed.

Before parsing, ANTIOCH-VS compiles GrASP's Cython distance extension once in the
main process. This avoids an upstream race in which `parse_files.py` launches
multiple workers that otherwise try to compile the same extension concurrently.
The compiled extension cache is kept inside the private workspace, preventing
reuse of a binary built by a different Python version.

Do not point two concurrent runs to the same `output_dir`. The external GrASP
clone itself can be shared because ANTIOCH-VS writes inputs and outputs only in
the private workspace.

## P2Rank

P2Rank is a standalone Java application distributed by its
[official repository](https://github.com/rdk/p2rank). Download and unpack a
release separately, then point `p2rank.executable` to its `prank` launcher.
P2Rank requires Java 17 or newer.

Set `paths.input_structures` to one receptor PDB, mmCIF, or BinaryCIF file, or
a directory such as AF2RAVE's
`cluster_center_structures`. For AF2 or AF2RAVE structures, retain
`p2rank.profile: "alphafold"`; it uses P2Rank's B-factor-independent model.
Set it to `null` for conventional X-ray crystal structures.

The runner invokes P2Rank separately for each receptor so every raw prediction
and visualization remains attributable to its input structure. It preserves all
P2Rank output under `p2rank_raw/`, including each `*_predictions.csv`, and
writes these ANTIOCH-VS files:

- `pockets/p2rank_ranked_pockets.csv`: every P2Rank pocket, ranked within each
  receptor. `center_x`, `center_y`, and `center_z` are the P2Rank pocket center;
  `size_x`, `size_y`, and `size_z` come from `docking_box.size`; when it is
  `null`, P2Rank uses its default 25 A cube.
- `p2rank_input_manifest.csv`: receptor-to-raw-output mapping, including a
  `SUCCESS` or `FAILED` status and any error message.
- `p2rank_run.json`: resolved run settings, pocket counts, and successful and
  failed receptor counts.
- `docking_boxes/<receptor>_best_pocket.yaml` and
  `docking_boxes/best_pockets.csv`: method-independent rank-1 docking boxes.

P2Rank failures are isolated per receptor: successful receptors are exported
even when another input fails. Review the manifest and `p2rank_runner.log`
before using the resulting pocket ranking.

`p2rank.visualizations: true` retains P2Rank's generated PyMOL/ChimeraX
visualization files under the raw output directories. The runner does not pick
a single pocket: all ranked pockets are exported for downstream selection.

## Export All GrASP Pockets

Run GrASP first. Then provide its colored PDB files to
`postprocess_grasp_pockets.py`:

```bash
python pocket_detection_module/postprocess_grasp_pockets.py \
  --config pocket_detection_module/grasp_config.yaml
```

The script processes every scored receptor PDB and deliberately exports every
detected pocket. It does not choose a single docking site.

**Pocket boundaries, ranking, and centers come from GrASP's own postprocessing,
not a ANTIOCH-VS reimplementation.** Five functions (`sort_clusters`,
`cluster_atoms_average`, `hull_center`, `get_centroid`, `hulls_from_clusters`)
are vendored verbatim from GrASP's `site_metrics.py`
([tiwarylab/GrASP](https://github.com/tiwarylab/GrASP), MIT license; Smith et
al., *J. Chem. Inf. Model.* 2024) under a clearly marked attribution block in
`postprocess_grasp_pockets.py`. ANTIOCH-VS adds only the docking-support layer
around them (reading the colored PDB, writing one atom-selection PDB per site,
and the padded Vina/DOCK6 box). Running this postprocessor therefore reproduces
GrASP's published pockets exactly — no need to run GrASP's benchmark harness.

For each receptor, the output directory contains:

- `<receptor>/pocket_01.pdb`, `<receptor>/pocket_02.pdb`, and so on: the
  high-scoring atoms assigned to each spatial pocket.
- `grasp_ranked_pockets.csv`: all pockets ranked within each receptor by
  GrASP's site score (the sum of squared GrASP atom probabilities, `score_type:
  square`). `center_x`, `center_y`, and `center_z` are GrASP's site center (the
  convex-hull center of the cluster), directly usable as a Vina box center. By
  default, `size_x`, `size_y`, and `size_z` are derived from the pocket bounds
  with padding and a minimum size; set `docking_box.size` to three values to
  override them.
- `grasp_pocket_run.json`: parameters and output provenance.
- `grasp_postprocess_manifest.csv`: per-colored-PDB success or failure status.
- `../docking_boxes/<receptor>_best_pocket.yaml` and
  `../docking_boxes/best_pockets.csv`: method-independent rank-1 docking boxes.

Pockets are formed by GrASP's recipe: atoms with B-factor probability above
`score_threshold` (GrASP default 0.3) are grouped with average-linkage
clustering using `distance_threshold` Angstrom (GrASP default 15), ranked by
`score_type` (default `square`), and each site's center is its convex-hull
center. `n_sites` caps how many ranked sites are kept. `min_atoms` is a
ANTIOCH-VS-only filter that discards sites with fewer atoms (default 1 = keep all,
matching GrASP). All are exposed in `grasp_config.yaml` and should be evaluated
for a target before use in prospective screening.

`docking_box.size: null` retains GrASP's derived dimensions. Set it to a fixed
cube or rectangular box, for example `size: [25.0, 25.0, 25.0]`, to use those
dimensions for every GrASP pocket. The center remains the detected pocket
center in either mode.

## Runtime Contract

The postprocessor reads only fixed-width PDB records and does not require the
GrASP *environment* (it runs from `pocket_detection_env`), but it does import
the vendored GrASP functions, so it requires `numpy`, `scikit-learn`, and
`scipy`. PyYAML is required only when using `--config`; command-line arguments
work without it. The input PDBs must be the GrASP-colored PDB outputs, not the
original unscored receptor PDBs. A B-factor value of zero is treated as a zero
prediction score.
