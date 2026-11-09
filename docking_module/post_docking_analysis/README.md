# Post-Docking Analysis

This folder contains result extractors for ANTIOCH-VS docking runs.

- `extract_docking_results.py` is the master entry point. Set `--software`
  to `vina` or `dock6`, and all remaining flags are forwarded to the matching
  extractor. It passes through the agent-readable JSON summary path to the
  selected extractor.
- `extract_vina_results.py` parses Vina `.pdbqt` outputs and ranks by
  `Vina_Score`. It can write a compact JSON summary for agents.
- `extract_dock6_results.py` parses DOCK6 scored `.mol2` outputs and ranks by
  `Grid_Score`. It can write a compact JSON summary for agents.
- `calculate_enrichment.py` combines active and decoy summary parquet files and
  reports DUD-E-compatible ROC AUC, adjusted LogAUC, and EF1 with explicit
  units and optional conventional cutoff diagnostics.
- `build_unweighted_rrf_consensus.py` fuses one ranking table per receptor into
  an unweighted reciprocal-rank-fusion consensus ranking. Active/decoy labels are
  optional: a benchmark run passes `enrichment_ranked.parquet`, an unlabeled
  screen passes the extractor summary directly.
- `collect_top_consensus_poses.py` traces the top consensus ligands back to the
  selected pose and receptor that produced each one.
- `stage_consensus_ctmd_pairs.py` converts that collection into one CTMD-ready
  receptor-ligand project per selected consensus hit.
- `analyze_prolif_interactions.py` converts top-ranked Vina PDBQT poses from
  an extracted CSV table to SDF and calculates ProLIF protein-ligand
  interaction fingerprints against a prepared receptor PDB.
- `environment.yaml` creates a small environment with `pandas`, `pyarrow`,
  `tqdm`, Open Babel, Meeko, RDKit, MDAnalysis, and ProLIF for table writing,
  optional structure export, enrichment, and interaction analysis.

## Environment

Create the environment without modifying your main docking environment:

```bash
micromamba create -f docking_module/post_docking_analysis/environment.yaml
micromamba activate dockflow-post-docking
```

Check the key tools:

```bash
python -c "import pandas, pyarrow, prolif, MDAnalysis; print('post-docking dependencies OK')"
obabel -V
mk_export.py --help
```

`obabel` is needed for DOCK6 Mol2 to SDF export. `mk_export.py` is needed for
Vina PDBQT to SDF export.

## Consensus To CTMD

Run this workflow separately for one docking engine at a time. Raw DOCK6
`Grid_Score` values and Vina `Vina_Score` values are not commensurate, whereas
the reciprocal-rank-fusion consensus is intentionally score-scale independent.

1. Build a consensus. There are two routes in, and they differ only in what you
   pass to `--input`.

   **Unlabeled screening** — your own library, no actives and decoys. Pass the
   extractor output directly. Enrichment cannot run without labels and is not
   needed: fusion works on rank order, so labels never enter the calculation.

```bash
python build_unweighted_rrf_consensus.py \
  --input receptor_1=/absolute/path/receptor_1/docking_summary.parquet \
  --input receptor_2=/absolute/path/receptor_2/docking_summary.parquet \
  --out_dir /absolute/path/ensemble_consensus
```

   Use `dock6_summary.parquet` for DOCK6. Those tables carry one row per pose;
   consensus keeps each ligand's best-scoring pose.

   **Benchmark / DUD-E** — pass the enrichment output, as before. Behaviour on
   this route is unchanged.

```bash
python build_unweighted_rrf_consensus.py \
  --input receptor_1=/absolute/path/receptor_1/dude_enrichment/enrichment_ranked.parquet \
  --input receptor_2=/absolute/path/receptor_2/dude_enrichment/enrichment_ranked.parquet \
  --out_dir /absolute/path/ensemble_consensus
```

   Required columns are the ligand column and one score column. All receptors in
   a run must resolve to the **same** score column; a mismatch is an error rather
   than a silent coercion, since fusing a Vina affinity against a DOCK6 grid
   score would still compute.

   **Labels are optional, but all-or-nothing.** `Class` absent from every input is
   a normal unlabeled run. What is rejected: a *mixture* of labelled and unlabeled
   inputs (inferring the missing labels would invent data, dropping the present
   ones would discard the benchmark), a `Class` value that is not `active` or
   `decoy`, and a labelled table containing only one class — enrichment over a
   single class is undefined, not merely poor.

   **The output directory must be clean.** A run refuses an `--out_dir` that
   already holds consensus outputs, naming what it found. Pass `--overwrite` to
   replace them; it removes every artifact this stage owns, including ones the
   new run will not write, so an unlabeled rerun cannot leave a labelled
   `consensus_actives.parquet` sitting beside a manifest reading
   `labels_present: false`. Files the stage does not own are never touched.

### What a consensus run writes

| File | When | Notes |
|---|---|---|
| `consensus_all.parquet` | always | every ligand, all consensus columns |
| `consensus_top_hits.csv` | always | the published shortlist, truncated by `--csv_n` |
| `consensus_manifest.csv` | always | one row per receptor: score column, ligands scored, ligands failed |
| `manifest.json` | always | run-level facts (see below) |
| `consensus_actives.parquet`, `consensus_decoys.parquet` | labelled runs only | evaluation inputs for enrichment; **omitted**, not empty, when there are no labels |
| `failures.csv` | only with `--expected_ligands` | inventory members that reached no receptor at all |

### The absence penalty is one value for the whole run

A ligand missing from a receptor is ranked below every ligand that scored, rather
than being rewarded by omission. That penalty is **one shared number per run**:

| Basis | Penalty | When |
|---|---|---|
| `inventory_plus_one` | `len(--expected_ligands) + 1` | an inventory is supplied |
| `max_scored_plus_one` | `max(ligands scored by any receptor) + 1` | otherwise |

Both the value and the basis are recorded in `manifest.json`, because a shortlist
cannot be audited without knowing what an absence cost.

> **Changed 2026-09-09.** The penalty used to be the length of the individual
> receptor table the ligand was absent from. Receptors finish unevenly, so that
> gave a *different* penalty per receptor — and a smaller one where more ligands
> had failed, making the least complete evidence the most forgiving. Runs with
> uneven receptor completion re-rank slightly under the corrected rule. On the
> JAK2 reference set (5 receptors, 6,437–6,441 ligands each), 6 of 6,441 ligands
> moved by one place, each a partially-docked ligand falling below a
> fully-covered one. Top-25, top-100 and top-10% membership were unchanged.

Every ligand row carries `coverage_succeeded` and `coverage_attempted` — how many
receptors it actually docked against, out of how many were in the run. A ligand
sitting at the absence penalty in three of four receptors is a different finding
from one that genuinely ranked mid-table everywhere, and the ranking alone does
not distinguish them.

`manifest.json` records the run-level facts that do not fit a per-receptor table:
`engine`, `method`, `rrf_k`, `score_column`, `score_direction`, `labels_present`,
`tie_break`, `receptor_order`, `expected_ligand_inventory`,
`total_failure_reporting`, `artifacts[]`, and `refinement_ran`. `receptor_order`
is there as evidence, not decoration: fusion sums one float per receptor and
float addition is not associative, so the order the receptors are added in is
what makes the result reproducible.

### Failure reporting is conditional, and says so

Consensus sees only the union of ligands present in the per-receptor tables. **A
ligand that failed in every receptor is absent from all of them and is therefore
invisible here by construction.** It cannot be reported from this data alone.

So `failures.csv` is written **only** when you pass `--expected_ligands` — the
prepared-ligand list from ligprep, or your input SMILES identifiers, one per
line. Without it, no `failures.csv` is written and `manifest.json` records
`total_failure_reporting: "unavailable_no_inventory"`.

**An absent or empty `failures.csv` does not mean nothing failed.** The manifest
field is the authority on whether the question could be answered at all. Partial
failures are always reported either way: they show up as
`Rank_<receptor> = total_in_receptor + 1` and reduced `coverage_succeeded`.

```bash
python build_unweighted_rrf_consensus.py \
  --input receptor_1=/absolute/path/receptor_1/docking_summary.parquet \
  --input receptor_2=/absolute/path/receptor_2/docking_summary.parquet \
  --expected_ligands /absolute/path/ligprep_out/prepared_ligands.txt \
  --out_dir /absolute/path/ensemble_consensus
```

2. Collect the top fraction. `--receptor_template` must point to the prepared
   receptor **PDB** used for each docking receptor. DOCK6 source MOL2 poses are
   extracted exactly; Vina source PDBQT files are reduced to their best model.

   Pass the same per-receptor tables you gave to consensus — the extractor
   summaries for an unlabeled run, the enrichment directories for a benchmark
   one. An unlabeled collection manifest simply has no `Class` column, and the
   pose filenames drop the active/decoy segment.

```bash
python collect_top_consensus_poses.py \
  --consensus_all /absolute/path/ensemble_consensus/consensus_all.parquet \
  --input receptor_1=/absolute/path/receptor_1/dude_enrichment \
  --input receptor_2=/absolute/path/receptor_2/dude_enrichment \
  --receptor_template '/absolute/path/minimized_receptors/{receptor}.pdb' \
  --out_dir /absolute/path/top_1pct_collection \
  --top_fraction 0.01
```

3. Stage each selected pair for CTMD. The stager validates that every DOCK6
   pose is a single MOL2 molecule. For Vina, it converts the collected
   single-model PDBQT pose through Meeko and Open Babel to MOL2. CTMD then
   assigns AM1-BCC charges during its own setup.

```bash
python stage_consensus_ctmd_pairs.py \
  --collection_manifest /absolute/path/top_1pct_collection/collection_manifest.csv \
  --out_dir /absolute/path/ctmd_projects \
  --ctmd_module /absolute/path/antioch-vs/ctmd_module
```

This produces one independent project per receptor-ligand pair:

```text
ctmd_projects/<consensus-rank>_<ligand>_<receptor>/
  rec_prepped.pdb
  molecules/ligand_pose.mol2
  run_ctmd_setup.sh
```

Run `./run_ctmd_setup.sh` inside an activated AMBER environment. Do not merge
projects whose selected ligands came from different receptors: CTMD uses one
receptor PDB per topology setup.

## Master Command

DOCK6 example:

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

Vina example:

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

The master command writes `dockflow_post_docking_summary.json` in `--out_dir`
by default. Set a custom JSON path with:

```bash
--json_out /absolute/path/to/post_docking_summary.json
```

The JSON summary is intentionally compact and is written by the selected
extractor. It records status, input paths, output artifact paths, row/score
metrics, export counts, and a small top-hit preview. The full score table
remains in Parquet.

If you call the extractors directly, their default JSON names are:

- `extract_vina_results.py`: `dockflow_vina_extraction_summary.json`
- `extract_dock6_results.py`: `dockflow_dock6_extraction_summary.json`

## DUD-E-Compatible Enrichment

`calculate_enrichment.py` evaluates a docking ranking of known actives against
DUD-E decoys. It reports the same three concepts used by the DUD-E paper and
website, with units made explicit:

| Reported metric | Meaning | Random ranking |
| --- | --- | --- |
| `ROC_AUC_Percent` | Area under the full ROC curve: active recall versus decoy recall. | `50%` |
| `DUD_E_LogAUC_0.001_Percent` | DUD-E adjusted semilog ROC area over decoy recall `0.001` to `1.0`. | `0%` |
| `EF1_Percent` | Percent of actives recovered at 1% decoy recall. | `1%` |

The script takes active and decoy result folders, ranks all ligands by their
best docking score, and then calculates these metrics:

```bash
python docking_module/post_docking_analysis/calculate_enrichment.py \
  --actives_dir /absolute/path/to/actives_summary \
  --decoys_dir /absolute/path/to/decoys_summary \
  --out_dir /absolute/path/to/enrichment_summary
```

It auto-detects `Grid_Score` for DOCK6 or `Vina_Score` for Vina, assumes lower
scores are better, and collapses to the best pose **separately within each
class** for every ligand by default. Use `--higher_is_better` for a score where
larger is favorable, and `--keep_all_poses` only when pose-level rather than
ligand-level evaluation is intentional.

If ligand preparation created duplicate output names such as
`CHEMBL1087268-again1`, add `--strip_again_suffix`. This preserves the raw
name in `Raw_Ligand`, removes only a terminal `-againN` suffix for ranking, and
keeps the best score for the parent ligand separately in the active and decoy
classes. Use this for DUD-E runs with those duplicated preparation files.

### Mathematics

Let `A` be the number of scored active ligands and `D` the number of scored
decoys. At each score threshold, the script uses:

```text
active recall = recovered actives / A
decoy recall  = recovered decoys / D
```

`ROC_AUC_Percent` is `100 * integral(active recall d(decoy recall))` over the
complete ROC curve. `50%` is random ranking.

`DUD_E_LogAUC_0.001_Percent` follows Mysinger et al. exactly in convention: it
integrates active recall against `log10(decoy recall)` from `0.001` to `1.0`,
normalizes by the three log decades, then subtracts the random semilog area.
For the default lower bound, the unadjusted random baseline is `14.462%`; the
reported adjusted DUD-E LogAUC therefore has random performance at `0%`.

`EF1_Percent` is **not** a conventional fold enrichment. It is
`100 * active recall` when decoy recall is `0.01`: the percent of actives found
once 1% of decoys have been found. This is the DUD-E website's `Enrichment`
column. The default is `--ef_decoy_fraction 0.01`; changing that value changes
the metric label accordingly, for example to `EF2_Percent`.

Equal docking scores are processed as a single score group. This avoids an
arbitrary input-file order affecting AUC, LogAUC, or EF1; values within a tied
group are linearly interpolated on the ROC curve.

### Outputs

The script writes:

- `enrichment_ranked.parquet`: full ranked active/decoy table
- `enrichment_ranked.csv`: ranked table for quick inspection
- `enrichment_metrics.csv`: the three default DUD-E-compatible metrics and their definitions
- `enrichment_roc_curve.csv`: score-tie-aware ROC points used for ROC AUC and EF1
- `enrichment_logauc_curve.csv`: the semilog ROC points used for adjusted DUD-E LogAUC
- `enrichment_summary.json`: the metrics, scoring direction, and input/scored coverage counts

The summary is deliberately explicit. For example, the default keys are
`roc_auc_percent`, `dude_logauc_percent`, and `ef1_percent`; no raw LogAUC or
ambiguous generic `Enrichment_Factor` is reported.

### Comparing Targets

`plot_multi_target_enrichment.py` overlays score-tie-aware ROC curves from
multiple `calculate_enrichment.py` output directories. It plots decoy recall
and active recall as percentages, with each target's **adjusted DUD-E LogAUC**
in the legend:

```bash
python plot_multi_target_enrichment.py \
  --target JAK2=/path/to/jak2/enrichment_output \
  --target ABL1=/path/to/abl1/enrichment_output \
  --target CDK2=/path/to/cdk2/enrichment_output \
  --out multi_target_enrichment.png
```

The input directory must contain `enrichment_roc_curve.csv` and
`enrichment_summary.json`. If a non-default `--out_prefix` was used when
calculating enrichment, pass the same name through `--prefix`.

### Optional Diagnostics

Top-of-list summaries can be useful while choosing an experimental follow-up
set, but they are not DUD-E EF1. They are disabled by default and written to a
separate `enrichment_diagnostics.csv` only when requested:

```bash
--score_column Grid_Score
--ligand_column Ligand
--diagnostic_fractions 0.001,0.01,0.05
--diagnostic_top_counts 100,500,1000
--logauc_min_decoy_fraction 0.001
--keep_all_poses
--higher_is_better
```

Those diagnostics use `Conventional_Fold_Enrichment`, defined as observed
actives divided by the number expected in a random selection of the same number
of **all** ranked ligands. Do not compare that fold value directly with DUD-E
`EF1_Percent`, which uses a 1% **decoy-recall** cutoff and is reported as active
recovery percent.

### Comparing With a DUD-E Website Result

For a like-for-like target comparison, use the DUD-E clustered active set and
its matching decoys, score every ligand, and retain one best score per ligand.
The active and decoy counts in `enrichment_summary.json` should match the
benchmark set after accounting for any ligands that failed preparation or
docking. A different receptor preparation, score, ligand subset, missing scores,
or use of raw rather than clustered actives is a different experiment and should
not be expected to reproduce the DUD-E website values.

## ProLIF Interaction Analysis

After extracting and ranking Vina docking results, calculate protein-ligand
interaction fingerprints for the top hits:

```bash
python docking_module/post_docking_analysis/analyze_prolif_interactions.py \
  --scores_csv /absolute/path/to/vina_top_hits.csv \
  --receptor_pdb /absolute/path/to/prepared_receptor_with_hydrogens.pdb \
  --out_dir /absolute/path/to/prolif_summary \
  --top_n 10000 \
  -c 8
```

The input CSV must include ligand names and docked PDBQT file paths. By
default, the script expects the Vina extractor columns:

- `Ligand`
- `Vina_Score`
- `File_Path`

Use `--ligand_column`, `--score_column`, or `--file_path_column` if your CSV
uses different names. For each selected row, the script keeps only pose 1 from
the Vina PDBQT, converts it to SDF with `mk_export.py`, and runs ProLIF one
ligand at a time through the single-structure generator. The PDBQT-to-SDF
conversion step is parallelized with the same `-c/--cpus` value by default.
Existing SDFs are reused on reruns, so interrupted runs can be resumed without
repeating already converted ligands. By default, interaction cutoffs use the
`schrodinger_maestro` preset where ProLIF exposes equivalent parameters.
Fingerprinting intentionally avoids ProLIF's iterable batch mode because
validation showed it can miss interactions that are found by direct one-ligand
fingerprints for some ligand/residue pairs.
The script also explicitly initializes RDKit ring information before ProLIF
fingerprinting so fused aromatic systems are handled deterministically during
pi interaction matching.

Outputs:

- `prolif_ligand_interaction_counts.parquet`: one row per ligand with counts
  for each ProLIF interaction type plus derived totals such as `HBond_Count`,
  `SaltBridge_Count`, `PiStacking_Count`, and `Total_Interaction_Count`
- `prolif_ligand_interaction_counts.csv`: CSV version for quick inspection
- `prolif_interactions_long.parquet`: long-form ligand/residue/interaction
  records
- `prolif_interactions_long.csv`: CSV version of the long-form table
- `prolif_residue_interaction_frequencies.csv`: residue-level interaction
  frequencies across ligands
- `prolif_converted_ligands.csv`: mapping from input CSV rows to converted SDFs
- `prolif_interaction_summary.json`: compact run summary

Useful options:

```bash
--interaction_preset schrodinger_maestro
--interaction_preset schrodinger_glide
--interaction_preset prolif
--hbond_distance 3.2
--hbond_dha_angle 110 180
--ionic_distance 5.0
--pication_distance 6.6
--pication_angle 0 30
--custom_parameters_json custom_prolif_parameters.json
--interactions Hydrophobic HBDonor HBAcceptor PiStacking Anionic Cationic
--residues all
--residues ASP25.A TYR100.A
--conversion_cpus 8
--overwrite_sdf
--guess_bonds
--keep_pose_pdbqt
```

Use a receptor PDB that is already prepared for interaction analysis and
contains explicit hydrogens. ProLIF can run without hydrogens in some cases,
but hydrogen-bond and ionic interaction calls are more reliable when receptor
and ligand protonation are meaningful.

Schrodinger preset mapping:

- `schrodinger_maestro`: H-bonds use a ProLIF donor-acceptor heavy-atom
  distance approximation of 3.5 A with D-H...A angle 120-180 degrees.
  Schrodinger Maestro defines H-bonds with H...A distance 2.8 A, but ProLIF's
  explicit-hydrogen H-bond classes use donor-acceptor distance instead. Salt
  bridge distance is 5.0 A; pi-pi face-to-face centroid distance is 4.4 A with
  plane angle 0-30 degrees; pi-pi edge-to-face centroid distance is 5.5 A with
  plane angle 30-90 degrees; pi-cation distance is 6.6 A with angle 0-30
  degrees.
- `schrodinger_glide`: same as above except the H-bond D-H...A angle is
  relaxed to 90-180 degrees. Schrodinger Glide defines H-bonds with H...A
  distance 2.5 A, but the ProLIF preset still uses the 3.5 A donor-acceptor
  distance approximation.

Preset values can be overridden from the command line. Common overrides include
`--hbond_distance`, `--hbond_dha_angle`, `--ionic_distance`,
`--pication_distance`, `--pication_angle`, `--pistacking_ftf_distance`,
`--pistacking_ftf_plane_angle`, `--pistacking_etf_distance`,
`--pistacking_etf_plane_angle`, `--hydrophobic_distance`, and
`--vdw_tolerance`. For interaction parameters not covered by these flags, pass
a JSON string or JSON file with `--custom_parameters_json`. The JSON keys should
match ProLIF interaction names and parameter names:

```json
{
  "HBAcceptor": {"distance": 3.2, "DHA_angle": [110, 180]},
  "HBDonor": {"distance": 3.2, "DHA_angle": [110, 180]},
  "PiCation": {"distance": 6.0, "angle": [0, 35]},
  "CationPi": {"distance": 6.0, "angle": [0, 35]}
}
```

Two Schrodinger criteria are not exactly reproduced because ProLIF does not
expose matching parameters in the standard interaction classes: the H...A-B
acceptor angle for explicit-hydrogen H-bonds, and the good/bad/ugly contact
ratio classification for generic contacts. ProLIF `VdWContact` is still
reported, but it uses ProLIF's van der Waals contact definition rather than
Schrodinger's three-bin contact ratio scheme.

## Notes

For DOCK6 multi-Mol2 inputs, multiple scored poses or conformers can share the
same ligand name. The master parquet keeps every scored pose. Use
`--best_per_ligand` when you want the terminal table, CSV, and exported
structures to include only the best-scoring pose per ligand.
The enrichment script also collapses to one best-scoring row per ligand by
default so active/decoy recovery is ligand-based rather than pose-based. This
collapse is done separately within the active and decoy classes, and the summary
JSON reports input/scored row counts and unique ligand counts for both classes.
