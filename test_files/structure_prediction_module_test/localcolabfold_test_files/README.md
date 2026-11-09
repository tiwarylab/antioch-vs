# LocalColabFold Test Files

This directory contains three LocalColabFold test configurations:

- `example_msa_only.yaml`: generates one deep MSA from the bundled JAK2 query;
  this contacts the configured MSA service.
- `example_predict_from_sequence.yaml`: generates one deep MSA, then creates
  three reduced-MSA prediction ensembles from it.
- `example_precomputed_a3m_rmsa.yaml`: runs the same three reduced-MSA
  predictions from bundled `3LPB.a3m`, without an MSA-server request.

For each file, set `paths.output_dir`, `environment.localcolabfold_cmd`, and,
for a prediction configuration, `localcolabfold.extra_args` `--data` cache
root. The runner stages A3M input in an isolated A3M-only workspace before
launching LocalColabFold; do not place a FASTA file in that generated workspace.

The configured output directory must be new or have empty `runs/rmsa_*`
directories while `execution.overwrite` is `false`. Use a different
`paths.output_dir` for a repeat test; do not enable overwrite unless replacing
an existing run is intentional.

Run a selected configuration from the repository root:

```bash
python structure_prediction_module/prediction/run_localcolabfold.py \
  --config test_files/structure_prediction_module_test/localcolabfold_test_files/example_predict_from_sequence.yaml
```

The three output folders are:

```text
<output_dir>/runs/rmsa_8_16/
<output_dir>/runs/rmsa_16_32/
<output_dir>/runs/rmsa_32_64/
```

To use another alignment, replace `paths.sequence_file` with its A3M path. The
first sequence in every supplied A3M must be the query sequence.

To test `prediction_only` with an A3M produced by `example_msa_only.yaml`, copy
`example_precomputed_a3m_rmsa.yaml`, set `paths.sequence_file` to the generated
deep A3M, and choose a new prediction `paths.output_dir`.
