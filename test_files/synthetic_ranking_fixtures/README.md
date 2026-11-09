# Synthetic ranking fixtures

**Every file in this directory is constructed. None of it came from a docking run.**

There are no receptors, no ligands, and no poses behind these numbers. Scores are arithmetic —
`-12.0 + 0.1 * (rank - 1)` — chosen so that the rank structure each test depends on is stated
explicitly in `make_fixtures.py` rather than being an accident of real data. The ligand identifiers
are `SYN####` precisely so they cannot be mistaken for ZINC or ChEMBL entries.

Registered in `STALE_REGISTRY.md` (dockflow-speckit) under "Mock-only fixtures", per constitution
Principle X.

Regenerate from the repository root:

```bash
python test_files/synthetic_ranking_fixtures/make_fixtures.py
```

Output is deterministic — no RNG, no clock, no filesystem ordering — so a rerun reproduces the same
bytes. The T007 baseline digests depend on that.

## What each set imitates

| Directory | Imitates | Shape |
|---|---|---|
| `unlabeled_vina/` | `docking_summary.parquet` from `extract_vina_results.py` | `Ligand`, `Vina_Score`, `File_Path`. **No `Class` column** — this is the P1 unlabeled case. Four receptors over one shared 20-ligand set. |
| `labeled_benchmark/` | `enrichment_ranked.parquet` from `calculate_enrichment.py` | Adds `Class`, `Is_Active`, `Rank`, `Cumulative_*`, `Raw_Ligand`. Same 20 ligands, so labeled and unlabeled runs are directly comparable. Stored as `<receptor>/enrichment_ranked.parquet` directories, which is the form `resolve_parquet()` resolves. |
| `unlabeled_dock6/` | `dock6_summary.parquet` from `extract_dock6_results.py` | One row per **pose**: `Ligand`, `Grid_Score`, `Pose_Index`, `Source_File`, `Mol2_Start_Line`, `Mol2_End_Line`. Three poses per ligand, written worst-first so a collapse that ignores the sort keeps the wrong row. |
| `partial_failure/` | A run where one receptor failed for half the library | `rec_01` has all 20 ligands, `rec_02` has 10. `expected_ligands.txt` holds 22 identifiers — the 20 plus `SYN0021`/`SYN0022`, which appear in **no** receptor table. |

## Why four receptors, and why the ties

`unlabeled_vina/` carries more structure than a first reading of the tasks suggests, because the
determinism obligation (SC-011, FR-014) cannot be demonstrated without it. Two receptors are enough
for the labels-optional and score-column tests; the determinism test uses all four.

Ranks for the four ligands that matter:

| Ligand | rec_01 | rec_02 | rec_03 | rec_04 | Role |
|---|---|---|---|---|---|
| `SYN0011` | 6 | 13 | 20 | 9 | permuted pair with `SYN0012` |
| `SYN0012` | 20 | 13 | 6 | 9 | permuted pair with `SYN0011` |
| `SYN0013` | 10 | 15 | 10 | 11 | exact-tie pair with `SYN0014` |
| `SYN0014` | 10 | 15 | 10 | 11 | exact-tie pair with `SYN0013` |

`SYN0011` and `SYN0012` hold the **same rank multiset in a different order**. Their RRF sums are
equal in exact arithmetic, but `Σ 1/(k + rank)` accumulates in a different association depending on
`--input` order, so one of them shifts by a single ulp and the two swap places. Two receptors could
never show this: `a + b == b + a` exactly in IEEE-754. Three or more can, and four do. This is what
makes sorting the receptor inputs a correctness fix rather than a tidy-up.

`SYN0013` and `SYN0014` hold **identical** ranks everywhere, so their scores are exactly equal in
either direction. No accumulation order separates them; only an explicit `Ligand` tie-break gives
them a stable order.

Verified against the generated files: reversing `--input` order flips `SYN0011`/`SYN0012` and leaves
`SYN0013`/`SYN0014` tied. Both failure modes are live, which is the condition Principle X asks for —
the guard is checked against a defect that actually reproduces, not merely passed on clean input.
