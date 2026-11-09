"""Regression checks for multi-receptor consensus ranking.

`build_unweighted_rrf_consensus.py` had no test before this module. The fixtures
it runs against are constructed, not produced by a docking run — see
`test_files/synthetic_ranking_fixtures/README.md` and the entry in
`STALE_REGISTRY.md`. They exercise rank arithmetic, column handling, and
determinism, and say nothing about whether real Vina or DOCK6 output parses.
"""

import hashlib
import io
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from docking_module.post_docking_analysis import (
    build_unweighted_rrf_consensus as consensus,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "docking_module" / "post_docking_analysis" / "build_unweighted_rrf_consensus.py"
FIXTURES = REPO_ROOT / "test_files" / "synthetic_ranking_fixtures"

UNLABELED_VINA = FIXTURES / "unlabeled_vina"
LABELED = FIXTURES / "labeled_benchmark"
DOCK6 = FIXTURES / "unlabeled_dock6"
PARTIAL = FIXTURES / "partial_failure"


def vina_input(receptor: str) -> str:
    return f"{receptor}={UNLABELED_VINA / f'{receptor}_docking_summary.parquet'}"


def run_consensus(*arguments: str, out_dir: Path) -> subprocess.CompletedProcess:
    """Invoke the CLI the way a user does, rather than importing main().

    Argument handling is part of what these tests cover, so the subprocess
    boundary is deliberate.
    """
    command = [sys.executable, str(SCRIPT), *arguments, "--out_dir", str(out_dir)]
    return subprocess.run(command, capture_output=True, text=True, cwd=REPO_ROOT)


def load_table(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)


def projected_digest(frame: pd.DataFrame, columns: list[str]) -> str:
    """SHA-256 of the named columns, in order, as canonical CSV.

    Not a digest of the file's bytes, for two reasons. Adding `coverage_*` and
    `Failed_Ligands` means no labeled artifact can stay byte-identical and still
    carry this slice, so the invariant worth asserting is that every column that
    existed before holds the same values in the same row order. And Parquet
    stores the writing library's version in its schema metadata, so raw bytes
    would report a pandas upgrade as a change in the science.
    """
    return hashlib.sha256(
        frame[columns].to_csv(index=False, lineterminator="\n").encode()
    ).hexdigest()


class ConsensusFixtureCase(unittest.TestCase):
    """Shared scaffolding: a scratch out_dir per test, and fixture loaders."""

    def setUp(self) -> None:
        self._temporary = tempfile.TemporaryDirectory()
        self.out_dir = Path(self._temporary.name) / "consensus"
        self.addCleanup(self._temporary.cleanup)

    def run_ok(self, *arguments: str) -> subprocess.CompletedProcess:
        result = run_consensus(*arguments, out_dir=self.out_dir)
        self.assertEqual(
            result.returncode,
            0,
            msg=f"consensus exited {result.returncode}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        )
        return result

    def manifest(self) -> dict:
        return json.loads((self.out_dir / "manifest.json").read_text())


class UnlabeledPathTests(ConsensusFixtureCase):
    """T009 (FR-006, FR-008): a run with no Class column anywhere."""

    def test_unlabeled_extractor_parquet_reaches_a_shortlist(self):
        self.run_ok("--input", vina_input("rec_01"), "--input", vina_input("rec_02"))

        shortlist = self.out_dir / "consensus_top_hits.csv"
        self.assertTrue(shortlist.is_file(), "the published shortlist must be produced")
        self.assertTrue((self.out_dir / "consensus_all.parquet").is_file())

        frame = pd.read_csv(shortlist)
        self.assertEqual(len(frame), 20)
        self.assertNotIn("Class", frame.columns)

        self.assertFalse(
            (self.out_dir / "consensus_actives.parquet").exists(),
            "label-derived evaluation input must be omitted, not empty",
        )
        self.assertFalse((self.out_dir / "consensus_decoys.parquet").exists())

        self.assertIs(self.manifest()["labels_present"], False)


class MultiPoseTests(ConsensusFixtureCase):
    """T010 (FR-006): DOCK6 tables carry one row per pose."""

    def test_multi_pose_table_collapses_to_best_pose_per_ligand(self):
        frame, score_column = consensus.read_ranking(
            "rec_01",
            DOCK6 / "rec_01_dock6_summary.parquet",
            "Ligand",
            None,
            lower_is_better=True,
        )

        self.assertEqual(score_column, "Grid_Score")
        self.assertEqual(len(frame), 20, "three poses per ligand must collapse to one row")
        self.assertFalse(frame["Ligand"].duplicated().any())

        source = pd.read_parquet(DOCK6 / "rec_01_dock6_summary.parquet")
        best = source.groupby("Ligand")["Grid_Score"].min()
        kept = frame.set_index("Ligand")["Grid_Score"]
        pd.testing.assert_series_equal(
            kept.sort_index(), best.sort_index(), check_names=False
        )


class ScoreColumnTests(ConsensusFixtureCase):
    """T011 (FR-015): receptors must agree on one score column."""

    def test_mismatched_score_columns_are_rejected(self):
        result = run_consensus(
            "--input",
            vina_input("rec_01"),
            "--input",
            f"rec_d6={DOCK6 / 'rec_01_dock6_summary.parquet'}",
            out_dir=self.out_dir,
        )

        self.assertNotEqual(result.returncode, 0, "a silent coercion is the failure mode")
        combined = result.stdout + result.stderr
        self.assertIn("Vina_Score", combined)
        self.assertIn("Grid_Score", combined)


class DeterminismTests(ConsensusFixtureCase):
    """T012 (FR-014, SC-011): reversing --input order changes nothing.

    Four receptors, not two. With two, `a + b == b + a` exactly in IEEE-754 and
    the receptor-ordering defect cannot surface at all. The fixture also holds a
    pair with identical ranks in every receptor, whose RRF scores are exactly
    equal and can only be ordered by the explicit tie-break.
    """

    RECEPTORS = ("rec_01", "rec_02", "rec_03", "rec_04")

    def _shortlist_bytes(self, receptors, out_name: str) -> bytes:
        out_dir = self.out_dir.parent / out_name
        arguments = []
        for receptor in receptors:
            arguments += ["--input", vina_input(receptor)]
        result = run_consensus(*arguments, out_dir=out_dir)
        self.assertEqual(
            result.returncode, 0, msg=f"{result.stdout}\n{result.stderr}"
        )
        return (out_dir / "consensus_top_hits.csv").read_bytes()

    def test_reversed_input_order_yields_byte_identical_shortlist(self):
        forward = self._shortlist_bytes(self.RECEPTORS, "forward")
        reverse = self._shortlist_bytes(tuple(reversed(self.RECEPTORS)), "reverse")

        self.assertEqual(
            forward,
            reverse,
            "consensus_top_hits.csv must not depend on the order --input was supplied in",
        )

    def test_exactly_tied_ligands_hold_a_stable_order(self):
        """SYN0013 and SYN0014 hold identical ranks everywhere, so their RRF
        scores are exactly equal and no accumulation order can separate them.

        Honest bound on this test: it does NOT reproduce the defect it guards.
        numpy's quicksort falls back to insertion sort below ~16 elements, and
        this fixture has 20 rows, so the unstable sort happens to preserve input
        order here — the assertion passed with the tie-break reverted. It is
        regression protection for a structural guarantee, not a defect shown
        failing first. Reproducing it would need a fixture large enough to reach
        a real quicksort partition. Recorded rather than papered over, per
        Principle X.
        """
        for label, receptors in (
            ("tie_forward", self.RECEPTORS),
            ("tie_reverse", tuple(reversed(self.RECEPTORS))),
        ):
            with self.subTest(input_order=label):
                frame = pd.read_csv(io.BytesIO(self._shortlist_bytes(receptors, label)))
                order = [
                    ligand
                    for ligand in frame["Ligand"].tolist()
                    if ligand in {"SYN0013", "SYN0014"}
                ]
                self.assertEqual(
                    order,
                    ["SYN0013", "SYN0014"],
                    "exactly-tied ligands must fall in Ligand order, not arbitrary order",
                )


class CoverageTests(ConsensusFixtureCase):
    """T013 (FR-006): a receptor that failed for half the library."""

    def _run_partial(self, *extra: str):
        return self.run_ok(
            "--input",
            f"rec_01={PARTIAL / 'rec_01_docking_summary.parquet'}",
            "--input",
            f"rec_02={PARTIAL / 'rec_02_docking_summary.parquet'}",
            *extra,
        )

    def test_absent_ligands_are_penalised_and_retained(self):
        self._run_partial()
        frame = pd.read_csv(self.out_dir / "consensus_top_hits.csv").set_index("Ligand")

        self.assertEqual(len(frame), 20, "partially-docked ligands stay in the consensus")

        docked = set(
            pd.read_parquet(PARTIAL / "rec_02_docking_summary.parquet")["Ligand"]
        )
        missing = sorted(set(frame.index) - docked)
        self.assertEqual(len(missing), 10)

        # The shared penalty, not this receptor's own length. rec_01 scored 20
        # and rec_02 scored 10, so the superseded per-table rule gave 11 here —
        # a lighter penalty precisely because rec_02 failed more often. This
        # assertion changed with the correction and is meant to.
        penalty = 21
        for ligand in missing:
            self.assertEqual(
                frame.loc[ligand, "Rank_rec_02"],
                penalty,
                "a docking failure ranks one worse than every ligand that scored "
                "anywhere in the run, on the same denominator for every receptor",
            )

    def test_coverage_columns_count_receptors_reached(self):
        self._run_partial()
        frame = pd.read_csv(self.out_dir / "consensus_top_hits.csv").set_index("Ligand")

        self.assertIn("coverage_succeeded", frame.columns)
        self.assertIn("coverage_attempted", frame.columns)
        self.assertTrue((frame["coverage_attempted"] == 2).all())

        docked = set(
            pd.read_parquet(PARTIAL / "rec_02_docking_summary.parquet")["Ligand"]
        )
        for ligand in frame.index:
            expected = 2 if ligand in docked else 1
            self.assertEqual(frame.loc[ligand, "coverage_succeeded"], expected)

    def test_per_receptor_manifest_reports_failure_counts(self):
        self._run_partial()
        manifest = pd.read_csv(self.out_dir / "consensus_manifest.csv").set_index(
            "Receptor"
        )

        self.assertEqual(list(manifest.columns[:3]), ["Score_Column", "Scored_Ligands", "Failed_Ligands"])
        self.assertEqual(manifest.loc["rec_01", "Failed_Ligands"], 0)
        self.assertEqual(manifest.loc["rec_02", "Failed_Ligands"], 10)


class FailureReportingTests(ConsensusFixtureCase):
    """T014 (FR-006): the difference between "none failed" and "cannot tell"."""

    INPUTS = None

    def setUp(self) -> None:
        super().setUp()
        self.INPUTS = [
            "--input",
            f"rec_01={PARTIAL / 'rec_01_docking_summary.parquet'}",
            "--input",
            f"rec_02={PARTIAL / 'rec_02_docking_summary.parquet'}",
        ]

    def test_without_an_inventory_the_gap_is_declared_not_implied(self):
        self.run_ok(*self.INPUTS)

        self.assertFalse(
            (self.out_dir / "failures.csv").exists(),
            "no inventory means the question cannot be answered; an empty file would imply it was",
        )
        self.assertEqual(
            self.manifest()["total_failure_reporting"], "unavailable_no_inventory"
        )
        self.assertIsNone(self.manifest()["expected_ligand_inventory"])

    def test_an_absent_failures_csv_is_never_a_clean_bill_of_health(self):
        self.run_ok(*self.INPUTS)

        manifest = self.manifest()
        self.assertFalse((self.out_dir / "failures.csv").exists())
        self.assertNotEqual(
            manifest["total_failure_reporting"],
            "none",
            "absence of the file must not be recorded as absence of failures",
        )

    def test_with_an_inventory_zero_coverage_members_are_listed(self):
        inventory = PARTIAL / "expected_ligands.txt"
        self.run_ok(*self.INPUTS, "--expected_ligands", str(inventory))

        failures_path = self.out_dir / "failures.csv"
        self.assertTrue(failures_path.is_file())

        failures = pd.read_csv(failures_path)
        self.assertEqual(
            sorted(failures["Ligand"].tolist()),
            ["SYN0021", "SYN0022"],
            "only inventory members invisible to every receptor are total failures",
        )

        manifest = self.manifest()
        self.assertEqual(manifest["total_failure_reporting"], "complete")
        self.assertEqual(manifest["expected_ligand_inventory"], str(inventory))


class ManifestTests(ConsensusFixtureCase):
    """T024 (FR-008, SC-012 manifest half): the run-level record."""

    def test_manifest_records_the_run_level_facts(self):
        self.run_ok("--input", vina_input("rec_02"), "--input", vina_input("rec_01"))
        manifest = self.manifest()

        self.assertEqual(manifest["method"], "reciprocal_rank_fusion")
        self.assertEqual(manifest["rrf_k"], 60.0)
        self.assertEqual(manifest["score_column"], "Vina_Score")
        self.assertEqual(manifest["score_direction"], "lower_is_better")
        self.assertEqual(manifest["tie_break"], ["Ligand"])
        self.assertIs(manifest["refinement_ran"], False)

        self.assertEqual(
            manifest["receptor_order"],
            ["rec_01", "rec_02"],
            "receptor_order is the evidence for SC-011, so it records the order actually used",
        )

        artifacts = {entry["path"]: entry for entry in manifest["artifacts"]}
        self.assertIn("consensus_top_hits.csv", artifacts)
        shortlist = artifacts["consensus_top_hits.csv"]
        self.assertEqual(shortlist["ordered_by"], "Consensus_RRF_Score")
        self.assertEqual(shortlist["direction"], "higher_is_better")
        self.assertEqual(shortlist["row_count"], 20)

    def test_higher_is_better_is_recorded_not_inferred(self):
        self.run_ok(
            "--input", vina_input("rec_01"), "--input", vina_input("rec_02"),
            "--higher_is_better",
        )
        self.assertEqual(self.manifest()["score_direction"], "higher_is_better")


UNEVEN = FIXTURES / "uneven_completion"


def uneven_input(receptor: str) -> str:
    return f"{receptor}={UNEVEN / f'{receptor}_docking_summary.parquet'}"


class SharedAbsencePenaltyTests(ConsensusFixtureCase):
    """T034, T035 (FR-006, SC-013): one penalty for the run, not one per table.

    The fixture's three receptors scored 12, 9, and 7 ligands. Under the
    superseded rule the penalty was the length of the table the ligand was
    absent from — 13, 10, 8 — so absence from the receptor that failed most
    cost the least. The three numbers differ from the shared value on purpose:
    if the fixture were even, both rules would agree and this would prove
    nothing.
    """

    RECEPTORS = ("unv_a", "unv_b", "unv_c")
    ABSENT_FROM_TWO = "SYN0007"

    def _run(self, *extra: str):
        arguments = []
        for receptor in self.RECEPTORS:
            arguments += ["--input", uneven_input(receptor)]
        return self.run_ok(*arguments, *extra)

    def test_absent_ligand_carries_the_same_penalty_in_every_receptor(self):
        self._run()
        frame = pd.read_csv(self.out_dir / "consensus_top_hits.csv").set_index("Ligand")

        penalties = {
            receptor: frame.loc[self.ABSENT_FROM_TWO, f"Rank_{receptor}"]
            for receptor in ("unv_b", "unv_c")
        }
        self.assertEqual(
            penalties["unv_b"],
            penalties["unv_c"],
            "the same failure must cost the same in every receptor; per-table "
            f"penalties would give 10 and 8 here, got {penalties}",
        )

    def test_penalty_is_max_scored_plus_one_without_an_inventory(self):
        self._run()
        frame = pd.read_csv(self.out_dir / "consensus_top_hits.csv").set_index("Ligand")

        # 12 is the largest per-receptor count in the fixture.
        self.assertEqual(frame.loc[self.ABSENT_FROM_TWO, "Rank_unv_b"], 13)
        self.assertEqual(frame.loc[self.ABSENT_FROM_TWO, "Rank_unv_c"], 13)

        manifest = self.manifest()
        self.assertEqual(manifest["absence_penalty"], 13)
        self.assertEqual(manifest["absence_penalty_basis"], "max_scored_plus_one")

    def test_an_inventory_moves_the_penalty_and_the_basis_says_so(self):
        inventory = UNEVEN / "expected_ligands.txt"
        self._run("--expected_ligands", str(inventory))
        frame = pd.read_csv(self.out_dir / "consensus_top_hits.csv").set_index("Ligand")

        # 15 identifiers in the inventory, deliberately more than any receptor saw.
        self.assertEqual(frame.loc[self.ABSENT_FROM_TWO, "Rank_unv_b"], 16)
        self.assertEqual(frame.loc[self.ABSENT_FROM_TWO, "Rank_unv_c"], 16)

        manifest = self.manifest()
        self.assertEqual(manifest["absence_penalty"], 16)
        self.assertEqual(
            manifest["absence_penalty_basis"],
            "inventory_plus_one",
            "the basis distinguishes counting what the run intended to dock "
            "from counting only what survived",
        )

    def test_coverage_still_reflects_receptors_actually_reached(self):
        self._run()
        frame = pd.read_csv(self.out_dir / "consensus_top_hits.csv").set_index("Ligand")
        self.assertEqual(frame.loc[self.ABSENT_FROM_TWO, "coverage_succeeded"], 1)
        self.assertEqual(frame.loc[self.ABSENT_FROM_TWO, "coverage_attempted"], 3)


class OutputDirectoryTests(ConsensusFixtureCase):
    """T036, T037 (FR-031, SC-014): a directory holds one run's artifacts."""

    OWNED = (
        "consensus_all.parquet",
        "consensus_top_hits.csv",
        "consensus_actives.parquet",
        "consensus_decoys.parquet",
        "consensus_manifest.csv",
        "manifest.json",
        "failures.csv",
    )

    def labeled_inputs(self) -> list[str]:
        return [
            "--input", f"rec_01={LABELED / 'rec_01'}",
            "--input", f"rec_02={LABELED / 'rec_02'}",
        ]

    def unlabeled_inputs(self) -> list[str]:
        return ["--input", vina_input("rec_01"), "--input", vina_input("rec_02")]

    def test_a_second_run_into_a_used_directory_is_refused(self):
        self.run_ok(*self.unlabeled_inputs())

        result = run_consensus(*self.unlabeled_inputs(), out_dir=self.out_dir)
        self.assertNotEqual(
            result.returncode, 0, "silently overwriting is the failure mode"
        )
        combined = result.stdout + result.stderr
        self.assertIn("consensus_all.parquet", combined, "the message must name what it found")
        self.assertIn("--overwrite", combined, "and say how to proceed deliberately")

    def test_overwrite_clears_label_artifacts_a_later_unlabeled_run_will_not_write(self):
        self.run_ok(*self.labeled_inputs())
        self.assertTrue((self.out_dir / "consensus_actives.parquet").is_file())

        self.run_ok(*self.unlabeled_inputs(), "--overwrite")

        self.assertFalse(
            (self.out_dir / "consensus_actives.parquet").exists(),
            "a labelled artifact surviving an unlabeled rerun makes the directory "
            "assert two contradictory things about one run",
        )
        self.assertFalse((self.out_dir / "consensus_decoys.parquet").exists())
        self.assertIs(self.manifest()["labels_present"], False)

    def test_overwrite_clears_a_failures_csv_the_new_run_will_not_write(self):
        self.run_ok(
            *self.unlabeled_inputs(),
            "--expected_ligands", str(PARTIAL / "expected_ligands.txt"),
        )
        self.assertTrue((self.out_dir / "failures.csv").is_file())

        self.run_ok(*self.unlabeled_inputs(), "--overwrite")

        self.assertFalse(
            (self.out_dir / "failures.csv").exists(),
            "a stale failures.csv would be read as this run's failure report",
        )
        self.assertEqual(
            self.manifest()["total_failure_reporting"], "unavailable_no_inventory"
        )

    def test_unrelated_files_are_left_alone(self):
        self.out_dir.mkdir(parents=True, exist_ok=True)
        bystander = self.out_dir / "notes.md"
        bystander.write_text("not ours\n")

        self.run_ok(*self.unlabeled_inputs())
        self.assertTrue(bystander.is_file(), "--out_dir is not the stage's to clean")
        self.assertEqual(bystander.read_text(), "not ours\n")

    def test_an_owned_name_existing_as_a_directory_fails_clearly(self):
        (self.out_dir / "consensus_all.parquet").mkdir(parents=True)

        result = run_consensus(
            *self.unlabeled_inputs(), "--overwrite", out_dir=self.out_dir
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("consensus_all.parquet", result.stdout + result.stderr)


class LabelIntegrityTests(ConsensusFixtureCase):
    """T038, T039, T040 (FR-032): labels are all-or-nothing, and complete."""

    def _write(self, name: str, frame: pd.DataFrame) -> str:
        path = self.out_dir.parent / "label_inputs" / f"{name}.parquet"
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_parquet(path, index=False)
        return f"{name}={path}"

    @staticmethod
    def _table(classes, first: int = 1) -> pd.DataFrame:
        """`first` shifts the ligand IDs.

        Tables that must not overlap use disjoint IDs on purpose. Reusing the
        same ligands with different labels trips the existing cross-receptor
        class-conflict check, and a test that fails there has not exercised the
        rule it names.
        """
        count = len(classes)
        return pd.DataFrame(
            {
                "Ligand": [f"SYN{index:04d}" for index in range(first, first + count)],
                "Vina_Score": [-10.0 + 0.1 * index for index in range(count)],
                "Class": classes,
            }
        )

    def test_mixed_labeled_and_unlabeled_inputs_are_rejected(self):
        labeled = self._write("lab", self._table(["active", "decoy", "decoy", "active"]))
        unlabeled_frame = self._table(["active", "decoy", "decoy", "active"]).drop(
            columns=["Class"]
        )
        unlabeled = self._write("unl", unlabeled_frame)

        result = run_consensus(
            "--input", labeled, "--input", unlabeled, out_dir=self.out_dir
        )
        self.assertNotEqual(result.returncode, 0)
        combined = result.stdout + result.stderr
        self.assertIn("lab", combined, "the message names which tables carry labels")
        self.assertIn("unl", combined, "and which do not")

    def test_a_class_value_outside_the_two_recognised_ones_is_rejected(self):
        bad = self._write("bad", self._table(["active", "decoy", "unknown", "decoy"]))
        other = self._write("ok", self._table(["active", "decoy", "decoy", "active"]))

        result = run_consensus("--input", bad, "--input", other, out_dir=self.out_dir)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "unknown",
            result.stdout + result.stderr,
            "the offending value must be named, not just counted",
        )

    def test_a_single_class_labeled_table_is_rejected(self):
        # Disjoint ligand IDs, so the only thing that can reject this run is the
        # both-classes requirement — not the class-conflict check.
        one_class = self._write("solo", self._table(["active", "active", "active"]))
        other = self._write("both", self._table(["active", "decoy", "decoy"], first=100))

        result = run_consensus(
            "--input", one_class, "--input", other, out_dir=self.out_dir
        )
        self.assertNotEqual(
            result.returncode,
            0,
            "enrichment over one class is undefined, not merely poor",
        )
        self.assertIn("solo", result.stdout + result.stderr)


class BenchmarkInvarianceTests(ConsensusFixtureCase):
    """T026, T027 (FR-006, Principle VIII): the labeled path must not drift."""

    BASELINE = LABELED / "baseline"

    def labeled_inputs(self) -> list[str]:
        return [
            "--input", f"rec_01={LABELED / 'rec_01'}",
            "--input", f"rec_02={LABELED / 'rec_02'}",
        ]

    def test_labeled_run_reproduces_the_pre_change_baseline(self):
        self.run_ok(*self.labeled_inputs())
        expected = json.loads((self.BASELINE / "digests.json").read_text())

        for name, record in sorted(expected["artifacts"].items()):
            with self.subTest(artifact=name):
                produced = self.out_dir / name
                self.assertTrue(
                    produced.is_file(), f"{name} disappeared from the labeled path"
                )
                frame = load_table(produced)

                missing = [c for c in record["columns"] if c not in frame.columns]
                self.assertEqual(
                    missing, [], f"{name} dropped pre-existing column(s) {missing}"
                )
                self.assertEqual(
                    projected_digest(frame, record["columns"]),
                    record["projected_sha256"],
                    f"{name} changed against the baseline captured before any edit",
                )

    def test_new_columns_are_additive_and_the_baseline_ones_keep_their_order(self):
        self.run_ok(*self.labeled_inputs())
        expected = json.loads((self.BASELINE / "digests.json").read_text())

        for name, record in sorted(expected["artifacts"].items()):
            with self.subTest(artifact=name):
                columns = list(load_table(self.out_dir / name).columns)
                kept = [c for c in columns if c in record["columns"]]
                self.assertEqual(
                    kept,
                    record["columns"],
                    f"{name} reordered or removed columns that existed before this slice",
                )

    def test_labeled_run_still_writes_label_derived_artifacts(self):
        self.run_ok(*self.labeled_inputs())

        self.assertTrue((self.out_dir / "consensus_actives.parquet").is_file())
        self.assertTrue((self.out_dir / "consensus_decoys.parquet").is_file())

        actives = pd.read_parquet(self.out_dir / "consensus_actives.parquet")
        decoys = pd.read_parquet(self.out_dir / "consensus_decoys.parquet")
        self.assertEqual(len(actives), 6)
        self.assertEqual(len(decoys), 14)

        self.assertIs(self.manifest()["labels_present"], True)


if __name__ == "__main__":
    unittest.main()
