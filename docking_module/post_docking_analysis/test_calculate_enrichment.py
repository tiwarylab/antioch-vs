"""Regression checks for DUD-E-compatible enrichment calculations."""

import unittest

import pandas as pd

from docking_module.post_docking_analysis import calculate_enrichment as enrichment


class EnrichmentMetricTests(unittest.TestCase):
    def test_perfect_ranking_reports_perfect_auc_and_ef1(self):
        frame = pd.DataFrame(
            {
                "Score": [0.0, 0.1] + list(range(1, 101)),
                "Ligand": ["a1", "a2"] + [f"d{index}" for index in range(100)],
                "Is_Active": [True, True] + [False] * 100,
            }
        )
        ranked = enrichment.rank_table(frame, "Score", "Ligand", lower_is_better=True, best_per_ligand=True)
        metrics, summary, log_curve, _, _ = enrichment.calculate_metrics(
            ranked,
            fractions=[],
            top_counts=[],
            score_column="Score",
            lower_is_better=True,
            ef_decoy_fraction=0.01,
            logauc_min_decoy_fraction=0.001,
        )

        values = metrics.set_index("Metric")["Value"]
        self.assertEqual(values["ROC_AUC_Percent"], 100.0)
        self.assertEqual(values["EF1_Percent"], 100.0)
        self.assertEqual(log_curve.iloc[0]["Active_Recall"], 1.0)
        _, random_baseline, _ = enrichment.calculate_log_auc(enrichment.build_roc_curve(ranked, "Score"))
        self.assertAlmostEqual(random_baseline, 14.462006247378284)

    def test_equal_scores_are_chance_level_for_auc(self):
        frame = pd.DataFrame(
            {
                "Score": [1.0, 1.0],
                "Ligand": ["active", "decoy"],
                "Is_Active": [True, False],
            }
        )
        ranked = enrichment.rank_table(frame, "Score", "Ligand", lower_is_better=True, best_per_ligand=True)
        roc_curve = enrichment.build_roc_curve(ranked, "Score")

        self.assertAlmostEqual(enrichment.calculate_auc(roc_curve), 0.5)

    def test_custom_ef_cutoff_is_labeled_by_decoy_percentage(self):
        frame = pd.DataFrame(
            {
                "Score": [0.0, 1.0, 2.0],
                "Ligand": ["active", "decoy1", "decoy2"],
                "Is_Active": [True, False, False],
            }
        )
        ranked = enrichment.rank_table(frame, "Score", "Ligand", lower_is_better=True, best_per_ligand=True)
        metrics, summary, _, _, _ = enrichment.calculate_metrics(
            ranked,
            fractions=[],
            top_counts=[],
            score_column="Score",
            lower_is_better=True,
            ef_decoy_fraction=0.02,
            logauc_min_decoy_fraction=0.001,
        )

        self.assertIn("EF2_Percent", set(metrics["Metric"]))
        self.assertIn("ef2_percent", summary)

    def test_strip_again_suffix_collapses_duplicate_preparations(self):
        frame = pd.DataFrame(
            {
                "Score": [-7.0, -8.0, -6.0],
                "Ligand": ["CHEMBL1", "CHEMBL1-again2", "CHEMBL2-again1"],
                "Is_Active": [True, True, True],
            }
        )

        normalized = enrichment.strip_again_suffix(frame, "Ligand")
        ranked = enrichment.rank_table(
            normalized,
            "Score",
            "Ligand",
            lower_is_better=True,
            best_per_ligand=True,
        )

        self.assertEqual(ranked["Ligand"].tolist(), ["CHEMBL1", "CHEMBL2"])
        self.assertEqual(ranked["Score"].tolist(), [-8.0, -6.0])
        self.assertEqual(ranked.loc[0, "Raw_Ligand"], "CHEMBL1-again2")


if __name__ == "__main__":
    unittest.main()
