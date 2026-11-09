"""Tests for Vina CPU/GPU ligand PDBQT compatibility classification."""

import csv
import json
import tempfile
import unittest
from pathlib import Path

from docking_module.vina_pipeline import classify_vina_ligands as classifier


def pdbqt_atom(atom_type: str) -> str:
    """Return a minimal, fixed-width PDBQT ATOM record."""
    return f"ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00    0.000 {atom_type:>2}\n"


class ClassifyVinaLigandsTests(unittest.TestCase):
    def test_classifies_gpu_compatible_and_cpu_only_atom_types(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            compatible = root / "compatible.pdbqt"
            cpu_only = root / "cpu_only.pdbqt"
            compatible.write_text(
                pdbqt_atom("C") + pdbqt_atom("NA") + pdbqt_atom("OA") + pdbqt_atom("Se")
            )
            cpu_only.write_text(pdbqt_atom("C") + pdbqt_atom("CG0") + pdbqt_atom("G0") + pdbqt_atom("Si"))

            gpu_result = classifier.classify_pdbqt(
                compatible, classifier.VINA_GPU_2_1_LEGACY_TYPES, classifier.VINA_CPU_CURRENT_TYPES
            )
            cpu_result = classifier.classify_pdbqt(
                cpu_only, classifier.VINA_GPU_2_1_LEGACY_TYPES, classifier.VINA_CPU_CURRENT_TYPES
            )

            self.assertEqual(gpu_result.status, "gpu_compatible")
            self.assertEqual(cpu_result.status, "cpu_only")
            self.assertEqual(cpu_result.unsupported_gpu_types, ("CG0", "G0", "Si"))

    def test_boron_and_nonstandard_types_are_not_offered_to_vina(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            unsupported = root / "unsupported.pdbqt"
            unsupported.write_text(pdbqt_atom("C") + pdbqt_atom("B") + pdbqt_atom("NS"))

            result = classifier.classify_pdbqt(
                unsupported, classifier.VINA_GPU_2_1_LEGACY_TYPES, classifier.VINA_CPU_CURRENT_TYPES
            )

            self.assertEqual(result.status, "vina_unsupported")
            self.assertEqual(result.unsupported_cpu_types, ("B", "NS"))

    def test_default_writes_overlapping_cpu_and_gpu_input_directories(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_dir = root / "prepared"
            input_dir.mkdir()
            (input_dir / "standard.pdbqt").write_text(pdbqt_atom("C"))
            (input_dir / "macrocycle.pdbqt").write_text(pdbqt_atom("CG1"))
            (input_dir / "boron.pdbqt").write_text(pdbqt_atom("B"))
            output_dir = root / "classified"

            result = classifier.main(["--input", str(input_dir), "--output-dir", str(output_dir)])

            self.assertEqual(result, 0)
            self.assertTrue((output_dir / "vina_gpu_compatible" / "standard.pdbqt").is_file())
            self.assertFalse((output_dir / "vina_gpu_compatible" / "macrocycle.pdbqt").exists())
            self.assertTrue((output_dir / "vina_cpu_compatible" / "standard.pdbqt").is_file())
            self.assertTrue((output_dir / "vina_cpu_compatible" / "macrocycle.pdbqt").is_file())
            self.assertFalse((output_dir / "vina_cpu_compatible" / "boron.pdbqt").exists())
            self.assertTrue((output_dir / "vina_unsupported" / "boron.pdbqt").is_file())

            with (output_dir / "vina_ligand_compatibility.csv").open(newline="") as handle:
                rows = {row["filename"]: row for row in csv.DictReader(handle)}
            self.assertEqual(rows["standard.pdbqt"]["status"], "gpu_compatible")
            self.assertEqual(rows["macrocycle.pdbqt"]["unsupported_gpu_types"], "CG1")
            self.assertEqual(rows["boron.pdbqt"]["unsupported_cpu_types"], "B")
            self.assertEqual(
                rows["standard.pdbqt"]["engine_destinations"],
                "vina_gpu_compatible;vina_cpu_compatible",
            )

            summary = json.loads((output_dir / "vina_ligand_compatibility_summary.json").read_text())
            self.assertEqual(summary["cpu_ligand_set"], "all_compatible")
            self.assertEqual(summary["counts"]["gpu_cpu_overlap"], 1)

    def test_cpu_only_set_writes_disjoint_engine_input_directories(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_dir = root / "prepared"
            input_dir.mkdir()
            (input_dir / "standard.pdbqt").write_text(pdbqt_atom("C"))
            (input_dir / "macrocycle.pdbqt").write_text(pdbqt_atom("CG1"))
            output_dir = root / "classified"

            result = classifier.main(
                [
                    "--input",
                    str(input_dir),
                    "--output-dir",
                    str(output_dir),
                    "--cpu-ligand-set",
                    "cpu_only",
                ]
            )

            self.assertEqual(result, 0)
            self.assertTrue((output_dir / "vina_gpu_compatible" / "standard.pdbqt").is_file())
            self.assertFalse((output_dir / "vina_cpu_compatible" / "standard.pdbqt").exists())
            self.assertTrue((output_dir / "vina_cpu_compatible" / "macrocycle.pdbqt").is_file())

            with (output_dir / "vina_ligand_compatibility.csv").open(newline="") as handle:
                rows = {row["filename"]: row for row in csv.DictReader(handle)}
            self.assertEqual(rows["standard.pdbqt"]["engine_destinations"], "vina_gpu_compatible")
            self.assertEqual(rows["macrocycle.pdbqt"]["engine_destinations"], "vina_cpu_compatible")

            summary = json.loads((output_dir / "vina_ligand_compatibility_summary.json").read_text())
            self.assertEqual(summary["cpu_ligand_set"], "cpu_only")
            self.assertEqual(summary["counts"]["vina_cpu_input"], 1)
            self.assertEqual(summary["counts"]["gpu_cpu_overlap"], 0)

    def test_invalid_file_is_not_offered_to_either_engine(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_dir = root / "prepared"
            input_dir.mkdir()
            (input_dir / "broken.pdbqt").write_text("REMARK no atoms\n")
            output_dir = root / "classified"

            classifier.main(["--input", str(input_dir), "--output-dir", str(output_dir)])

            self.assertTrue((output_dir / "invalid_pdbqt" / "broken.pdbqt").is_file())
            self.assertFalse((output_dir / "vina_cpu_compatible" / "broken.pdbqt").exists())
            self.assertFalse((output_dir / "vina_gpu_compatible" / "broken.pdbqt").exists())
