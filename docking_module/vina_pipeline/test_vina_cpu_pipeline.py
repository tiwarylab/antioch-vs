"""Logging behavior for the Vina CPU runner."""

import contextlib
import io
import logging
import tempfile
import unittest
from pathlib import Path

from docking_module.dock6_pipeline import dock6_pipeline
from docking_module.vina_pipeline import vina_cpu_pipeline
from docking_module.vina_pipeline import vina_gpu_pipeline


class DockingLoggingTests(unittest.TestCase):
    def test_info_uses_stdout_and_errors_use_stderr(self):
        for pipeline in (vina_cpu_pipeline, vina_gpu_pipeline, dock6_pipeline):
            with self.subTest(pipeline=pipeline.__name__):
                root_logger = logging.getLogger()
                original_handlers = root_logger.handlers[:]
                original_level = root_logger.level
                stdout = io.StringIO()
                stderr = io.StringIO()

                try:
                    with tempfile.TemporaryDirectory() as temporary:
                        log_file = Path(temporary) / "pipeline.log"
                        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                            pipeline.configure_logging(log_file)
                            logging.info("normal progress")
                            logging.warning("normal warning")
                            logging.error("real failure")

                        log_text = log_file.read_text()
                finally:
                    for handler in root_logger.handlers[:]:
                        root_logger.removeHandler(handler)
                        handler.close()
                    for handler in original_handlers:
                        root_logger.addHandler(handler)
                    root_logger.setLevel(original_level)

                self.assertIn("normal progress", stdout.getvalue())
                self.assertIn("normal warning", stdout.getvalue())
                self.assertNotIn("real failure", stdout.getvalue())
                self.assertIn("real failure", stderr.getvalue())
                self.assertIn("normal progress", log_text)
                self.assertIn("normal warning", log_text)
                self.assertIn("real failure", log_text)

    def test_dock6_overwrite_warning_uses_stdout(self):
        root_logger = logging.getLogger()
        original_handlers = root_logger.handlers[:]
        original_level = root_logger.level
        stdout = io.StringIO()
        stderr = io.StringIO()

        try:
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                receptor_dms = root / "receptor.dms"
                receptor_mol2 = root / "receptor.mol2"
                reference_mol2 = root / "reference.mol2"
                ligand_mol2 = root / "ligands.mol2"
                for path in (receptor_dms, receptor_mol2, reference_mol2, ligand_mol2):
                    path.write_text("placeholder\n")

                out_dir = root / "dock6_output"
                stale_dir = out_dir / "01_spheres"
                stale_dir.mkdir(parents=True)
                (stale_dir / "stale_file").write_text("old\n")
                sif_image = root / "dock6.sif"
                sif_image.write_text("placeholder\n")

                with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                    dock6_pipeline.Dock6Pipeline(
                        receptor_dms=str(receptor_dms),
                        receptor_mol2=str(receptor_mol2),
                        ref_ligand_mol2=str(reference_mol2),
                        ligand_file=str(ligand_mol2),
                        out_dir=str(out_dir),
                        overwrite=True,
                        sif_image=str(sif_image),
                    )
        finally:
            for handler in root_logger.handlers[:]:
                root_logger.removeHandler(handler)
                handler.close()
            for handler in original_handlers:
                root_logger.addHandler(handler)
            root_logger.setLevel(original_level)

        self.assertIn("Overwriting existing directory", stdout.getvalue())
        self.assertEqual(stderr.getvalue(), "")
