"""Tests for the DOCK6 name-length guards.

DOCK6 reads file names and MOL2 title lines into fixed-width buffers. Over the
limit, sphgen truncates the name it was given and grid segfaults on the receptor
title — neither failure names its own cause, so these tests pin the behaviour
that keeps the generated inputs short.
"""

import tempfile
import unittest
from pathlib import Path

from docking_module.dock6_pipeline import dock6_pipeline as d6

RECEPTOR_MOL2 = """@<TRIPOS>MOLECULE
{title}
 3 2 1 0 0
PROTEIN
AMBER ff14SB

@<TRIPOS>ATOM
      1 N          1.0000     2.0000     3.0000 N.4       1 ALA1       0.1000
      2 CA         2.0000     3.0000     4.0000 C.3       1 ALA1       0.2000
      3 C          3.0000     4.0000     5.0000 C.2       1 ALA1       0.3000
@<TRIPOS>BOND
     1    1    2 1
     2    2    3 1
@<TRIPOS>SUBSTRUCTURE
     1 ALA1        1 RESIDUE           4 A     ALA     1 ROOT
"""


def write_receptor(path: Path, title: str) -> Path:
    path.write_text(RECEPTOR_MOL2.format(title=title))
    return path


class Mol2TitleTests(unittest.TestCase):
    def test_reads_the_title_line_of_the_first_record(self):
        with tempfile.TemporaryDirectory() as temporary:
            receptor = write_receptor(Path(temporary) / "r.mol2", "some_model.pdb")
            self.assertEqual(d6.read_mol2_title(receptor), "some_model.pdb")

    def test_returns_none_when_there_is_no_molecule_record(self):
        with tempfile.TemporaryDirectory() as temporary:
            empty = Path(temporary) / "empty.mol2"
            empty.write_text("@<TRIPOS>ATOM\n")
            self.assertIsNone(d6.read_mol2_title(empty))

    def test_rewrites_every_title_and_leaves_the_rest_byte_identical(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            long_title = "x" * 120
            source = write_receptor(root / "r.mol2", long_title)
            # Two records, so a multi-record file cannot silently keep one title.
            source.write_text(source.read_text() * 2)
            destination = root / "staged.mol2"

            d6.write_mol2_with_safe_title(source, destination, "receptor")

            written = destination.read_text().splitlines()
            titles = [
                written[index + 1]
                for index, line in enumerate(written)
                if line.startswith(d6.MOLECULE_MARKER)
            ]
            self.assertEqual(titles, ["receptor", "receptor"])
            self.assertNotIn(long_title, destination.read_text())

            # Everything that is not a title line survives the copy unchanged.
            original = source.read_text().splitlines()
            kept_source = [
                line
                for index, line in enumerate(original)
                if not (index and original[index - 1].startswith(d6.MOLECULE_MARKER))
            ]
            kept_written = [
                line
                for index, line in enumerate(written)
                if not (index and written[index - 1].startswith(d6.MOLECULE_MARKER))
            ]
            self.assertEqual(kept_source, kept_written)


class GeneratedNameTests(unittest.TestCase):
    def test_generated_names_stay_under_the_limit(self):
        """The names DOCK6 is handed are fixed, so no receptor name can lengthen them."""
        generated = [
            d6.SPHERE_FILE,
            d6.BOX_FILE,
            d6.GRID_PREFIX,
            f"../02_gridbox/{d6.GRID_PREFIX}",
            "../01_spheres/selected_spheres.sph",
        ]
        for name in generated:
            with self.subTest(name=name):
                self.assertLessEqual(len(name), d6.DOCK6_NAME_LIMIT)

    def test_a_pathological_receptor_name_does_not_reach_the_generated_names(self):
        pathological = "a" * 200
        for name in (d6.SPHERE_FILE, d6.BOX_FILE, d6.GRID_PREFIX):
            with self.subTest(name=name):
                self.assertNotIn(pathological, name)


if __name__ == "__main__":
    unittest.main()
