"""Tests for critical-sphere labelling.

DOCK6 discards any orientation that does not use a sphere from every labelled cluster, so the
labelling has to be exact: the wrong column, a silently empty cluster, or a corrupted sphere line
would cost a whole screening run with nothing naming the cause.
"""

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np

from docking_module.dock6_pipeline import dock6_pipeline as d6

SPHERE_HEADER = "DOCK spheres within 10.0 ang of ligands\ncluster     1   number of spheres in cluster    4\n"
# id, x, y, z, radius, surface atom, critical cluster, colour — sphgen's fixed-width layout.
SPHERES = [
    (274, 0.0, 0.0, 0.0, 1.400, 291),
    (286, 2.0, 0.0, 0.0, 1.401, 262),
    (290, 10.0, 0.0, 0.0, 2.265, 2432),
    (291, 20.0, 0.0, 0.0, 3.939, 1646),
]

MOL2 = """@<TRIPOS>MOLECULE
receptor
 4 0 2 0 0
PROTEIN
AMBER ff14SB

@<TRIPOS>ATOM
      1 N          0.5000    0.0000    0.0000 N.am      93 GLU    -0.4157
      2 O          0.0000    0.5000    0.0000 O.2       93 GLU    -0.5679
      3 CB         0.0000    0.0000   -8.0000 C.3       93 GLU     0.0145
      4 N         19.5000    0.0000    0.0000 N.am      95 LEU    -0.4157
@<TRIPOS>BOND
@<TRIPOS>SUBSTRUCTURE
"""


def write_spheres(path, spheres=SPHERES):
    with open(path, "w") as handle:
        handle.write(SPHERE_HEADER)
        for i, x, y, z, r, surf in spheres:
            handle.write(f"{i:5d}{x:10.5f}{y:10.5f}{z:10.5f}{r:8.3f}{surf:5d}{0:2d}{0:3d}\n")
    return path


def read_clusters(path):
    return [int(line.split()[6]) for line in open(path)
            if len(line.split()) == 8 and line.split()[0].isdigit()]


class LabelSphereTests(unittest.TestCase):
    def test_labels_only_spheres_within_the_radius(self):
        with tempfile.TemporaryDirectory() as tmp:
            sph = write_spheres(Path(tmp) / "selected_spheres.sph")
            counts, closest = d6.label_critical_spheres(sph, [np.array([[0.0, 0.0, 0.0]])], 3.0)
            self.assertEqual(counts, [2])            # the spheres at 0 and 2 Å
            self.assertAlmostEqual(closest[0], 0.0)
            self.assertEqual(read_clusters(sph), [1, 1, 0, 0])

    def test_preserves_every_other_field_and_the_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            sph = write_spheres(Path(tmp) / "s.sph")
            before = [line for line in open(sph)]
            d6.label_critical_spheres(sph, [np.array([[0.0, 0.0, 0.0]])], 3.0)
            after = [line for line in open(sph)]
            self.assertEqual(before[:2], after[:2])  # header untouched
            for old, new in zip(before[2:], after[2:]):
                o, n = old.split(), new.split()
                self.assertEqual(o[:6], n[:6])       # id, coordinates, radius, surface atom
                self.assertEqual(o[7], n[7])         # colour
                self.assertEqual(len(new.split()), 8)

    def test_several_clusters_number_from_one_and_do_not_overlap(self):
        with tempfile.TemporaryDirectory() as tmp:
            sph = write_spheres(Path(tmp) / "s.sph")
            counts, _ = d6.label_critical_spheres(
                sph, [np.array([[0.0, 0.0, 0.0]]), np.array([[20.0, 0.0, 0.0]])], 3.0)
            self.assertEqual(counts, [2, 1])
            self.assertEqual(read_clusters(sph), [1, 1, 0, 2])

    def test_a_sphere_claimed_by_two_clusters_goes_to_the_lower_one(self):
        with tempfile.TemporaryDirectory() as tmp:
            sph = write_spheres(Path(tmp) / "s.sph")
            counts, _ = d6.label_critical_spheres(
                sph, [np.array([[0.0, 0.0, 0.0]]), np.array([[0.0, 0.0, 0.0]])], 3.0)
            self.assertEqual(counts, [2, 0])         # DOCK6 allows one cluster per sphere

    def test_reports_the_nearest_sphere_when_nothing_is_in_range(self):
        with tempfile.TemporaryDirectory() as tmp:
            sph = write_spheres(Path(tmp) / "s.sph")
            counts, closest = d6.label_critical_spheres(sph, [np.array([[50.0, 0.0, 0.0]])], 3.0)
            self.assertEqual(counts, [0])
            self.assertAlmostEqual(closest[0], 30.0)


class SelectionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.mol2 = Path(self.tmp.name) / "receptor.mol2"
        self.mol2.write_text(MOL2)

    def tearDown(self):
        self.tmp.cleanup()

    def test_selection_with_an_atom_filter_picks_only_those_atoms(self):
        points, description = d6.resolve_critical_points(
            {"selection": "resid 93 and name N O"}, str(self.mol2))
        self.assertEqual(len(points), 2)             # not the CB 8 Å away
        self.assertIn("GLU93:N", description)
        self.assertNotIn("GLU93:CB", description)

    def test_selection_without_an_atom_filter_takes_the_whole_residue(self):
        points, _ = d6.resolve_critical_points({"selection": "resid 93"}, str(self.mol2))
        self.assertEqual(len(points), 3)

    def test_explicit_center_and_points_file(self):
        points, _ = d6.resolve_critical_points({"center": [1.0, 2.0, 3.0]}, str(self.mol2))
        np.testing.assert_allclose(points, [[1.0, 2.0, 3.0]])
        points, _ = d6.resolve_critical_points({"points": str(self.mol2)}, str(self.mol2))
        self.assertEqual(len(points), 4)

    def test_a_selection_matching_nothing_is_an_error(self):
        with self.assertRaises(ValueError) as caught:
            d6.resolve_critical_points({"selection": "resid 999"}, str(self.mol2))
        self.assertIn("matched no heavy atoms", str(caught.exception))

    def test_exactly_one_source_is_required(self):
        for spec in ({}, {"selection": "resid 93", "center": [0, 0, 0]}):
            with self.assertRaises(ValueError):
                d6.resolve_critical_points(spec, str(self.mol2))

    def test_one_mapping_or_a_list_of_mappings_is_accepted(self):
        parse = d6.Dock6Pipeline._parse_critical_spheres
        self.assertEqual(parse(None), [])
        self.assertEqual(len(parse({"selection": "resid 93"})), 1)
        self.assertEqual(len(parse([{"selection": "resid 93"}, {"center": [0, 0, 0]}])), 2)
        with self.assertRaises(ValueError):
            parse(["resid 93"])                      # a bare string is not a cluster
        with self.assertRaises(ValueError):
            parse([{"center": [0, 0, 0]}] * (d6.CRITICAL_MAX_CLUSTERS + 1))


class ExecutionEnvironmentTests(unittest.TestCase):
    """DOCK6 is reproducible within an execution environment but not across them.

    The same config run with an empty launcher (one rank) and with srun over many ranks returns
    different poses, and nothing in the parameters records which was used. These fields make that
    difference visible in the summary.
    """

    def test_no_launcher_is_one_rank(self):
        env = d6.describe_execution_environment("")
        self.assertEqual(env["mpi_ranks"], 1)
        self.assertEqual(env["dock6_launcher"], "")
        env = d6.describe_execution_environment(None)
        self.assertEqual(env["mpi_ranks"], 1)

    def test_rank_count_is_read_from_the_launcher_after_expansion(self):
        os.environ["SLURM_NTASKS"] = "32"
        try:
            env = d6.describe_execution_environment("srun --mpi=pmix -n ${SLURM_NTASKS}")
        finally:
            del os.environ["SLURM_NTASKS"]
        self.assertEqual(env["mpi_ranks"], 32)
        self.assertEqual(env["dock6_launcher_expanded"], "srun --mpi=pmix -n 32")
        self.assertEqual(env["dock6_launcher"], "srun --mpi=pmix -n ${SLURM_NTASKS}")

    def test_other_launcher_spellings(self):
        self.assertEqual(d6.describe_execution_environment("mpirun -n 8")["mpi_ranks"], 8)
        self.assertEqual(d6.describe_execution_environment("mpiexec --ntasks 4")["mpi_ranks"], 4)

    def test_an_unparsable_launcher_reports_no_rank_count_rather_than_guessing(self):
        env = d6.describe_execution_environment("some-site-wrapper")
        self.assertIsNone(env["mpi_ranks"])
        self.assertEqual(env["dock6_launcher"], "some-site-wrapper")

    def test_host_is_recorded(self):
        env = d6.describe_execution_environment("")
        self.assertTrue(env["hostname"])
        self.assertIn("cpu", env)


class FailFastTests(unittest.TestCase):
    """A bad critical-sphere spec must fail before any DOCK6 executable runs.

    Resolution happens in __init__, so an unusable selection costs nothing: without this, the run
    would build spheres with sphgen first and only then discover the problem.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        dock6 = root / "dock6" / "parameters"
        dock6.mkdir(parents=True)
        (root / "dock6" / "bin").mkdir()
        for name in ("vdw_AMBER_parm99.defn", "flex.defn", "flex_drive.tbl"):
            (dock6 / name).write_text("x\n")
        (root / "receptor.mol2").write_text(MOL2)
        for name in ("rec.dms", "ref.mol2", "ligands.mol2"):
            (root / name).write_text("x\n")
        self.root = root
        self.kwargs = dict(
            receptor_dms=str(root / "rec.dms"), receptor_mol2=str(root / "receptor.mol2"),
            ref_ligand_mol2=str(root / "ref.mol2"), ligand_file=str(root / "ligands.mol2"),
            out_dir=str(root / "out"), dock6_dir=str(root / "dock6"),
        )

    def tearDown(self):
        self.tmp.cleanup()

    def test_a_selection_matching_nothing_raises_at_construction(self):
        with self.assertRaises(ValueError) as caught:
            d6.Dock6Pipeline(critical_spheres={"selection": "resid 999"}, **self.kwargs)
        self.assertIn("matched no heavy atoms", str(caught.exception))
        # Nothing was handed to DOCK6: the sphere directory is still empty.
        self.assertEqual(list((self.root / "out" / "01_spheres").iterdir()), [])

    def test_a_missing_reference_raises_at_construction(self):
        with self.assertRaises(FileNotFoundError):
            d6.Dock6Pipeline(
                critical_spheres={"selection": "resid 93", "reference": str(self.root / "nope.pdb")},
                **self.kwargs)

    def test_a_valid_selection_is_resolved_once_up_front(self):
        pipeline = d6.Dock6Pipeline(critical_spheres={"selection": "resid 93 and name N O"},
                                    **self.kwargs)
        self.assertEqual(len(pipeline.critical_points), 1)
        coordinates, description = pipeline.critical_points[0]
        self.assertEqual(len(coordinates), 2)
        self.assertIn("GLU93:N", description)

    def test_without_critical_spheres_nothing_is_resolved(self):
        pipeline = d6.Dock6Pipeline(**self.kwargs)
        self.assertEqual(pipeline.critical_clusters, [])
        self.assertEqual(pipeline.critical_points, [])


if __name__ == "__main__":
    unittest.main()
