import argparse
import json
import logging
import os
import platform
import re
import shlex
import shutil
import subprocess
import sys
import warnings
from datetime import datetime
from pathlib import Path

import yaml

# DOCK6 reads file names and MOL2 title lines into fixed-width buffers. Past
# roughly 80 characters sphgen silently truncates the name it was given, and
# grid segfaults on an over-long receptor title. Neither failure names its own
# cause, so the pipeline keeps every generated name short rather than leaving
# the limit for a user to rediscover.
DOCK6_NAME_LIMIT = 80

# Intermediates live in a per-run working directory, so they do not need unique
# names. Fixed short ones keep the generated inputs under the limit no matter
# how long the receptor file is called. Only the final flex output and the log
# carry the descriptive prefix, because those are the documented artifacts.
SPHERE_FILE = "receptor.sph"
BOX_FILE = "receptor.box.pdb"
GRID_PREFIX = "receptor_grid"
SAFE_MOL2_TITLE = "receptor"

# Spheres are pseudo-atoms filling the pocket, so their centers sit about a hydrogen bond away
# from the atoms that define a site. 3 Å keeps the critical region to the spheres a ligand atom
# could actually occupy; much more and it reaches spheres on the far side of the named residues.
CRITICAL_RADIUS = 3.0
CRITICAL_MAX_CLUSTERS = 10      # DOCK's sphere format reserves one digit for the cluster number
CRITICAL_LOOSE_FRACTION = 0.5   # above this, the filter is not really restricting anything

MOLECULE_MARKER = "@<TRIPOS>MOLECULE"

DOCK6_BOOLEAN_PARAMETERS = {
    "automated_matching",
    "bump_filter",
    "calculate_rmsd",
    "chemical_matching",
    "contact_score_primary",
    "critical_points",
    "grid_lig_efficiency",
    "grid_score_primary",
    "limit_max_ligands",
    "minimize_anchor",
    "minimize_flexible_growth",
    "minimize_flexible_growth_ramp",
    "minimize_ligand",
    "orient_ligand",
    "rank_ligands",
    "read_mol_solvation",
    "score_molecules",
    "simplex_final_min",
    "simplex_restraint_min",
    "skip_molecule",
    "use_advanced_simplex_parameters",
    "use_clash_overlap",
    "use_database_filter",
    "use_internal_energy",
    "use_ligand_spheres",
    "use_rmsd_reference_mol",
    "write_growth_trees",
    "write_mol_solvation",
    "write_orientations",
}


def normalize_dock6_boolean_parameters(parameters):
    """Render YAML booleans as the yes/no values required by DOCK6."""
    normalized = dict(parameters)
    for key in DOCK6_BOOLEAN_PARAMETERS:
        if key not in normalized:
            continue

        value = normalized[key]
        if isinstance(value, bool):
            normalized[key] = "yes" if value else "no"
        elif isinstance(value, str) and value.lower() in {"yes", "no", "true", "false"}:
            normalized[key] = "yes" if value.lower() in {"yes", "true"} else "no"
        else:
            raise ValueError(
                f"DOCK6 parameter '{key}' must be a boolean or one of yes/no, not {value!r}."
            )
    return normalized


def dock6_parameter_is_yes(value):
    """Interpret a YAML/DOCK6 boolean value without changing the input mapping."""
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"yes", "true", "1"}


class MaxLogLevelFilter(logging.Filter):
    """Pass records at or below one logging level."""

    def __init__(self, maximum_level):
        super().__init__()
        self.maximum_level = maximum_level

    def filter(self, record):
        return record.levelno <= self.maximum_level


def configure_logging(log_file):
    formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s', datefmt='%H:%M:%S')
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
        handler.close()

    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.INFO)
    stdout_handler.addFilter(MaxLogLevelFilter(logging.WARNING))
    stdout_handler.setFormatter(formatter)
    root_logger.addHandler(stdout_handler)

    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(logging.ERROR)
    stderr_handler.setFormatter(formatter)
    root_logger.addHandler(stderr_handler)

    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)


def write_config_snapshot(log_file, config_data):
    with open(log_file, "a") as handle:
        handle.write("\n" + "=" * 60 + "\n")
        handle.write("--- CONFIG SNAPSHOT ---\n")
        handle.write("=" * 60 + "\n")
        handle.write(yaml.safe_dump(config_data, sort_keys=False))
        handle.write("\n")


def _json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _json_safe(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def write_json_summary(summary, json_out):
    json_path = Path(json_out)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    with json_path.open("w") as handle:
        json.dump(_json_safe(summary), handle, indent=2, sort_keys=False)
        handle.write("\n")


def load_pocket_center(pocket_definition):
    """Read a ANTIOCH-VS rank-1 pocket definition and return its center."""
    if not pocket_definition:
        return None
    path = Path(pocket_definition).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"pocket_definition does not exist: {path}")
    with path.open() as handle:
        definition = yaml.safe_load(handle) or {}
    center = definition.get("box", {}).get("center")
    if not isinstance(center, (list, tuple)) or len(center) != 3:
        raise ValueError(f"{path} must contain box.center with three numeric coordinates.")
    try:
        return path, tuple(float(value) for value in center)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{path} contains non-numeric box.center coordinates.") from exc


def read_mol2_title(path):
    """Return the title line of the first MOL2 record, or None if there is none."""
    with open(path, errors="replace") as handle:
        for line in handle:
            if line.startswith(MOLECULE_MARKER):
                return next(handle, "").rstrip("\n")
    return None


def write_mol2_with_safe_title(src, dest, title):
    """Copy a MOL2, replacing every record's title line with `title`.

    Only the line after each @<TRIPOS>MOLECULE is touched; atom, bond and
    substructure blocks are copied through unchanged. Never use this on a
    ligand file: there the title line is the ligand ID that
    extract_dock6_results.py reads back as `Ligand`.
    """
    with open(src, errors="replace") as source, open(dest, "w") as target:
        pending_title = False
        for line in source:
            if pending_title:
                target.write(f"{title}\n")
                pending_title = False
                continue
            target.write(line)
            if line.startswith(MOLECULE_MARKER):
                pending_title = True


def _load_structure(path):
    """Open a receptor/reference structure with MDAnalysis (MOL2 and PDB both work)."""
    import MDAnalysis as mda

    mda_logger = logging.getLogger("MDAnalysis")
    previous = mda_logger.level
    # MDAnalysis warns about missing elements and logs its mass/type guessing at INFO. Neither
    # matters here, since this only reads coordinates and residue identity, and both would
    # otherwise land in the run's own pipeline log.
    mda_logger.setLevel(logging.WARNING)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return mda.Universe(str(path))
    finally:
        mda_logger.setLevel(previous)


def map_residue_numbers(reference_path, receptor_path):
    """Map reference residue numbers onto a receptor numbered differently.

    Receptors in one benchmark rarely share a numbering scheme — a crystal may carry UniProt
    numbers while a predicted model starts at 1 — so a residue list written against one structure
    cannot be read literally against another. Sequences come from MDAnalysis and are aligned with
    Biopython; only positions whose residue type matches are mapped, so a questionable alignment
    drops a residue instead of silently pointing at a different one.
    """
    from Bio.Align import PairwiseAligner

    ref = _load_structure(reference_path).select_atoms("protein").residues
    rec = _load_structure(receptor_path).select_atoms("protein").residues
    ref_seq = ref.sequence(format="string")
    rec_seq = rec.sequence(format="string")
    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                              open_gap_score=-10, extend_gap_score=-0.5)
    alignment = aligner.align(ref_seq, rec_seq)[0]
    mapping = {}
    for (ref_start, ref_end), (rec_start, _) in zip(*alignment.aligned):
        for offset in range(ref_end - ref_start):
            if ref_seq[ref_start + offset] == rec_seq[rec_start + offset]:
                mapping[int(ref.resids[ref_start + offset])] = int(rec.resids[rec_start + offset])
    return mapping


def resolve_critical_points(spec, receptor_mol2):
    """Coordinates that define one critical cluster, plus a human-readable description.

    A cluster is given as exactly one of:
      selection      an MDAnalysis selection against the receptor, e.g. "resid 930 932 and name N O"
                     (the selection language carries the atom filter: `backbone`, `name OD1 OD2`,
                     `not backbone`, ...)
      selection + reference
                     the same, with residue numbers written in `reference`'s numbering and mapped
                     onto this receptor by sequence alignment
      points         a MOL2/PDB whose heavy atoms mark the site (a fragment, ion, or water)
      center         explicit [x, y, z]
    """
    import numpy as np

    given = [k for k in ("selection", "points", "center") if spec.get(k) is not None]
    if len(given) != 1:
        raise ValueError(
            "Each critical cluster needs exactly one of 'selection', 'points' or 'center'; "
            f"got {given or 'none'}."
        )

    if spec.get("center") is not None:
        center = spec["center"]
        if not isinstance(center, (list, tuple)) or len(center) != 3:
            raise ValueError("critical cluster 'center' must be three numeric coordinates.")
        return np.array([[float(v) for v in center]]), f"center {tuple(float(v) for v in center)}"

    if spec.get("points") is not None:
        path = os.path.abspath(os.path.expanduser(spec["points"]))
        if not os.path.isfile(path):
            raise FileNotFoundError(f"critical cluster 'points' file does not exist: {path}")
        atoms = _load_structure(path).select_atoms("not name H*")
        if atoms.n_atoms == 0:
            raise ValueError(f"No heavy atoms in critical cluster 'points' file: {path}")
        return atoms.positions, f"{atoms.n_atoms} atoms from {os.path.basename(path)}"

    selection = str(spec["selection"])
    reference = spec.get("reference")
    described = selection
    if reference:
        reference = os.path.abspath(os.path.expanduser(reference))
        if not os.path.isfile(reference):
            raise FileNotFoundError(f"critical cluster 'reference' does not exist: {reference}")
        mapping = map_residue_numbers(reference, receptor_mol2)
        resids = sorted({int(tok) for tok in re.findall(r"(?<=\bresid\s)([\d\s]+)", selection + " ")
                         for tok in tok.split()})
        if not resids:
            raise ValueError(
                "critical cluster gives a 'reference' but its selection has no 'resid' to map: "
                f"{selection!r}"
            )
        missing = [r for r in resids if r not in mapping]
        if missing:
            raise ValueError(
                f"residues {missing} from {os.path.basename(reference)} could not be mapped onto "
                f"{os.path.basename(receptor_mol2)} by sequence alignment (different construct?)."
            )
        for original in sorted(resids, reverse=True):  # longest first, so 93 is not hit by 9
            selection = re.sub(rf"(?<=\bresid\s)((?:\d+\s+)*){original}\b",
                               rf"\g<1>{mapping[original]}", selection)
        described = (f"{spec['selection']!r} in {os.path.basename(reference)} numbering "
                     f"-> {selection!r}")

    atoms = _load_structure(receptor_mol2).select_atoms(f"({selection}) and not name H*")
    if atoms.n_atoms == 0:
        raise ValueError(
            f"critical cluster selection matched no heavy atoms in "
            f"{os.path.basename(receptor_mol2)}: {selection!r}"
        )
    labels = sorted({f"{a.resname}{a.resid}:{a.name}" for a in atoms})
    return atoms.positions, f"{described} -> {atoms.n_atoms} atoms ({', '.join(labels[:6])}" \
                            f"{', ...' if len(labels) > 6 else ''})"


def describe_execution_environment(launcher):
    """What a rerun has to match to reproduce a run's poses.

    DOCK6 is reproducible for a given receptor, parameter set and execution environment, but not
    across environments: the same config run with an empty `dock6_launcher` (a single rank) and
    with `srun -n $SLURM_NTASKS` can return different poses, and so can a different CPU. None of
    that is visible in the parameters, so it is recorded here to make a surprising difference
    between two runs diagnosable rather than mysterious.
    """
    expanded = os.path.expandvars(launcher) if launcher else ""
    ranks = None
    match = re.search(r"(?:^|\s)-(?:n|-ntasks[= ])\s*(\d+)", expanded)
    if match:
        ranks = int(match.group(1))
    elif launcher:
        ranks = None            # a launcher we cannot parse: do not guess a number
    else:
        ranks = 1               # no launcher: dock6.mpi is invoked directly, as one rank
    return {
        "dock6_launcher": launcher or "",
        "dock6_launcher_expanded": expanded,
        "mpi_ranks": ranks,
        "slurm_ntasks": os.environ.get("SLURM_NTASKS"),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "hostname": platform.node(),
        "cpu": platform.processor() or platform.machine(),
    }


def read_sphere_count(path):
    """(header line count, number of sphere records) of a DOCK sphere file."""
    header = spheres = 0
    with open(path, errors="replace") as handle:
        lines = handle.readlines()
    for line in lines:
        fields = line.split()
        if len(fields) == 8 and fields[0].isdigit():
            spheres += 1
        else:
            header += 1
    return header, spheres


def label_critical_spheres(path, clusters, radius):
    """Write cluster numbers into column 7 of a DOCK sphere file.

    `clusters` is one array of points per critical cluster; a sphere whose CENTER lies within
    `radius` of any point of cluster i is labelled i+1. Spheres are pseudo-atoms filling the
    pocket, so their centers sit roughly a hydrogen bond away from the atoms that define the site:
    a sphere 2.9 Å from a backbone carbonyl is where a ligand's donor would sit.

    DOCK6 requires a valid orientation to use at least one sphere from EVERY labelled cluster, and
    a sphere belongs to at most one cluster, so where clusters overlap the lower-numbered one wins.
    Returns per-cluster sphere counts and the distance from each cluster to its nearest sphere.
    """
    import numpy as np

    header, spheres = [], []
    with open(path, errors="replace") as handle:
        lines = handle.readlines()
    for line in lines:
        fields = line.split()
        (spheres if len(fields) == 8 and fields[0].isdigit() else header).append((fields, line))
    if not spheres:
        raise RuntimeError(f"No spheres found in {path}; cannot label critical clusters.")

    centers = np.array([[float(f[1]), float(f[2]), float(f[3])] for f, _ in spheres])
    assigned = np.zeros(len(spheres), dtype=int)
    counts, closest = [], []
    for index, points in enumerate(clusters, start=1):
        distance = np.linalg.norm(centers[:, None] - np.asarray(points, float)[None], axis=2).min(1)
        hit = (distance <= radius) & (assigned == 0)
        assigned[hit] = index
        counts.append(int(hit.sum()))
        closest.append(float(distance.min()))

    with open(path, "w") as handle:
        handle.writelines(line for _, line in header)
        for (fields, _), cluster in zip(spheres, assigned):
            # Rebuilt in sphgen's fixed-width layout so DOCK6's reader stays happy.
            handle.write(f"{int(fields[0]):5d}{float(fields[1]):10.5f}{float(fields[2]):10.5f}"
                         f"{float(fields[3]):10.5f}{float(fields[4]):8.3f}{int(fields[5]):5d}"
                         f"{cluster:2d}{int(fields[7]):3d}\n")
    return counts, closest


def write_pseudo_reference_ligand(path, center):
    """Write a one-atom MOL2 reference ligand at a predicted pocket center."""
    x, y, z = center
    with open(path, "w") as handle:
        handle.write("@<TRIPOS>MOLECULE\nPOCKET_CENTER\n 1 0 0 0 0\nSMALL\nNO_CHARGES\n\n")
        handle.write("@<TRIPOS>ATOM\n")
        handle.write(f"      1 CEN {x:10.4f} {y:10.4f} {z:10.4f} C.3       1 PCK       0.0000\n")
        handle.write("@<TRIPOS>BOND\n@<TRIPOS>SUBSTRUCTURE\n")
        handle.write("     1 PCK         1 TEMP              0 ****  ****    0 ROOT\n")


class Dock6Pipeline:
    def __init__(self, receptor_dms, receptor_mol2, ref_ligand_mol2, ligand_file, out_dir, **kwargs):
        # 1. Resolve Absolute Paths
        self.receptor_dms = os.path.abspath(receptor_dms)
        self.receptor_mol2 = os.path.abspath(receptor_mol2)
        pocket_definition = kwargs.get("pocket_definition")
        self.pocket_definition = load_pocket_center(pocket_definition) if pocket_definition else None
        self.ref_ligand_mol2 = os.path.abspath(ref_ligand_mol2) if ref_ligand_mol2 else None
        self.ligand_file = os.path.abspath(ligand_file)
        self.out_dir = os.path.abspath(out_dir)
        rmsd_reference_mol2 = kwargs.get("rmsd_reference_mol2")
        self.rmsd_reference_mol2 = (
            os.path.abspath(rmsd_reference_mol2)
            if rmsd_reference_mol2
            else self.ref_ligand_mol2
        )

        # --- NEW: Verify All Input Files Exist ---
        required_files = {
            "receptor_dms": self.receptor_dms,
            "receptor_mol2": self.receptor_mol2,
            "ligand_file": self.ligand_file
        }
        if not self.pocket_definition:
            required_files["ref_ligand_mol2"] = self.ref_ligand_mol2
        if self.rmsd_reference_mol2:
            required_files["rmsd_reference_mol2"] = self.rmsd_reference_mol2

        missing_files = []
        for name, filepath in required_files.items():
            if not os.path.isfile(filepath):
                missing_files.append(f"{name} ({filepath})")

        if missing_files:
            error_msg = "\n[!] FATAL ERROR: The following required input files do not exist:\n"
            for mf in missing_files:
                error_msg += f"    - {mf}\n"
            raise FileNotFoundError(error_msg)
        # ------------------------------------------

        self.prefix = os.path.splitext(os.path.basename(self.receptor_mol2))[0]
        self.log_file = os.path.join(self.out_dir, f"{self.prefix}_dock6_pipeline.log")
        self.overwrite = kwargs.get("overwrite", False)
        self.run_config = kwargs.get("run_config", {})

        os.makedirs(self.out_dir, exist_ok=True)
        if self.overwrite and os.path.exists(self.log_file):
            os.remove(self.log_file)
        configure_logging(self.log_file)

        # 2. Setup Subdirectories & Handle Overwrite
        self.dirs = {
            "spheres": os.path.join(self.out_dir, "01_spheres"),
            "gridbox": os.path.join(self.out_dir, "02_gridbox"),
            "docking": os.path.join(self.out_dir, "03_flex_dock")
        }

        for d_name, d_path in self.dirs.items():
            if os.path.exists(d_path) and os.listdir(d_path):
                if self.overwrite:
                    logging.warning(f"Overwriting existing directory: {d_path}")
                    shutil.rmtree(d_path)
                else:
                    raise FileExistsError(
                        f"Directory {d_path} is not empty. DOCK6 will fail if output files already exist. "
                        f"Please use the --overwrite flag to clear previous results."
                    )
            os.makedirs(d_path, exist_ok=True)

        # 3. Extract Configuration Parameters
        self.config = {
            "sphere_radius": kwargs.get("sphere_radius", 15.0),
            "box_margin": kwargs.get("box_margin", 8.0),
            "dock6_dir": kwargs.get("dock6_dir", "/usr/local/dock6"),
            "dock6_executable": kwargs.get("dock6_executable"),
            "dock6_launcher": kwargs.get("dock6_launcher"),
            "max_orientations": kwargs.get("max_orientations", 1000),
            "sif_image": kwargs.get("sif_image"),
            "singularity_binds": kwargs.get("singularity_binds", []),
            "container_commands": kwargs.get("container_commands", {}),
            "dock6_name_limit": kwargs.get("dock6_name_limit", DOCK6_NAME_LIMIT),
            "critical_spheres": kwargs.get("critical_spheres"),
            "critical_radius": float(kwargs.get("critical_radius", CRITICAL_RADIUS)),
        }
        self.name_warnings = []
        self.critical_report = None
        self.critical_clusters = self._parse_critical_spheres(self.config["critical_spheres"])
        # Resolved now, not at sphere generation, so bad input fails before any DOCK6 run.
        self.critical_points = self._resolve_critical_spheres()

        self.use_singularity = bool(self.config["sif_image"])
        if self.use_singularity:
            self.config["sif_image"] = os.path.abspath(self.config["sif_image"])
            if not os.path.isfile(self.config["sif_image"]):
                raise FileNotFoundError(
                    f"Singularity image not found: {self.config['sif_image']}"
                )

            automatic_binds = [
                self.out_dir,
                os.path.dirname(self.receptor_dms),
                os.path.dirname(self.receptor_mol2),
                os.path.dirname(self.ligand_file),
            ]
            if self.ref_ligand_mol2:
                automatic_binds.append(os.path.dirname(self.ref_ligand_mol2))
            if self.rmsd_reference_mol2:
                automatic_binds.append(os.path.dirname(self.rmsd_reference_mol2))
            configured_binds = [os.path.abspath(path) for path in self.config["singularity_binds"]]
            self.singularity_binds = list(dict.fromkeys(configured_binds + automatic_binds))
            logging.info(
                "Using DOCK6 Singularity image: %s", self.config["sif_image"]
            )
        else:
            self.singularity_binds = []

        self.extra_flex_args = kwargs.get("extra_flex_args", {})
        if self.pocket_definition and not self.rmsd_reference_mol2:
            rmsd_requested = (
                dock6_parameter_is_yes(self.extra_flex_args.get("calculate_rmsd", "no"))
                and dock6_parameter_is_yes(self.extra_flex_args.get("use_rmsd_reference_mol", "no"))
            )
            if rmsd_requested:
                raise ValueError(
                    "A pocket_definition provides a pseudo-ligand for sphere selection, "
                    "not a valid pose-RMSD reference. Set calculate_rmsd and "
                    "use_rmsd_reference_mol to 'no', or provide paths.rmsd_reference_mol2."
                )

        # 4. Resolve Parameter Files
        params_dir = os.path.join(self.config["dock6_dir"], "parameters")
        self.vdw_defn_src = os.path.join(params_dir, "vdw_AMBER_parm99.defn")
        self.flex_defn_src = os.path.join(params_dir, "flex.defn")
        self.flex_drive_src = os.path.join(params_dir, "flex_drive.tbl")

        if not self.use_singularity:
            for p_file in [self.vdw_defn_src, self.flex_defn_src, self.flex_drive_src]:
                if not os.path.exists(p_file):
                    raise FileNotFoundError(f"DOCK6 parameter file missing: {p_file}")

        # 5. Define Executables
        bin_dir = os.path.join(self.config["dock6_dir"], "bin")
        dock6_executable = self.config["dock6_executable"]
        if dock6_executable:
            native_dock6 = (
                os.path.abspath(dock6_executable)
                if os.path.sep in dock6_executable
                else dock6_executable
            )
        else:
            native_dock6 = os.path.join(bin_dir, "dock6")

        native_executables = {
            "sphgen": os.path.join(bin_dir, "sphgen"),
            "sphere_selector": os.path.join(bin_dir, "sphere_selector"),
            "showbox": os.path.join(bin_dir, "showbox"),
            "grid": os.path.join(bin_dir, "grid"),
            "dock6": native_dock6
        }
        if self.use_singularity:
            self.exe = {
                name: self.config["container_commands"].get(name, name)
                for name in native_executables
            }
        else:
            self.exe = native_executables

    def _command(self, executable, arguments=""):
        """Build a native or Singularity-wrapped DOCK6 command."""
        command = [executable]
        if self.use_singularity:
            command = ["singularity", "exec"]
            if self.singularity_binds:
                command.extend(["-B", ",".join(self.singularity_binds)])
            command.extend([self.config["sif_image"], executable])

        rendered = shlex.join(command)
        return f"{rendered} {arguments}" if arguments else rendered

    def _dock6_command(self, arguments=""):
        """Build the final flexible-docking command, with optional MPI launcher."""
        command = self._command(self.exe["dock6"], arguments)
        launcher = self.config["dock6_launcher"]
        return f"{launcher} {command}" if launcher else command

    def _create_symlink(self, target, link_name, target_dir):
        """Creates a local symlink in the specified subdirectory."""
        link_path = os.path.join(target_dir, link_name)
        if os.path.exists(link_path) or os.path.islink(link_path):
            os.remove(link_path)
        os.symlink(target, link_path)
        return link_name

    @staticmethod
    def _parse_critical_spheres(spec):
        """Normalize the critical_spheres config into a list of cluster specs.

        Accepts a single mapping (one cluster, the common case) or a list of mappings.
        """
        if spec in (None, {}, []):
            return []
        clusters = spec if isinstance(spec, list) else [spec]
        for cluster in clusters:
            if not isinstance(cluster, dict):
                raise ValueError(
                    "critical_spheres must be a mapping, or a list of mappings for several "
                    f"clusters; got {type(cluster).__name__}."
                )
        if len(clusters) > CRITICAL_MAX_CLUSTERS:
            raise ValueError(
                f"DOCK6 allows at most {CRITICAL_MAX_CLUSTERS} critical clusters; got {len(clusters)}."
            )
        return clusters

    def _resolve_critical_spheres(self):
        """Resolve every critical cluster to coordinates, before any DOCK6 executable runs.

        Done during __init__ rather than at sphere generation so that a missing optional
        dependency, a selection that matches nothing, or an unmappable reference residue fails
        while the run is still free, instead of after sphgen has built a surface.
        """
        resolved = []
        for index, cluster in enumerate(self.critical_clusters, start=1):
            coordinates, description = resolve_critical_points(cluster, self.receptor_mol2)
            resolved.append((coordinates, description))
            logging.info("Critical cluster %d: %s", index, description)
        return resolved

    def _label_critical_spheres(self, spheres_path):
        """Mark the spheres of each critical cluster, and report what was marked."""
        radius = self.config["critical_radius"]
        points = [coordinates for coordinates, _ in self.critical_points]
        descriptions = [description for _, description in self.critical_points]

        counts, closest = label_critical_spheres(spheres_path, points, radius)
        total = sum(counts)

        empty = [i + 1 for i, n in enumerate(counts) if n == 0]
        if empty:
            # DOCK6 discards every orientation that misses a cluster, so an unlabelled cluster
            # would return zero poses for every ligand with nothing naming the cause.
            raise RuntimeError(
                "No spheres lie within {:.2f} Å of critical cluster(s) {} (nearest sphere: {}). "
                "DOCK6 would reject every orientation. Check the selection, or raise "
                "parameters.critical_radius.".format(
                    radius, ", ".join(map(str, empty)),
                    ", ".join(f"{closest[i - 1]:.2f} Å" for i in empty))
            )

        _, spheres = read_sphere_count(spheres_path)
        fraction = total / spheres if spheres else 0.0
        if fraction > CRITICAL_LOOSE_FRACTION:
            message = (f"{total} of {spheres} selected spheres ({fraction:.0%}) were marked "
                       f"critical; the filter restricts orientations very little at "
                       f"critical_radius {radius:.2f} Å.")
            logging.warning(message)
            self.name_warnings.append(message)

        self.critical_report = {
            "radius": radius,
            "n_clusters": len(counts),
            "spheres_per_cluster": counts,
            "nearest_sphere_per_cluster": [round(d, 3) for d in closest],
            "selected_spheres_total": spheres,
            "clusters": descriptions,
        }
        logging.info("Critical spheres: %s of %d selected spheres marked (radius %.2f Å)",
                     "+".join(map(str, counts)), spheres, radius)

    def _stage_receptor_mol2(self, link_name, target_dir):
        """Make the receptor MOL2 available locally with a title DOCK6 can read.

        A title within the limit is left alone and symlinked as before, so a
        normal run still points at the Chimera output. An over-long title is
        the documented `grid` segfault, so in that case a copy is staged with a
        short title instead; the original file is not modified.
        """
        title = read_mol2_title(self.receptor_mol2)
        limit = self.config["dock6_name_limit"]
        if title is None or len(title) <= limit:
            return self._create_symlink(self.receptor_mol2, link_name, target_dir)

        dest = os.path.join(target_dir, link_name)
        if os.path.islink(dest):
            os.remove(dest)
        write_mol2_with_safe_title(self.receptor_mol2, dest, SAFE_MOL2_TITLE)
        logging.warning(
            "Receptor MOL2 title is %d characters, over the DOCK6 limit of %d. "
            "Staged a copy titled %r for grid generation; %s is unchanged.",
            len(title), limit, SAFE_MOL2_TITLE, self.receptor_mol2,
        )
        self.name_warnings.append(
            f"Receptor MOL2 title ({len(title)} characters) exceeded the DOCK6 limit "
            f"of {limit}; a copy titled {SAFE_MOL2_TITLE!r} was staged for grid generation."
        )
        return link_name

    def _check_name_limits(self):
        """Record any generated DOCK6 input name that is still over the limit.

        Every name below is fixed, so this cannot fail today. It is here so
        that a future change which reintroduces the receptor prefix into a
        generated name shows up as a warning instead of a truncated sphgen
        path or a segfault.
        """
        limit = self.config["dock6_name_limit"]
        generated = {
            "sphere_file": SPHERE_FILE,
            "box_file": BOX_FILE,
            "grid_prefix": GRID_PREFIX,
            "flex_grid_prefix": f"../02_gridbox/{GRID_PREFIX}",
            "flex_spheres": "../01_spheres/selected_spheres.sph",
        }
        for name, value in generated.items():
            if len(value) > limit:
                message = (
                    f"Generated DOCK6 name {name}={value!r} is {len(value)} characters, "
                    f"over the limit of {limit}. sphgen will truncate it."
                )
                logging.warning(message)
                self.name_warnings.append(message)

    def _run_command(self, cmd, step_name, work_dir, use_shell=False):
        """Executes a command safely inside a specific subfolder."""
        logging.info(f"STARTING: {step_name} (in {os.path.basename(work_dir)})")
        cwd = os.getcwd()
        os.chdir(work_dir)

        try:
            result = subprocess.run(cmd, shell=use_shell, capture_output=True, text=True, executable='/bin/bash' if use_shell else None)

            with open(self.log_file, 'a') as f:
                f.write(f"\n{'='*50}\n--- {step_name.upper()} ---\n{'='*50}\n")
                f.write(result.stdout)
                if result.stderr:
                    f.write("\n[STDERR]:\n" + result.stderr)

            if result.returncode != 0:
                logging.error(f"FAILED: {step_name}. Check {self.log_file} for details.")
                raise subprocess.CalledProcessError(result.returncode, cmd)

        finally:
            os.chdir(cwd)

        logging.info(f"FINISHED: {step_name}")

    def step_1_spheres(self):
        logging.info("\n--- STEP 1: GENERATING SPHERES ---")
        wd = self.dirs["spheres"]

        local_dms = self._create_symlink(self.receptor_dms, "receptor.dms", wd)
        if self.pocket_definition:
            _, center = self.pocket_definition
            ref_path = os.path.join(wd, "pocket_center.mol2")
            write_pseudo_reference_ligand(ref_path, center)
            local_ref = os.path.basename(ref_path)
            logging.info("Using pocket-definition pseudo-reference ligand at %.3f %.3f %.3f", *center)
        else:
            local_ref = self._create_symlink(self.ref_ligand_mol2, "ref_ligand.mol2", wd)
        sph_name = SPHERE_FILE

        self._check_name_limits()

        with open(os.path.join(wd, "INSPH"), 'w') as f:
            f.write(f"{local_dms}\nR\nX\n0.0\n4.0\n1.4\n{sph_name}\n")

        self._run_command(self._command(self.exe['sphgen'], "-i INSPH -o OUTSPH"), "Sphgen", wd, use_shell=True)
        self._run_command(
            self._command(
                self.exe['sphere_selector'],
                f"{shlex.quote(sph_name)} {shlex.quote(local_ref)} {self.config['sphere_radius']}"
            ),
            "Sphere Selector",
            wd,
            use_shell=True
        )

        selected = os.path.join(wd, "selected_spheres.sph")
        if not os.path.exists(selected):
            raise RuntimeError("Sphere selection failed: selected_spheres.sph not found.")

        if self.critical_clusters:
            self._label_critical_spheres(selected)

    def step_2_box(self):
        logging.info("\n--- STEP 2: GENERATING BOX ---")
        wd = self.dirs["gridbox"]
        rel_spheres = "../01_spheres/selected_spheres.sph"
        box_name = BOX_FILE

        with open(os.path.join(wd, "showbox.in"), 'w') as f:
            f.write(f"Y\n{self.config['box_margin']}\n{rel_spheres}\n1\n{box_name}\n")

        self._run_command(self._command(self.exe['showbox'], "< showbox.in"), "ShowBox", wd, use_shell=True)

    def step_3_grid(self):
        logging.info("\n--- STEP 3: GENERATING GRID ---")
        wd = self.dirs["gridbox"]

        local_rec = self._stage_receptor_mol2("receptor.mol2", wd)
        local_vdw = (
            self.vdw_defn_src
            if self.use_singularity
            else self._create_symlink(self.vdw_defn_src, "vdw.defn", wd)
        )

        box_name = BOX_FILE
        grid_prefix = GRID_PREFIX

        content = f"""allow_non_integral_charges yes
compute_grids yes
grid_spacing 0.3
output_molecule no
contact_score no
energy_score yes
energy_cutoff_distance 10
atom_model a
attractive_exponent 6
repulsive_exponent 9
distance_dielectric yes
dielectric_factor 4
bump_filter yes
bump_overlap 0.75
receptor_file {local_rec}
box_file {box_name}
vdw_definition_file {local_vdw}
score_grid_prefix {grid_prefix}
"""
        with open(os.path.join(wd, "grid.in"), 'w') as f:
            f.write(content)

        self._run_command(self._command(self.exe['grid'], "-i grid.in -o grid.out"), "Grid Generation", wd, use_shell=True)

    def step_4_flex_dock(self):
        logging.info("\n--- STEP 4: FLEXIBLE DOCKING ---")
        wd = self.dirs["docking"]

        local_lig = self._create_symlink(self.ligand_file, "ligands.mol2", wd)
        if self.use_singularity:
            local_vdw = self.vdw_defn_src
            local_flex = self.flex_defn_src
            local_drive = self.flex_drive_src
        else:
            local_vdw = self._create_symlink(self.vdw_defn_src, "vdw.defn", wd)
            local_flex = self._create_symlink(self.flex_defn_src, "flex.defn", wd)
            local_drive = self._create_symlink(self.flex_drive_src, "flex_drive.tbl", wd)

        rel_spheres = "../01_spheres/selected_spheres.sph"
        rel_grid_prefix = f"../02_gridbox/{GRID_PREFIX}"
        out_prefix = f"{self.prefix}_flex"

        # A. Set the baseline parameters from the tutorial
        flex_params = {
            "conformer_search_type": "flex",
            "pruning_orient_score_cutoff": 1000.0,
            "pruning_conformer_score_cutoff": 100.0,
            "pruning_conformer_score_scaling_factor": 1.0,
            "use_clash_overlap": "no",
            "write_growth_trees": "no",
            "use_internal_energy": "yes",
            "internal_energy_rep_exp": 12,
            "internal_energy_cutoff": 100.0,
            "calculate_rmsd": "no",
            "use_rmsd_reference_mol": "no",
            "use_database_filter": "no",
            "orient_ligand": "yes",
            "automated_matching": "yes",
            "max_orientations": self.config['max_orientations'],
            # Turned on automatically when critical_spheres is configured; an explicit
            # extra_flex_args entry still wins, so a run can label spheres without filtering on them.
            "critical_points": "yes" if self.critical_clusters else "no",
            "chemical_matching": "no",
            "use_ligand_spheres": "no",
            "bump_filter": "no",
            "score_molecules": "yes",
            "contact_score_primary": "no",
            "grid_score_primary": "yes",
            "grid_score_rep_rad_scale": 1,
            "grid_score_vdw_scale": 1,
            "grid_score_es_scale": 1,
            "grid_lig_efficiency": "no",
            "minimize_ligand": "yes",
            "minimize_anchor": "yes",
            "minimize_flexible_growth": "yes",
            "use_advanced_simplex_parameters": "no",
            "minimize_flexible_growth_ramp": "yes",
            "simplex_max_cycles": 1,
            "simplex_score_converge": 0.1,
            "simplex_cycle_converge": 1.0,
            "simplex_trans_step": 1.0,
            "simplex_rot_step": 0.1,
            "simplex_tors_step": 10.0,
            "simplex_anchor_max_iterations": 500,
            "simplex_grow_max_iterations": 500,
            "simplex_grow_tors_premin_iterations": 0,
            "simplex_final_min": "no",
            "simplex_random_seed": 0,
            "simplex_restraint_min": "no",
            "atom_model": "all",
            "write_mol_solvation": "no",
            "write_orientations": "no",
            "num_final_scored_poses": 1,
            "num_preclustered_conformers": 1,
            "score_threshold": 100.0,
            "rank_ligands": "no"
        }

        # B. Override defaults with any parameters the user provided in the YAML
        flex_params.update(self.extra_flex_args)
        flex_params = normalize_dock6_boolean_parameters(flex_params)

        # C. Forcibly set the dynamic file paths (Protects pipeline integrity)
        flex_params.pop("rmsd_reference_filename", None)
        flex_params.update({
            "ligand_atom_file": local_lig,
            "receptor_site_file": rel_spheres,
            "grid_score_grid_prefix": rel_grid_prefix,
            "vdw_defn_file": local_vdw,
            "flex_defn_file": local_flex,
            "flex_drive_file": local_drive,
            "ligand_outfile_prefix": out_prefix
        })

        if flex_params["calculate_rmsd"] == "yes" and flex_params["use_rmsd_reference_mol"] == "yes":
            flex_params["rmsd_reference_filename"] = self._create_symlink(
                self.rmsd_reference_mol2, "rmsd_reference.mol2", wd
            )

        # Construct the flex.in file format natively
        content = ""
        for key, val in flex_params.items():
            # Formats it beautifully: Key left-aligned to 45 chars, value right after
            content += f"{key:<45} {val}\n"

        with open(os.path.join(wd, "flex.in"), 'w') as f:
            f.write(content)

        self._run_command(self._dock6_command("-i flex.in -o flex.out"), "Flexible Docking", wd, use_shell=True)

    def build_summary(self, status="success", error=None):
        docking_dir = self.dirs["docking"]
        flex_outputs = sorted(str(path) for path in Path(docking_dir).glob(f"{self.prefix}_flex*"))
        selected_spheres = Path(self.dirs["spheres"]) / "selected_spheres.sph"
        box_file = Path(self.dirs["gridbox"]) / BOX_FILE
        grid_files = sorted(str(path) for path in Path(self.dirs["gridbox"]).glob(f"{GRID_PREFIX}*"))

        warnings = list(self.name_warnings)
        if status == "success" and not flex_outputs:
            warnings.append("No final DOCK6 flex output files were found.")

        return {
            "stage": "dock6",
            "status": status,
            "timestamp": datetime.now().isoformat(),
            "inputs": {
                "receptor_dms": self.receptor_dms,
                "receptor_mol2": self.receptor_mol2,
                "reference_ligand_mol2": self.ref_ligand_mol2,
                "pocket_definition": str(self.pocket_definition[0]) if self.pocket_definition else None,
                "ligand_file": self.ligand_file,
            },
            "configuration": {
                "sphere_radius": self.config.get("sphere_radius"),
                "box_margin": self.config.get("box_margin"),
                "max_orientations": self.config.get("max_orientations"),
                "use_singularity": self.use_singularity,
                "dock6_name_limit": self.config.get("dock6_name_limit"),
                "critical_spheres": self.critical_report,
                "execution": describe_execution_environment(self.config.get("dock6_launcher")),
                "extra_flex_args": self.extra_flex_args,
                "pocket_center": self.pocket_definition[1] if self.pocket_definition else None,
            },
            "outputs": {
                "output_dir": self.out_dir,
                "log_file": self.log_file,
                "spheres_dir": self.dirs["spheres"],
                "gridbox_dir": self.dirs["gridbox"],
                "docking_dir": docking_dir,
                "selected_spheres": str(selected_spheres),
                "box_file": str(box_file),
                "grid_files": grid_files,
                "flex_output_files": flex_outputs,
            },
            "metrics": {
                "grid_file_count": len(grid_files),
                "flex_output_file_count": len(flex_outputs),
            },
            "warnings": warnings,
            "error": error,
        }

    def run_all(self):
        logging.info("!!! DOCK6 PIPELINE INITIALIZED !!!")

        with open(self.log_file, 'a') as f:
            f.write("Starting DOCK6 Pipeline\n")
        write_config_snapshot(self.log_file, self.run_config)

        self.step_1_spheres()
        self.step_2_box()
        self.step_3_grid()
        self.step_4_flex_dock()

        logging.info("!!! ALL STEPS COMPLETED SUCCESSFULLY !!!")
        logging.info(f"Final docked models saved as: {self.prefix}_flex_scored.mol2 in {self.dirs['docking']}")
        return self.build_summary()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DOCK6 Batch Preparation and Docking Pipeline")
    parser.add_argument("--config", required=True, help="Path to the YAML configuration file")
    parser.add_argument("--overwrite", action="store_true", help="Force overwrite of existing output directories")
    parser.add_argument("--json-out", help="Optional path for a machine-readable JSON run summary")
    args = parser.parse_args()

    try:
        with open(args.config) as file:
            config_data = yaml.safe_load(file)
    except Exception as e:
        raise ValueError(f"Error reading the YAML file:\n{e}")

    paths = config_data.get('paths', {})
    params = config_data.get('parameters', {})
    env = config_data.get('environment', {})

    pipeline = None
    json_out = args.json_out or config_data.get("reporting", {}).get("json_out")
    try:
        pipeline = Dock6Pipeline(
            receptor_dms=paths.get('receptor_dms'),
            receptor_mol2=paths.get('receptor_mol2'),
            ref_ligand_mol2=paths.get('ref_ligand_mol2'),
            rmsd_reference_mol2=paths.get('rmsd_reference_mol2'),
            ligand_file=paths.get('ligand_file'),
            out_dir=paths.get('output_dir'),

            pocket_definition=paths.get('pocket_definition'),

            overwrite=args.overwrite,
            sphere_radius=params.get('sphere_radius', 15.0),
            box_margin=params.get('box_margin', 8.0),
            max_orientations=params.get('max_orientations', 1000),
            dock6_name_limit=params.get('dock6_name_limit', DOCK6_NAME_LIMIT),
            critical_spheres=params.get('critical_spheres'),
            critical_radius=params.get('critical_radius', CRITICAL_RADIUS),

            extra_flex_args=params.get('extra_flex_args', {}),

            dock6_dir=env.get('dock6_dir', '/usr/local/dock6'),
            dock6_executable=env.get('dock6_executable'),
            dock6_launcher=env.get('dock6_launcher'),
            sif_image=env.get('sif_image'),
            singularity_binds=env.get('singularity_binds', []),
            container_commands=env.get('container_commands', {}),
            run_config=config_data
        )
        summary = pipeline.run_all()
        write_json_summary(summary, json_out or Path(pipeline.out_dir) / "dockflow_dock6_summary.json")
    except Exception as exc:
        if pipeline:
            failure = pipeline.build_summary(status="failed", error=str(exc))
            write_json_summary(failure, json_out or Path(pipeline.out_dir) / "dockflow_dock6_summary.json")
        raise
