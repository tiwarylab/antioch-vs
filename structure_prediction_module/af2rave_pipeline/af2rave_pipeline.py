#!/usr/bin/env python3
"""
AF2RAVE structure prediction and feature-selection pipeline.

This script is an orchestration layer around external LocalColabFold and
AF2RAVE installations. It does not install or vendor those tools.
"""

import argparse
import csv
import hashlib
import inspect
import json
import logging
import math
import platform
import re
import shutil
from datetime import datetime, timezone
from importlib import import_module, metadata
from pathlib import Path

import yaml

PLDDT_BFACTOR_SCALES = {"0_100": (0.0, 100.0, 1.0), "0_1": (0.0, 1.0, 100.0)}


class CollectedStructures(list):
    """Collected PDB paths plus their source-relative ordering fingerprint input."""

    def __init__(self, paths: list[Path], ordering_entries: list[str]):
        super().__init__(paths)
        self.ordering_entries = ordering_entries


def load_config(config_path: str) -> dict:
    with open(config_path) as handle:
        return yaml.safe_load(handle)


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def configure_logging(log_file: Path) -> None:
    formatter = logging.Formatter("[%(asctime)s] %(message)s", datefmt="%H:%M:%S")
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
        handler.close()

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)


def write_config_snapshot(config: dict) -> None:
    """Write the configuration through the configured logging handlers."""
    logging.info("%s", "=" * 60)
    logging.info("--- CONFIG SNAPSHOT ---")
    logging.info("%s", "=" * 60)
    for line in yaml.safe_dump(config, sort_keys=False).rstrip().splitlines():
        logging.info("%s", line)


def _json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _json_safe(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    try:
        import numpy as np

        if isinstance(value, np.integer):
            return int(value)
        if isinstance(value, np.floating):
            return float(value)
        if isinstance(value, np.ndarray):
            return value.tolist()
    except ImportError:
        pass
    return value


def write_json_summary(summary: dict, json_out: str | Path) -> None:
    json_path = Path(json_out)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    with json_path.open("w") as handle:
        json.dump(_json_safe(summary), handle, indent=2, sort_keys=False)
        handle.write("\n")


def installed_package_version(distribution: str, module_name: str | None = None) -> str | None:
    """Return package metadata or a source-install identity without a hard dependency."""
    try:
        return metadata.version(distribution)
    except metadata.PackageNotFoundError:
        try:
            module = import_module(module_name or distribution.replace("-", "_"))
        except ImportError:
            return None
        version = getattr(module, "__version__", None)
        if version is not None:
            return str(version)
        module_file = getattr(module, "__file__", None)
        if module_file is None:
            return "source:module_file_unavailable"
        source_path = Path(module_file).resolve()
        if source_path.suffix == ".pyc":
            source_candidate = source_path.with_suffix(".py")
            if source_candidate.is_file():
                source_path = source_candidate
        source_digest = hashlib.sha256(source_path.read_bytes()).hexdigest()
        return f"source:{module_name or distribution.replace('-', '_')}:{source_digest}"


def structure_ordering_metadata(structures: list[Path]) -> dict:
    """Describe the ordered collected-structure list consumed by AF2RAVE."""
    ordering_entries = getattr(structures, "ordering_entries", None)
    if ordering_entries is None:
        ordering_entries = [structure.name for structure in structures]
    digest = hashlib.sha256("\n".join(ordering_entries).encode("utf-8")).hexdigest()
    return {
        "ordering_rule": "source directory order, then natural sort of source-relative PDB paths",
        "structure_count": len(ordering_entries),
        "source_relative_path_sha256": digest,
        "ordering_entries_source": "in_memory_collection" if hasattr(structures, "ordering_entries") else "pdb_names",
    }


def runtime_versions() -> dict[str, str | None]:
    """Capture versions that can affect AF2RAVE feature analysis and clustering."""
    return {
        "python": platform.python_version(),
        "af2rave": installed_package_version("af2rave"),
        "numpy": installed_package_version("numpy"),
        "scikit_learn": installed_package_version("scikit-learn", module_name="sklearn"),
        "mdtraj": installed_package_version("mdtraj"),
        "biopython": installed_package_version("biopython", module_name="Bio"),
        "matplotlib": installed_package_version("matplotlib"),
    }


def callable_kwargs(function, config: dict, exclude: set[str] | None = None, extra_key: str | None = None) -> dict:
    exclude = exclude or set()
    signature = inspect.signature(function)
    parameters = signature.parameters
    accepts_extra = any(parameter.kind == inspect.Parameter.VAR_KEYWORD for parameter in parameters.values())
    accepted_names = {
        name
        for name, parameter in parameters.items()
        if parameter.kind in (inspect.Parameter.POSITIONAL_OR_KEYWORD, inspect.Parameter.KEYWORD_ONLY)
    }

    kwargs = {
        key: value
        for key, value in config.items()
        if key in accepted_names and key not in exclude and value is not None
    }

    if extra_key and config.get(extra_key):
        extra_kwargs = config[extra_key]
        if not isinstance(extra_kwargs, dict):
            raise ValueError(f"{extra_key} must be a dictionary of keyword arguments.")
        kwargs.update({key: value for key, value in extra_kwargs.items() if value is not None})

    if not accepts_extra:
        unexpected = set(kwargs) - accepted_names
        if unexpected:
            names = ", ".join(sorted(unexpected))
            raise ValueError(f"Unexpected arguments for {function.__name__}: {names}")

    return kwargs


def natural_sort_key(path: Path | str) -> list[object]:
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r"(\d+)", str(path))]


def resolve_structure_dirs(paths: dict) -> list[Path]:
    configured_dirs = []
    if paths.get("structures_dir"):
        configured_dirs.append(paths["structures_dir"])

    structure_dirs = paths.get("structure_dirs", [])
    if structure_dirs:
        if not isinstance(structure_dirs, list):
            raise ValueError("paths.structure_dirs must be a YAML list of directory paths.")
        configured_dirs.extend(structure_dirs)

    if paths.get("structure_dirs_file"):
        dirs_file = Path(paths["structure_dirs_file"]).expanduser().resolve()
        if not dirs_file.is_file():
            raise FileNotFoundError(f"Structure-directory list file not found: {dirs_file}")
        for line in dirs_file.read_text().splitlines():
            path_text = line.split("#", 1)[0].strip()
            if path_text:
                configured_dirs.append(path_text)

    if not configured_dirs:
        raise ValueError(
            "Provide paths.structures_dir, paths.structure_dirs, or paths.structure_dirs_file."
        )

    resolved_dirs = []
    seen = set()
    for value in configured_dirs:
        directory = Path(value).expanduser().resolve()
        if not directory.is_dir():
            raise NotADirectoryError(f"Predicted-structure directory not found: {directory}")
        if directory not in seen:
            seen.add(directory)
            resolved_dirs.append(directory)
    return resolved_dirs


def collect_structures(source_dirs: list[Path], output_dir: Path) -> list[Path]:
    structures_dir = ensure_dir(output_dir / "structures")
    discovered = []
    seen_files = set()
    for source_index, source_dir in enumerate(source_dirs, start=1):
        for pdb_file in sorted(
            source_dir.rglob("*.pdb"),
            key=lambda path: natural_sort_key(path.relative_to(source_dir)),
        ):
            resolved_pdb = pdb_file.resolve()
            if resolved_pdb not in seen_files:
                seen_files.add(resolved_pdb)
                discovered.append((source_index, source_dir, pdb_file))
    if not discovered:
        directories = ", ".join(str(path) for path in source_dirs)
        raise FileNotFoundError(f"No PDB files found in: {directories}")

    copied = []
    manifest_rows = []
    ordering_entries = []
    collected_targets = set()
    for collection_index, (source_index, source_dir, pdb_file) in enumerate(discovered, start=1):
        relative_path = pdb_file.relative_to(source_dir).as_posix()
        safe_stem = re.sub(r"[^A-Za-z0-9._-]+", "_", pdb_file.stem)[:80] or "structure"
        target_name = f"s{source_index:02d}_{collection_index:06d}__{safe_stem}.pdb"
        target = structures_dir / target_name
        if target in collected_targets:
            raise RuntimeError(f"Collected structure filename collision: {target}")
        collected_targets.add(target)
        if pdb_file.resolve() != target.resolve():
            shutil.copy2(pdb_file, target)
        copied.append(target)
        manifest_rows.append(
            {
                "collection_index": collection_index,
                "source_directory": str(source_dir),
                "source_pdb": str(pdb_file),
                "collected_pdb": str(target),
            }
        )
        ordering_entries.append(f"{source_index}\t{source_dir.name}\t{relative_path}")

    if len(set(copied)) != len(copied):
        raise RuntimeError("Collected structure paths are not unique.")

    with (output_dir / "structure_collection_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_rows[0].keys())
        writer.writeheader()
        writer.writerows(manifest_rows)
    return CollectedStructures(copied, ordering_entries)


def parse_plddt_residue_ranges(value) -> list[tuple[str | None, int, int]]:
    """Parse optional `CHAIN:START-END` pLDDT averaging ranges."""
    if value in (None, []):
        return []
    if not isinstance(value, list):
        raise ValueError("paths.plddt_residue_ranges must be a YAML list, such as ['A:41-187'].")

    ranges = []
    pattern = re.compile(r"(?:(?P<chain>[^:\s]+):)?(?P<start>-?\d+)\s*-\s*(?P<end>-?\d+)$")
    for item in value:
        match = pattern.fullmatch(str(item).strip())
        if match is None:
            raise ValueError(
                f"Invalid pLDDT residue range {item!r}. Use START-END or CHAIN:START-END."
            )
        start = int(match.group("start"))
        end = int(match.group("end"))
        if start > end:
            raise ValueError(f"Invalid pLDDT residue range {item!r}: start must not exceed end.")
        ranges.append((match.group("chain"), start, end))
    return ranges


def residue_is_in_plddt_ranges(chain_id: str, residue_number: int, ranges: list[tuple[str | None, int, int]]) -> bool:
    return not ranges or any(
        (range_chain is None or range_chain == chain_id) and start <= residue_number <= end
        for range_chain, start, end in ranges
    )


def pdb_mean_plddt(
    pdb_file: Path,
    bfactor_scale: str = "0_100",
    residue_ranges: list[tuple[str | None, int, int]] | None = None,
) -> float | None:
    """Return a residue-averaged pLDDT score normalized to the 0-100 scale."""
    if bfactor_scale not in PLDDT_BFACTOR_SCALES:
        raise ValueError("paths.plddt_bfactor_scale must be '0_100' or '0_1'.")
    lower_bound, upper_bound, multiplier = PLDDT_BFACTOR_SCALES[bfactor_scale]
    residue_ranges = residue_ranges or []
    residue_plddt = {}
    with pdb_file.open() as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if len(line) < 66:
                continue
            try:
                plddt = float(line[60:66])
            except ValueError:
                continue

            chain_id = line[21].strip()
            residue_number = line[22:26].strip()
            insertion_code = line[26].strip()
            try:
                numeric_residue_number = int(residue_number)
            except ValueError:
                continue
            if not residue_is_in_plddt_ranges(chain_id, numeric_residue_number, residue_ranges):
                continue
            if plddt < lower_bound or plddt > upper_bound:
                raise ValueError(
                    f"B-factor value {plddt:.3f} in {pdb_file} is outside the configured "
                    f"pLDDT {bfactor_scale} range [{lower_bound:g}, {upper_bound:g}]."
                )
            residue_key = (chain_id, residue_number, insertion_code)
            residue_plddt.setdefault(residue_key, plddt)

    if not residue_plddt:
        return None
    if bfactor_scale == "0_100" and max(residue_plddt.values()) <= 1.0:
        raise ValueError(
            f"B-factor values in {pdb_file} appear to use a 0-1 pLDDT scale. "
            "Set paths.plddt_bfactor_scale: '0_1' to use and normalize them."
        )
    return multiplier * sum(residue_plddt.values()) / len(residue_plddt)


def calculate_plddt_scores(
    structures: list[Path],
    bfactor_scale: str,
    residue_ranges: list[tuple[str | None, int, int]],
    skip_invalid_structures: bool = False,
) -> dict[Path, float | None]:
    if bfactor_scale not in PLDDT_BFACTOR_SCALES:
        raise ValueError("paths.plddt_bfactor_scale must be '0_100' or '0_1'.")
    scores = {}
    for structure in structures:
        try:
            scores[structure] = pdb_mean_plddt(
                structure, bfactor_scale=bfactor_scale, residue_ranges=residue_ranges
            )
        except ValueError as exc:
            if not skip_invalid_structures:
                raise
            logging.warning("Skipping invalid B-factor pLDDT score for %s: %s", structure, exc)
            scores[structure] = None
    return scores


def write_plddt_scores(output_dir: Path, plddt_scores: dict[Path, float | None]) -> None:
    with (output_dir / "plddt_scores.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["pdb", "residue_averaged_plddt"])
        for structure, score in sorted(
            plddt_scores.items(), key=lambda item: (item[1] is not None, item[1] or float("-inf")), reverse=True
        ):
            writer.writerow([structure, "" if score is None else f"{score:.4f}"])


def resolve_reference_structure(reference: str, structures: list[Path], source_dirs: list[Path]) -> Path:
    if not structures:
        raise ValueError("Cannot resolve a reference structure because no input PDBs were collected.")
    candidates = []
    reference_path = Path(reference).expanduser()
    if reference_path.is_absolute():
        candidates.append(reference_path)
    else:
        candidates.extend(
            [
                Path.cwd() / reference_path,
                structures[0].parent / reference_path,
                structures[0].parent / reference_path.name,
            ]
        )
        for source_dir in source_dirs:
            candidates.extend([source_dir / reference_path, source_dir / reference_path.name])

    for candidate in candidates:
        if candidate.exists():
            resolved = candidate.resolve()
            for structure in structures:
                if structure.resolve() == resolved or structure.name == resolved.name:
                    return structure
            return resolved

    raise FileNotFoundError(f"Reference structure not found: {reference}")


def select_reference_structure(
    paths: dict,
    structures: list[Path],
    source_dirs: list[Path],
    plddt_scores: dict[Path, float | None],
) -> Path | None:
    explicit_reference = paths.get("reference_structure")
    if explicit_reference:
        reference = resolve_reference_structure(str(explicit_reference), structures, source_dirs)
        logging.info("Using user-provided AF2RAVE reference structure: %s", reference)
        return reference

    reference_selection = str(paths.get("reference_selection", "best_plddt")).lower()
    if reference_selection in {"first", "first_loaded", "af2rave_default", "none"}:
        logging.info("Using AF2RAVE default reference structure: first loaded PDB.")
        return None
    if reference_selection not in {"best_plddt", "auto", "highest_plddt"}:
        raise ValueError("paths.reference_selection must be 'best_plddt' or 'first'.")

    scored = [(score, structure) for structure, score in plddt_scores.items() if score is not None]

    if not scored:
        logging.warning(
            "Could not read AF2-style pLDDT values from PDB B-factor columns. "
            "Using AF2RAVE default reference structure: first loaded PDB."
        )
        return None

    mean_plddt, reference = max(scored, key=lambda item: item[0])
    logging.info("Using highest residue-averaged pLDDT reference structure: %s (%.2f)", reference, mean_plddt)
    return reference


def import_af2rave_modules():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import af2rave.feature as af2feature
        import af2rave.simulation as af2sim
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        raise ImportError(
            "AF2RAVE analysis requires af2rave, numpy, and matplotlib in the active environment. "
            "Install af2rave in a separate environment and run this script from that environment."
        ) from exc
    return np, plt, af2feature, af2sim


def configured_steric_clash_cutoff(filtering: dict):
    """Return the canonical steric cutoff, legacy alias, or AF2RAVE default."""
    if "steric_clash_cutoff" in filtering:
        return filtering["steric_clash_cutoff"]
    if "nonbonded_cutoff" in filtering:
        return filtering["nonbonded_cutoff"]
    return 1.0


def write_feature_stats(rmsa, np, plt, output_dir: Path, filtering: dict, analysis: dict) -> None:
    rmsd_selection = filtering.get("rmsd_selection", "name CA")
    rmsd_cutoff = filtering.get("rmsd_cutoff")
    nonbonded_cutoff = configured_steric_clash_cutoff(filtering)
    peptide_bond_mean = filtering.get("peptide_bond_mean")
    peptide_bond_std = filtering.get("peptide_bond_std")
    write_full_stats = bool(analysis.get("write_full_feature_stats", False))
    write_nonbonded_stats = write_full_stats or nonbonded_cutoff is not None
    write_peptide_stats = (
        write_full_stats
        or peptide_bond_mean is not None
        or peptide_bond_std is not None
    )

    rmsd = rmsa.get_rmsd(rmsd_selection)
    nonbonded_distances = rmsa.minimum_nonbonded_distance if write_nonbonded_stats else {}
    pb_stats = rmsa.peptide_bond_stats if write_peptide_stats else {}

    rows = []
    for pdb_name in rmsa.pdb_name:
        rows.append(
            {
                "pdb": pdb_name,
                "rmsd": rmsd.get(pdb_name, ""),
                "minimum_nonbonded_distance": nonbonded_distances.get(pdb_name, ""),
                "peptide_bond_mean": pb_stats.get(pdb_name, ("", ""))[0],
                "peptide_bond_std": pb_stats.get(pdb_name, ("", ""))[1],
            }
        )

    fieldnames = [
        "pdb",
        "rmsd",
        "minimum_nonbonded_distance",
        "peptide_bond_mean",
        "peptide_bond_std",
    ]
    with (output_dir / "feature_stats.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    n_panels = 1 + int(write_nonbonded_stats) + (2 if write_peptide_stats else 0)
    fig, axes = plt.subplots(figsize=(6, 3 * n_panels), nrows=n_panels, ncols=1, tight_layout=True)
    axes = np.atleast_1d(axes)
    panel_index = 0

    axes[panel_index].set_title(f"RMSD: {rmsd_selection}")
    rmsd_values = np.array(list(rmsd.values()), dtype=float)
    axes[panel_index].hist(rmsd_values, bins=20)
    if rmsd_cutoff is not None:
        axes[panel_index].axvline(x=float(rmsd_cutoff), color="r", linestyle="--")
    panel_index += 1

    if write_nonbonded_stats:
        axes[panel_index].set_title("Nearest heavy atom distance")
        clash_values = np.array(list(nonbonded_distances.values()), dtype=float)
        axes[panel_index].hist(clash_values, bins=20)
        if nonbonded_cutoff is not None:
            axes[panel_index].axvline(x=float(nonbonded_cutoff), color="r", linestyle="--")
        panel_index += 1

    if write_peptide_stats:
        mean = np.array([pb_stats[name][0] for name in pb_stats], dtype=float)
        std = np.array([pb_stats[name][1] for name in pb_stats], dtype=float)
        axes[panel_index].set_title("Peptide bond length mean")
        axes[panel_index].hist(mean, bins=20)
        if peptide_bond_mean is not None:
            axes[panel_index].axvline(x=float(peptide_bond_mean), color="r", linestyle="--")
        panel_index += 1

        axes[panel_index].set_title("Peptide bond length std")
        axes[panel_index].hist(std, bins=20)
        if peptide_bond_std is not None:
            axes[panel_index].axvline(x=float(peptide_bond_std), color="r", linestyle="--")

    fig.savefig(output_dir / "feature_stats.png", dpi=160)
    plt.close(fig)

    skipped_stats = []
    if not write_nonbonded_stats:
        skipped_stats.append("minimum_nonbonded_distance")
    if not write_peptide_stats:
        skipped_stats.extend(["peptide_bond_mean", "peptide_bond_std"])
    if skipped_stats:
        logging.info(
            "Skipped optional feature stats not used for filtering: %s. "
            "Set analysis.write_full_feature_stats: true to calculate them.",
            ", ".join(skipped_stats),
        )


def validate_filter_parameter_precedence(
    filtering: dict, extra_key: str, canonical_parameters: dict[str, str]
) -> None:
    """Reject alternate filter parameter names that would disagree with canonical YAML keys."""
    direct_conflicts = [
        parameter
        for parameter, config_key in canonical_parameters.items()
        if parameter != config_key and parameter in filtering
    ]
    if direct_conflicts:
        details = ", ".join(
            f"filtering.{parameter} (use filtering.{canonical_parameters[parameter]})"
            for parameter in direct_conflicts
        )
        raise ValueError(f"Unsupported direct filter parameter(s): {details}.")

    extra_kwargs = filtering.get(extra_key, {})
    if extra_kwargs is None:
        return
    if not isinstance(extra_kwargs, dict):
        raise ValueError(f"filtering.{extra_key} must be a dictionary of keyword arguments.")
    extra_conflicts = sorted(set(extra_kwargs) & set(canonical_parameters))
    if extra_conflicts:
        details = ", ".join(
            f"{parameter} (use filtering.{canonical_parameters[parameter]})"
            for parameter in extra_conflicts
        )
        raise ValueError(f"filtering.{extra_key} must not override canonical filter parameters: {details}.")


def apply_filters(rmsa, filtering: dict) -> None:
    filters = []

    rmsd_cutoff = filtering.get("rmsd_cutoff")
    if rmsd_cutoff is not None:
        validate_filter_parameter_precedence(
            filtering,
            "rmsd_filter_kwargs",
            {"selection": "rmsd_selection", "rmsd_cutoff": "rmsd_cutoff"},
        )
        rmsd_kwargs = callable_kwargs(
            rmsa.rmsd_filter,
            filtering,
            exclude={
                "selection",
                "rmsd_selection",
                "rmsd_cutoff",
                "steric_clash_cutoff",
                "nonbonded_cutoff",
                "peptide_bond_mean",
                "peptide_bond_std",
                "rmsd_filter_kwargs",
                "steric_clash_filter_kwargs",
                "peptide_bond_filter_kwargs",
            },
            extra_key="rmsd_filter_kwargs",
        )
        rmsd_kwargs.update(
            {"selection": filtering.get("rmsd_selection", "name CA"), "rmsd_cutoff": float(rmsd_cutoff)}
        )
        filters.append(rmsa.rmsd_filter(**rmsd_kwargs))

    nonbonded_cutoff = configured_steric_clash_cutoff(filtering)
    if nonbonded_cutoff is not None:
        validate_filter_parameter_precedence(
            filtering,
            "steric_clash_filter_kwargs",
            {"min_non_bonded_cutoff": "steric_clash_cutoff"},
        )
        steric_kwargs = callable_kwargs(
            rmsa.steric_clash_filter,
            filtering,
            exclude={
                "min_non_bonded_cutoff",
                "rmsd_selection",
                "rmsd_cutoff",
                "steric_clash_cutoff",
                "nonbonded_cutoff",
                "peptide_bond_mean",
                "peptide_bond_std",
                "rmsd_filter_kwargs",
                "steric_clash_filter_kwargs",
                "peptide_bond_filter_kwargs",
            },
            extra_key="steric_clash_filter_kwargs",
        )
        steric_kwargs["min_non_bonded_cutoff"] = float(nonbonded_cutoff)
        filters.append(rmsa.steric_clash_filter(**steric_kwargs))

    peptide_bond_mean = filtering.get("peptide_bond_mean")
    peptide_bond_std = filtering.get("peptide_bond_std")
    if peptide_bond_mean is not None or peptide_bond_std is not None:
        if peptide_bond_mean is None or peptide_bond_std is None:
            raise ValueError(
                "Set both peptide_bond_mean and peptide_bond_std to enable peptide-bond filtering."
            )
        validate_filter_parameter_precedence(
            filtering,
            "peptide_bond_filter_kwargs",
            {"mean_cutoff": "peptide_bond_mean", "std_cutoff": "peptide_bond_std"},
        )
        peptide_kwargs = callable_kwargs(
            rmsa.peptide_bond_filter,
            filtering,
            exclude={
                "mean_cutoff",
                "std_cutoff",
                "rmsd_selection",
                "rmsd_cutoff",
                "steric_clash_cutoff",
                "nonbonded_cutoff",
                "peptide_bond_mean",
                "peptide_bond_std",
                "rmsd_filter_kwargs",
                "steric_clash_filter_kwargs",
                "peptide_bond_filter_kwargs",
            },
            extra_key="peptide_bond_filter_kwargs",
        )
        peptide_kwargs.update({"mean_cutoff": float(peptide_bond_mean), "std_cutoff": float(peptide_bond_std)})
        filters.append(rmsa.peptide_bond_filter(**peptide_kwargs))

    if filters:
        rmsa.apply_filter(*filters)


def resolve_output_path(path_value: str | Path | None, output_dir: Path) -> Path | None:
    if path_value is None:
        return None
    path = Path(path_value).expanduser()
    if not path.is_absolute():
        path = output_dir / path
    return path


def feature_cache_path(output_dir: Path, feature_config: dict) -> Path:
    return resolve_output_path(feature_config.get("feature_cache_npz", "ranked_feature_values.npz"), output_dir)


def cluster_centers_path(output_dir: Path, clustering: dict) -> Path:
    configured_path = clustering.get("cluster_centers_dir", clustering.get("centers_dir", "cluster_center_structures"))
    return resolve_output_path(configured_path, output_dir)


def write_selected_features_csv(output_dir: Path, names: list[str], cv: list[float]) -> None:
    with (output_dir / "selected_features.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["rank", "feature_name", "coefficient_of_variation"])
        for rank, (name, value) in enumerate(zip(names, cv), start=1):
            writer.writerow([rank, name, value])


def feature_cache_signature(selection: str, rank_kwargs: dict) -> str:
    """Fingerprint the configuration that defines cached AF2RAVE features."""
    payload = json.dumps(
        _json_safe({"selection": selection, "rank_feature_kwargs": rank_kwargs}),
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def write_feature_cache(
    rmsa,
    np,
    cache_path: Path,
    names: list[str],
    cv: list[float],
    selection: str,
    rank_kwargs: dict,
) -> None:
    ensure_dir(cache_path.parent)
    feature_matrix = np.array([rmsa.features[name] for name in names], dtype=np.float64).T
    atom_pairs = np.array([rmsa.atom_pairs[name] for name in names], dtype=np.int_)
    np.savez(
        cache_path,
        names=np.array(names, dtype=object),
        coefficient_of_variation=np.array(cv, dtype=np.float64),
        atom_pairs=atom_pairs,
        feature_matrix=feature_matrix,
        pdb_names=np.array([str(name) for name in rmsa.pdb_name], dtype=object),
        feature_selection=np.array(selection),
        rank_feature_kwargs_json=np.array(json.dumps(_json_safe(rank_kwargs), sort_keys=True)),
        feature_definition_sha256=np.array(feature_cache_signature(selection, rank_kwargs)),
        af2rave_version=np.array(installed_package_version("af2rave") or ""),
    )
    logging.info("Wrote ranked feature cache: %s", cache_path)


def selected_feature_cv_cutoff(cv: list[float], n_features: int) -> float:
    """Validate AF2RAVE ranking order and return the last selected feature's CV."""
    if n_features < 1:
        raise ValueError("features.n_features must be at least 1.")
    if not cv:
        raise ValueError("AF2RAVE rank_feature returned no feature CV values.")
    if any(following > current + 1e-12 for current, following in zip(cv, cv[1:])):
        raise ValueError("AF2RAVE rank_feature returned CV values that are not non-increasing.")
    return float(cv[min(n_features, len(cv)) - 1])


def load_feature_cache(
    rmsa, np, cache_path: Path, selection: str, rank_kwargs: dict
) -> tuple[list[str], list[float]]:
    with np.load(cache_path, allow_pickle=True) as cache:
        required_metadata = {
            "feature_selection",
            "rank_feature_kwargs_json",
            "feature_definition_sha256",
            "af2rave_version",
        }
        missing_metadata = sorted(required_metadata - set(cache.files))
        if missing_metadata:
            raise ValueError(
                f"Feature cache lacks compatibility metadata ({', '.join(missing_metadata)}): {cache_path}. "
                "Regenerate it with reuse_feature_cache: false."
            )
        names = [str(name) for name in cache["names"].tolist()]
        cv = cache["coefficient_of_variation"].astype(np.float64).tolist()
        atom_pairs = cache["atom_pairs"].astype(np.int_)
        feature_matrix = cache["feature_matrix"].astype(np.float64)
        cached_pdb_names = [str(name) for name in cache["pdb_names"].tolist()] if "pdb_names" in cache else None
        cached_selection = str(cache["feature_selection"].item())
        cached_rank_kwargs = str(cache["rank_feature_kwargs_json"].item())
        cached_signature = str(cache["feature_definition_sha256"].item())
        cached_af2rave_version = str(cache["af2rave_version"].item()) or None

    if feature_matrix.shape[0] != len(rmsa):
        raise ValueError(
            f"Feature cache structure count ({feature_matrix.shape[0]}) does not match "
            f"the current filtered structure count ({len(rmsa)}): {cache_path}"
        )
    current_pdb_names = [str(name) for name in rmsa.pdb_name]
    if cached_pdb_names is not None and cached_pdb_names != current_pdb_names:
        raise ValueError(
            "Feature cache PDB list/order does not match the current filtered structures. "
            f"Recalculate the cache or use the same filtered structure set: {cache_path}"
        )
    current_signature = feature_cache_signature(selection, rank_kwargs)
    if cached_signature != current_signature:
        raise ValueError(
            "Feature cache selection or rank_feature_kwargs do not match the current configuration. "
            f"Cached selection: {cached_selection!r}; cached rank_feature_kwargs: {cached_rank_kwargs}. "
            "Regenerate the cache with reuse_feature_cache: false."
        )
    current_af2rave_version = installed_package_version("af2rave")
    if cached_af2rave_version != current_af2rave_version:
        raise ValueError(
            "Feature cache AF2RAVE version does not match the current environment "
            f"({cached_af2rave_version or 'unknown'} != {current_af2rave_version or 'unknown'}). "
            "Regenerate the cache with reuse_feature_cache: false."
        )
    if feature_matrix.shape[1] != len(names):
        raise ValueError(
            f"Feature cache column count ({feature_matrix.shape[1]}) does not match "
            f"the number of cached feature names ({len(names)}): {cache_path}"
        )
    if atom_pairs.shape[0] != len(names):
        raise ValueError(
            f"Feature cache atom-pair count ({atom_pairs.shape[0]}) does not match "
            f"the number of cached feature names ({len(names)}): {cache_path}"
        )

    rmsa.features.clear()
    rmsa.atom_pairs.clear()
    for index, name in enumerate(names):
        rmsa.features[name] = feature_matrix[:, index]
        rmsa.atom_pairs[name] = atom_pairs[index]

    logging.info("Loaded ranked feature cache: %s", cache_path)
    return names, cv


def rank_features(
    rmsa, np, plt, output_dir: Path, feature_config: dict
) -> tuple[list[str], list[float], Path | None, Path | None]:
    selection = feature_config.get("selection", "name CA")
    n_features = int(feature_config.get("n_features", 200))
    cache_path = feature_cache_path(output_dir, feature_config)
    reuse_cache = feature_config.get("reuse_feature_cache", False)
    rank_kwargs = callable_kwargs(
        rmsa.rank_feature,
        feature_config,
        exclude={
            "selection",
            "n_features",
            "rank_feature_kwargs",
            "reuse_feature_cache",
            "write_feature_cache",
            "feature_cache_npz",
        },
        extra_key="rank_feature_kwargs",
    )
    loaded_from_cache = False
    written_cache = None
    reused_cache = None
    if reuse_cache and cache_path.exists():
        names, cv = load_feature_cache(rmsa, np, cache_path, selection, rank_kwargs)
        loaded_from_cache = True
        reused_cache = cache_path
    else:
        if reuse_cache:
            logging.info("Feature cache not found; calculating features: %s", cache_path)

        names, cv = rmsa.rank_feature(selection=selection, **rank_kwargs)
        names = list(names)
        cv = list(cv)

    if len(names) < n_features:
        logging.warning(
            "Requested %d ranked features, but AF2RAVE returned %d; clustering will use all %d available features.",
            n_features,
            len(names),
            len(names),
        )
    cv_cutoff = selected_feature_cv_cutoff(cv, n_features)
    write_selected_features_csv(output_dir, names, cv)
    if feature_config.get("write_feature_cache", False) and not loaded_from_cache:
        write_feature_cache(rmsa, np, cache_path, names, cv, selection, rank_kwargs)
        written_cache = cache_path

    fig, ax = plt.subplots(figsize=(4, 3), tight_layout=True)
    counts, _, _ = ax.hist(cv, bins=100, label=f"{len(names)}-dim", density=True)
    ax.vlines(
        cv_cutoff,
        0,
        np.max(counts) * 0.6,
        color="red",
        linestyle="--",
        linewidth=1,
        label=f"Cutoff={cv_cutoff:.3f}",
    )
    ax.set_ylim(0, np.max(counts) * 1.2)
    ax.set_xlabel("Coefficient of Variation")
    ax.set_ylabel("Density")
    ax.legend()
    fig.savefig(output_dir / "feature_cv_histogram.png", dpi=160)
    plt.close(fig)

    return list(names), list(cv), written_cache, reused_cache


def selected_feature_matrix(rmsa, np, selected_names: list[str]):
    if not selected_names:
        raise ValueError("Feature list cannot be empty.")
    return np.array([rmsa.features[name] for name in selected_names], dtype=np.float64).T


def kmeans_clustering(rmsa, np, output_dir: Path, selected_names: list[str], clustering: dict):
    from sklearn.cluster import KMeans

    z = selected_feature_matrix(rmsa, np, selected_names)
    npoints, ndim = z.shape
    n_clusters = int(
        clustering.get(
            "n_clusters",
            clustering.get("k", clustering.get("max_centers", 10)),
        )
    )
    if n_clusters < 1:
        raise ValueError("clustering.n_clusters must be at least 1 for k-means clustering.")
    if n_clusters > npoints:
        raise ValueError(
            f"clustering.n_clusters ({n_clusters}) cannot exceed the number of filtered structures ({npoints})."
        )

    kmeans_kwargs = {
        "n_clusters": n_clusters,
        "random_state": int(clustering.get("random_state", clustering.get("randomseed", 0))),
        "n_init": clustering.get("n_init", 10),
    }
    kmeans_kwargs.update(clustering.get("kmeans_kwargs", {}))
    model = KMeans(**kmeans_kwargs)
    labels = model.fit_predict(z)

    center_ids = []
    for cluster_id in range(n_clusters):
        members = np.where(labels == cluster_id)[0]
        if members.size == 0:
            continue
        distances = np.linalg.norm(z[members] - model.cluster_centers_[cluster_id], axis=1)
        best_member = int(members[int(np.argmin(distances))])
        center_ids.append(best_member)

    center_id = np.array(center_ids, dtype=np.int_)
    if len(center_ids) != n_clusters:
        logging.warning(
            "K-means requested %d clusters, but selected %d non-empty cluster representatives.",
            n_clusters,
            len(center_ids),
        )

    center_id_set = set(map(int, center_id))
    with (output_dir / "cluster_assignments.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "pdb",
                "cluster",
                "euclidean_distance_to_centroid",
                "rms_feature_distance_to_centroid",
                "is_center",
            ]
        )
        for index, (pdb_name, label) in enumerate(zip(rmsa.pdb_name, labels)):
            euclidean_distance = np.linalg.norm(z[index] - model.cluster_centers_[label])
            rms_feature_distance = euclidean_distance / np.sqrt(ndim)
            writer.writerow(
                [
                    pdb_name,
                    int(label),
                    float(euclidean_distance),
                    float(rms_feature_distance),
                    index in center_id_set,
                ]
            )

    return center_id


def regular_space_clustering(rmsa, selected_names: list[str], clustering: dict, output_dir: Path):
    initial_min_dist = float(clustering.get("min_dist", clustering.get("dist_cutoff", 5.0)))
    max_centers = int(clustering.get("max_centers", 100))
    target_centers = clustering.get("target_centers")
    max_iterations = int(clustering.get("target_search_max_iterations", 30))
    min_dist_tolerance = float(clustering.get("target_search_min_dist_tolerance", 0.05))
    relative_tolerance = float(clustering.get("target_search_relative_tolerance", 0.01))
    max_evaluations = int(clustering.get("target_search_max_evaluations", 1000))

    if initial_min_dist <= 0:
        raise ValueError("clustering.min_dist must be greater than zero for regular-space clustering.")
    if max_centers < 2:
        raise ValueError("clustering.max_centers must be at least 2 for regular-space clustering.")

    base_kwargs = callable_kwargs(
        rmsa.regular_space_clustering,
        clustering,
        exclude={
            "method",
            "feature_name",
            "min_dist",
            "dist_cutoff",
            "max_centers",
            "target_centers",
            "target_search_max_iterations",
            "target_search_min_dist_tolerance",
            "target_search_relative_tolerance",
            "target_search_max_evaluations",
            "n_clusters",
            "k",
            "random_state",
            "n_init",
            "kmeans_kwargs",
            "pca_components",
            "pca_kwargs",
        },
    )

    def run_once(min_dist: float):
        try:
            center_id = rmsa.regular_space_clustering(
                selected_names,
                min_dist=min_dist,
                max_centers=max_centers,
                **base_kwargs,
            )
        except ValueError as exc:
            message = str(exc)
            if "maximum number of cluster centers" not in message.lower():
                raise
            return None, max_centers, "capacity_exceeded"
        return center_id, len(center_id), "completed"

    if target_centers is None:
        return rmsa.regular_space_clustering(
            selected_names,
            min_dist=initial_min_dist,
            max_centers=max_centers,
            **base_kwargs,
        )

    target_centers = int(target_centers)
    if target_centers < 1:
        raise ValueError("clustering.target_centers must be at least 1.")
    if target_centers >= max_centers:
        raise ValueError(
            "clustering.target_centers must be smaller than clustering.max_centers. "
            "Keep max_centers as a higher AF2RAVE search capacity, such as 100."
        )
    if max_iterations < 1:
        raise ValueError("clustering.target_search_max_iterations must be at least 1.")
    if min_dist_tolerance <= 0:
        raise ValueError("clustering.target_search_min_dist_tolerance must be greater than zero.")
    if relative_tolerance <= 0:
        raise ValueError("clustering.target_search_relative_tolerance must be greater than zero.")
    if max_evaluations < 1:
        raise ValueError("clustering.target_search_max_evaluations must be at least 1.")

    trials = []

    def evaluate(min_dist: float, log_trial: bool = True):
        center_id, n_centers, status = run_once(min_dist)
        trials.append(
            {
                "trial": len(trials) + 1,
                "min_dist": min_dist,
                "n_centers": n_centers,
                "status": status,
            }
        )
        if log_trial or len(trials) % 1000 == 0:
            logging.info(
                "Regular-space target search: min_dist=%.6f, centers=%s, status=%s.",
                min_dist,
                n_centers if status == "completed" else f">={max_centers}",
                status,
            )
        return center_id, n_centers

    def write_trials():
        trial_path = output_dir / "regular_space_target_search.csv"
        with trial_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["trial", "min_dist", "n_centers", "status"])
            writer.writeheader()
            writer.writerows(trials)
        return trial_path

    if target_centers == 1:
        trials.append(
            {
                "trial": 1,
                "min_dist": "",
                "n_centers": 1,
                "status": "single_center_first_structure",
            }
        )
        logging.info(
            "Regular-space target search selected the first filtered structure for target_centers=1."
        )
        write_trials()
        return [0]

    def effective_tolerance(distance: float) -> float:
        return max(min_dist_tolerance, relative_tolerance * distance)

    def scan_highest_target(
        target_dist: float,
        upper_dist: float,
        target_center_id,
    ):
        """Find the largest target hit on a deterministic min_dist grid.

        Regular-space clustering is a greedy, order-sensitive procedure. Its
        center count can make small non-monotonic jumps, so binary refinement
        alone cannot safely identify the upper target boundary.
        """
        scan_tolerance = effective_tolerance(target_dist)
        steps = math.ceil((upper_dist - target_dist) / scan_tolerance)
        if steps > max_evaluations:
            trial_path = write_trials()
            raise ValueError(
                f"Finding the largest min_dist for target_centers={target_centers} requires "
                f"{steps} evaluations at effective tolerance={scan_tolerance:g}, "
                f"which exceeds target_search_max_evaluations={max_evaluations}. "
                f"Increase the evaluation limit or use a coarser tolerance. Trials: {trial_path}"
            )

        logging.info(
            "Regular-space target search scanning %d min_dist values from %.6f to %.6f "
            "at effective resolution %.6f.",
            steps,
            target_dist,
            upper_dist,
            scan_tolerance,
        )
        selected_dist = target_dist
        selected_center_id = target_center_id
        for step in range(1, steps):
            candidate_dist = target_dist + step * scan_tolerance
            center_id, n_centers = evaluate(candidate_dist, log_trial=False)
            if n_centers == target_centers:
                selected_dist = candidate_dist
                selected_center_id = center_id

        logging.info(
            "Regular-space target search selected largest sampled min_dist=%.6f for %d centers.",
            selected_dist,
            target_centers,
        )
        write_trials()
        return selected_center_id

    def find_upper_target_boundary(target_dist: float, target_center_id):
        """Bracket a target hit from above before refining its upper edge."""
        candidate_dist = target_dist
        for _ in range(max_iterations):
            candidate_dist *= 2.0
            center_id, n_centers = evaluate(candidate_dist)
            if n_centers == target_centers:
                target_dist = candidate_dist
                target_center_id = center_id
                continue
            # Greedy regular-space counts can jump in either direction. Any
            # non-target count brackets the highest known target hit for a scan.
            return scan_highest_target(target_dist, candidate_dist, target_center_id)

        logging.info(
            "Regular-space target search retained min_dist=%.6f for %d centers after %d expansions; "
            "no upper non-target boundary was observed.",
            target_dist,
            target_centers,
            max_iterations,
        )
        write_trials()
        return target_center_id

    center_id, n_centers = evaluate(initial_min_dist)
    if n_centers == target_centers:
        return find_upper_target_boundary(initial_min_dist, center_id)

    if n_centers > target_centers:
        lower_dist = initial_min_dist
        upper_dist = initial_min_dist
        for _ in range(max_iterations):
            upper_dist *= 2.0
            center_id, n_centers = evaluate(upper_dist)
            if n_centers == target_centers:
                return find_upper_target_boundary(upper_dist, center_id)
            if n_centers < target_centers:
                break
            lower_dist = upper_dist
        else:
            trial_path = write_trials()
            raise ValueError(
                f"Could not reduce regular-space centers to target_centers={target_centers} within "
                f"{max_iterations} min_dist expansions. Trials: {trial_path}"
            )
    else:
        upper_dist = initial_min_dist
        lower_dist = initial_min_dist
        for _ in range(max_iterations):
            lower_dist /= 2.0
            center_id, n_centers = evaluate(lower_dist)
            if n_centers == target_centers:
                return find_upper_target_boundary(lower_dist, center_id)
            if n_centers > target_centers:
                break
            upper_dist = lower_dist
        else:
            trial_path = write_trials()
            raise ValueError(
                f"Could not increase regular-space centers to target_centers={target_centers} within "
                f"{max_iterations} min_dist reductions. Trials: {trial_path}"
            )

    for _ in range(max_iterations):
        if upper_dist - lower_dist <= effective_tolerance(max(abs(lower_dist), abs(upper_dist))):
            break
        candidate_dist = (lower_dist + upper_dist) / 2.0
        center_id, n_centers = evaluate(candidate_dist)
        if n_centers == target_centers:
            return find_upper_target_boundary(candidate_dist, center_id)
        if n_centers > target_centers:
            lower_dist = candidate_dist
        else:
            upper_dist = candidate_dist

    trial_path = write_trials()
    raise ValueError(
        f"No min_dist found that produces exactly target_centers={target_centers}. "
        f"Searched {len(trials)} deterministic trials between {lower_dist:.6f} and {upper_dist:.6f}. "
        f"Review {trial_path} or change target_centers."
    )


def ca_atoms_by_residue(structure):
    """Return one protein CA atom and residue name per chain/residue from the first PDB model."""
    atoms = {}
    model = next(structure.get_models(), None)
    if model is None:
        return atoms
    for chain in model:
        for residue in chain:
            if residue.id[0] != " " or "CA" not in residue:
                continue
            atom = residue["CA"]
            if atom.is_disordered():
                atom = atom.selected_child
            atoms[(chain.id, residue.id)] = (atom, residue.get_resname())
    return atoms


def validate_common_ca_residue_names(moving_ca: dict, reference_ca: dict, common_residues: list) -> None:
    """Reject number-matched CA pairs whose residue identities differ."""
    mismatches = [
        (key, moving_ca[key][1], reference_ca[key][1])
        for key in common_residues
        if moving_ca[key][1] != reference_ca[key][1]
    ]
    if mismatches:
        details = ", ".join(
            f"{chain}:{residue[1]}{residue[2].strip() or ''} ({moving} != {reference})"
            for (chain, residue), moving, reference in mismatches[:5]
        )
        raise ValueError(
            "Cannot align PDBs with residue-name mismatches at shared chain/residue identifiers: "
            f"{details}."
        )


def pdb_ca_atom_count(pdb_file: Path) -> int:
    """Return the number of protein CA atoms in the first PDB model."""
    try:
        from Bio.PDB import PDBParser
    except ImportError as exc:
        raise ImportError(
            "Aligned cluster-center export requires Biopython. Install it with 'pip install biopython'."
        ) from exc

    structure = PDBParser(QUIET=True).get_structure("reference", str(pdb_file))
    return len(ca_atoms_by_residue(structure))


def align_pdb_to_reference(moving_pdb: Path, reference_pdb: Path, aligned_pdb: Path) -> tuple[int, float]:
    """Rigidly align a PDB to a reference with Biopython common-CA fitting."""
    try:
        from Bio.PDB import PDBIO, PDBParser, Superimposer
    except ImportError as exc:
        raise ImportError(
            "Aligned cluster-center export requires Biopython. Install it with 'pip install biopython'."
        ) from exc

    parser = PDBParser(QUIET=True)
    moving_structure = parser.get_structure("moving", str(moving_pdb))
    reference_structure = parser.get_structure("reference", str(reference_pdb))
    moving_ca = ca_atoms_by_residue(moving_structure)
    reference_ca = ca_atoms_by_residue(reference_structure)
    common_residues = sorted(set(moving_ca) & set(reference_ca))
    if len(common_residues) < 3:
        raise ValueError(
            f"Cannot align {moving_pdb.name} to {reference_pdb.name}: "
            f"only {len(common_residues)} common CA atoms were found."
        )
    validate_common_ca_residue_names(moving_ca, reference_ca, common_residues)

    superimposer = Superimposer()
    superimposer.set_atoms(
        [reference_ca[key][0] for key in common_residues],
        [moving_ca[key][0] for key in common_residues],
    )
    superimposer.apply(moving_structure.get_atoms())
    writer = PDBIO()
    writer.set_structure(moving_structure)
    writer.save(str(aligned_pdb))

    return len(common_residues), float(superimposer.rms)


def copy_cluster_centers(
    rmsa, center_id, output_dir: Path, clustering: dict, pipeline_reference: Path
) -> Path | None:
    enabled = clustering.get("copy_cluster_centers", clustering.get("copy_centers", False))
    if not enabled:
        return None

    centers_dir = cluster_centers_path(output_dir, clustering)
    ensure_dir(centers_dir)

    if clustering.get("clean_cluster_centers_dir", False):
        for old_pdb in centers_dir.glob("*.pdb"):
            old_pdb.unlink()

    align_centers = bool(clustering.get("align_cluster_centers", True))
    alignment_mode = str(clustering.get("center_alignment_reference", "pipeline_reference")).lower()
    if alignment_mode == "pipeline_reference":
        reference_source = Path(pipeline_reference)
    elif alignment_mode == "first_center":
        if len(center_id) == 0:
            raise ValueError("Cannot use first_center alignment with no selected cluster centers.")
        reference_source = Path(rmsa.pdb_name[int(center_id[0])])
    else:
        raise ValueError(
            "clustering.center_alignment_reference must be 'pipeline_reference' or 'first_center'."
        )

    manifest = centers_dir / "cluster_center_structures.csv"
    center_targets = {}
    with manifest.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "center_rank",
                "structure_index",
                "pdb",
                "source_path",
                "copied_path",
                "alignment_reference",
                "alignment_ca_atoms",
                "alignment_common_ca_rmsd_A",
            ]
        )
        for rank, index in enumerate(center_id, start=1):
            source = Path(rmsa.pdb_name[int(index)])
            target = centers_dir / source.name
            previous_source = center_targets.get(target)
            if previous_source is not None:
                raise RuntimeError(
                    f"Cluster-center filename collision: {source} and {previous_source} both map to {target}."
                )
            center_targets[target] = source
            if align_centers and source != reference_source:
                alignment_atoms, alignment_rmsd = align_pdb_to_reference(source, reference_source, target)
            else:
                shutil.copy2(source, target)
                alignment_atoms, alignment_rmsd = (
                    (pdb_ca_atom_count(reference_source), 0.0) if align_centers else ("", "")
                )
            writer.writerow(
                [
                    rank,
                    int(index),
                    source.name,
                    str(source),
                    str(target),
                    str(reference_source) if align_centers else "",
                    alignment_atoms,
                    f"{alignment_rmsd:.6f}" if isinstance(alignment_rmsd, float) else alignment_rmsd,
                ]
            )

    action = "Aligned and copied" if align_centers else "Copied"
    logging.info("%s %d cluster center structure(s) to %s.", action, len(center_id), centers_dir)
    return centers_dir


def cluster_and_plot(
    rmsa, np, plt, output_dir: Path, names: list[str], config: dict, pipeline_reference: Path
) -> object:
    feature_config = config.get("features", {})
    clustering = config.get("clustering", {})
    filtering = config.get("filtering", {})
    n_features = int(feature_config.get("n_features", 200))
    selected_names = names[:n_features]
    components = int(clustering.get("pca_components", clustering.get("n_components", 2)))
    if components < 2:
        raise ValueError("clustering.pca_components must be at least 2 for the PCA cluster plot.")

    method = str(clustering.get("method", "regular_space")).lower().replace("-", "_")
    if method in {"regular_space", "regular", "rsc"}:
        center_id = regular_space_clustering(rmsa, selected_names, clustering, output_dir)
    elif method in {"kmeans", "k_means"}:
        center_id = kmeans_clustering(rmsa, np, output_dir, selected_names, clustering)
    else:
        raise ValueError("clustering.method must be 'regular_space' or 'kmeans'.")

    with (output_dir / "cluster_centers.txt").open("w") as handle:
        for index in center_id:
            handle.write(f"{rmsa.pdb_name[index]}\n")

    copy_cluster_centers(rmsa, center_id, output_dir, clustering, pipeline_reference)

    feature_matrix = selected_feature_matrix(rmsa, np, selected_names)
    n_samples, n_selected_features = feature_matrix.shape
    pca_plot = output_dir / "pca_clusters.png"
    pca_diagnostic = output_dir / "pca_unavailable.txt"
    if min(n_samples, n_selected_features) < components:
        limitations = []
        if n_samples < components:
            limitations.append(f"the filtered structure count is {n_samples}")
        if n_selected_features < components:
            limitations.append(f"the selected feature count is {n_selected_features}")
        reason = " and ".join(limitations)
        pca_diagnostic.write_text(
            "PCA cluster plot was not created.\n\n"
            f"Requested PCA components: {components}\n"
            f"Filtered structures (samples): {n_samples}\n"
            f"Selected features: {n_selected_features}\n\n"
            f"Reason: A {components}-component PCA requires at least {components} filtered "
            f"structures and {components} selected features, but {reason}.\n\n"
            "Clustering completed, and downstream simulation-box setup can continue. "
            "Provide more input structures, relax filtering cutoffs, or select more features "
            "to create a 2D PCA plot.\n"
        )
        pca_plot.unlink(missing_ok=True)
        logging.warning("PCA cluster plot skipped: %s. See %s", reason, pca_diagnostic)
        return center_id
    pca_diagnostic.unlink(missing_ok=True)

    from sklearn.decomposition import PCA

    pca = PCA(
        n_components=components,
        **callable_kwargs(
            PCA,
            clustering,
            exclude={
                "method",
                "n_components",
                "pca_components",
                "dist_cutoff",
                "min_dist",
                "max_centers",
                "target_centers",
                "target_search_max_iterations",
                "target_search_min_dist_tolerance",
                "target_search_max_evaluations",
                "batch_size",
                "randomseed",
                "random_state",
                "n_clusters",
                "n_init",
                "kmeans_kwargs",
                "copy_cluster_centers",
                "copy_centers",
                "cluster_centers_dir",
                "centers_dir",
                "clean_cluster_centers_dir",
                "align_cluster_centers",
                "center_alignment_reference",
                "pca_kwargs",
            },
            extra_key="pca_kwargs",
        ),
    )
    pca_result = pca.fit_transform(feature_matrix)
    fig, ax = plt.subplots(figsize=(4, 3), layout="compressed")
    rmsd_by_pdb = rmsa.get_rmsd(filtering.get("rmsd_selection", "name CA"))
    missing_rmsd = [pdb_name for pdb_name in rmsa.pdb_name if pdb_name not in rmsd_by_pdb]
    if missing_rmsd:
        raise ValueError(
            "AF2RAVE RMSD results are missing entries for PCA plotting: "
            + ", ".join(map(str, missing_rmsd[:5]))
        )
    rmsd = np.array([rmsd_by_pdb[pdb_name] for pdb_name in rmsa.pdb_name], dtype=float)
    rmsd_vmax = filtering.get("rmsd_cutoff")
    if rmsd_vmax is None:
        rmsd_vmax = float(np.nanmax(rmsd)) if rmsd.size else 10.0
    im = ax.scatter(
        pca_result[:, 0],
        pca_result[:, 1],
        s=5,
        c=rmsd,
        cmap="gnuplot",
        norm="linear",
        edgecolors="none",
        vmin=0.0,
        vmax=float(rmsd_vmax),
    )
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    cbar = fig.colorbar(im, orientation="vertical")
    cbar.set_label("RMSD from reference structure")

    pca_centers = pca_result[center_id]
    ax.scatter(
        pca_centers[:, 0],
        pca_centers[:, 1],
        s=20,
        edgecolor="k",
        marker="H",
        facecolor="None",
        linewidth=0.5,
    )
    fig.savefig(output_dir / "pca_clusters.png", dpi=160)
    plt.close(fig)
    return center_id


def normalize_disulfide_bonds(disulfide_bonds) -> list[tuple[int, int]]:
    if disulfide_bonds in (None, False):
        return []
    if not isinstance(disulfide_bonds, list):
        raise ValueError("simulation_box.disulfide_bonds must be a list of two-residue index pairs.")

    normalized = []
    for pair in disulfide_bonds:
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            raise ValueError("Each simulation_box.disulfide_bonds entry must contain exactly two residue indices.")
        normalized.append((int(pair[0]), int(pair[1])))
    return normalized


def create_simulation_boxes(rmsa, np, af2sim, output_dir: Path, names: list[str], center_id, config: dict) -> None:
    box_config = config.get("simulation_box", {})
    if not box_config.get("enabled", False):
        return

    n_features = int(config.get("features", {}).get("n_features", 200))
    boxes_dir = ensure_dir(output_dir / "simulation_boxes")
    index = None

    pdb_names = [rmsa.pdb_name[i] for i in center_id]
    for pdb_name in pdb_names:
        prefix = Path(pdb_name).name
        box = af2sim.SimulationBox(pdb_name)
        create_box_kwargs = {
            "pH": float(box_config.get("pH", 7.0)),
            "padding": float(box_config.get("water_padding", 10.0)),
            "water_model": box_config.get("water_model", "tip3p"),
            "positiveIon": box_config.get("cation", "Na+"),
            "negativeIon": box_config.get("anion", "Cl-"),
            "ionicStrength": float(box_config.get("molar", 0.0)),
        }
        create_box_kwargs.update(
            callable_kwargs(
                box.create_box,
                box_config,
                exclude={
                    "enabled",
                    "pH",
                    "water_padding",
                    "water_model",
                    "cation",
                    "anion",
                    "molar",
                    "create_box_kwargs",
                    "disulfide_bonds",
                },
                extra_key="create_box_kwargs",
            )
        )
        box.create_box(**create_box_kwargs)

        disulfide_bonds = normalize_disulfide_bonds(box_config.get("disulfide_bonds"))
        if disulfide_bonds:
            box.add_disulfide_bond(disulfide_bonds)

        box.save_pdb(str(boxes_dir / f"boxed_{prefix}"))

        new_index = box.map_atom_index([tuple(rmsa.atom_pairs[name]) for name in names[:n_features]])
        if index is None:
            index = new_index
            np.save(boxes_dir / "index.npy", index)
        else:
            if not np.equal(index, new_index).all():
                raise ValueError("Feature atom index mapping is inconsistent across simulation boxes.")


class AF2RAVEPipeline:
    def __init__(self, config: dict):
        self.config = config
        self.paths = config.get("paths", {})
        output_dir = self.paths.get("output_dir")
        if not output_dir:
            raise ValueError("paths.output_dir is required.")
        self.output_dir = ensure_dir(Path(output_dir).expanduser().resolve())
        self.log_file = self.output_dir / "af2rave_pipeline.log"
        self.log_file.write_text("")
        configure_logging(self.log_file)
        logging.info("Starting AF2RAVE Pipeline")
        write_config_snapshot(self.config)

    def build_summary(
        self,
        status: str = "success",
        error: str | None = None,
        source_dir: Path | None = None,
        source_dirs: list[Path] | None = None,
        structures: list[Path] | None = None,
        reference_structure: Path | None = None,
        reference_score: float | None = None,
        filtered_count: int | None = None,
        cluster_center_count: int | None = None,
        feature_cache_written: Path | None = None,
        feature_cache_reused: Path | None = None,
    ) -> dict:
        structures = structures or []
        source_dirs = source_dirs or ([] if source_dir is None else [source_dir])
        clustering = self.config.get("clustering", {})
        clustering_method = str(clustering.get("method", "regular_space")).lower().replace("-", "_")
        output_candidates = {
            "output_dir": self.output_dir,
            "log_file": self.log_file,
        }
        if status == "success":
            output_candidates.update(
                {
                    "report": self.output_dir / "af2rave_report.txt",
                    "structures_dir": self.output_dir / "structures",
                    "structure_collection_manifest_csv": self.output_dir / "structure_collection_manifest.csv",
                    "plddt_scores_csv": self.output_dir / "plddt_scores.csv",
                    "feature_stats_csv": self.output_dir / "feature_stats.csv",
                    "feature_stats_png": self.output_dir / "feature_stats.png",
                    "filtered_structures_txt": self.output_dir / "filtered_structures.txt",
                    "selected_features_csv": self.output_dir / "selected_features.csv",
                    "feature_cv_histogram_png": self.output_dir / "feature_cv_histogram.png",
                    "cluster_centers_txt": self.output_dir / "cluster_centers.txt",
                    "pca_clusters_png": self.output_dir / "pca_clusters.png",
                    "pca_unavailable_txt": self.output_dir / "pca_unavailable.txt",
                }
            )
            if feature_cache_written is not None:
                output_candidates["ranked_feature_values_npz"] = feature_cache_written
            if clustering_method in {"kmeans", "k_means"}:
                output_candidates["cluster_assignments_csv"] = self.output_dir / "cluster_assignments.csv"
            if clustering_method in {"regular_space", "regular", "rsc"} and clustering.get("target_centers") is not None:
                output_candidates["regular_space_target_search_csv"] = (
                    self.output_dir / "regular_space_target_search.csv"
                )
            if clustering.get("copy_cluster_centers", clustering.get("copy_centers", False)):
                output_candidates["cluster_center_structures_dir"] = cluster_centers_path(self.output_dir, clustering)
            if self.config.get("simulation_box", {}).get("enabled", False):
                output_candidates["simulation_boxes_dir"] = self.output_dir / "simulation_boxes"

        outputs = {name: path for name, path in output_candidates.items() if path.exists()}
        warnings = []
        if status == "success" and not structures:
            warnings.append("No structures were collected.")

        if clustering_method in {"regular_space", "regular", "rsc"}:
            clustering_reproducibility = {
                "method": "regular_space",
                "random_seed": None,
                "seed_status": "not_applicable",
                "order_dependency": "deterministic greedy traversal of the recorded structure ordering",
            }
        else:
            clustering_reproducibility = {
                "method": clustering_method,
                "random_seed": int(clustering.get("random_state", clustering.get("randomseed", 0))),
                "seed_status": "configured",
            }

        return {
            "stage": "af2rave",
            "status": status,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "inputs": {
                "config_paths": self.paths,
                "source_dir": source_dir,
                "source_dirs": source_dirs,
                "reused_feature_cache": feature_cache_reused,
            },
            "reproducibility": {
                "runtime_versions": runtime_versions(),
                "structure_ordering": structure_ordering_metadata(structures),
                "clustering": clustering_reproducibility,
            },
            "configuration": {
                "input_structures": {
                    "structures_dir": self.paths.get("structures_dir"),
                    "structure_dirs": self.paths.get("structure_dirs", []),
                    "structure_dirs_file": self.paths.get("structure_dirs_file"),
                },
                "filtering": self.config.get("filtering", {}),
                "features": self.config.get("features", {}),
                "clustering": self.config.get("clustering", {}),
                "simulation_box": self.config.get("simulation_box", {}),
            },
            "outputs": outputs,
            "metrics": {
                "structures_collected": len(structures),
                "structures_after_filtering": filtered_count,
                "cluster_centers": cluster_center_count,
                "reference_structure": reference_structure,
                "reference_residue_averaged_plddt": reference_score,
            },
            "warnings": warnings,
            "error": error,
        }

    def run_all(self) -> dict:
        source_dirs = resolve_structure_dirs(self.paths)
        structures = collect_structures(source_dirs, self.output_dir)
        logging.info("Collected %d structures.", len(structures))
        plddt_scale = str(self.paths.get("plddt_bfactor_scale", "0_100"))
        if plddt_scale not in PLDDT_BFACTOR_SCALES:
            raise ValueError("paths.plddt_bfactor_scale must be '0_100' or '0_1'.")
        plddt_ranges = parse_plddt_residue_ranges(self.paths.get("plddt_residue_ranges", []))
        plddt_scores = calculate_plddt_scores(
            structures,
            plddt_scale,
            plddt_ranges,
            skip_invalid_structures=bool(self.paths.get("reference_structure")),
        )
        write_plddt_scores(self.output_dir, plddt_scores)
        logging.info("Wrote full pLDDT ranking: %s", self.output_dir / "plddt_scores.csv")

        np, plt, af2feature, af2sim = import_af2rave_modules()
        reference_structure = select_reference_structure(self.paths, structures, source_dirs, plddt_scores)
        rmsa_kwargs = {}
        if reference_structure is not None:
            rmsa_kwargs["ref_pdb"] = str(reference_structure)
        rmsa = af2feature.FeatureSelection([str(path) for path in structures], **rmsa_kwargs)

        write_feature_stats(
            rmsa,
            np,
            plt,
            self.output_dir,
            self.config.get("filtering", {}),
            self.config.get("analysis", {}),
        )
        apply_filters(rmsa, self.config.get("filtering", {}))
        if len(rmsa) == 0:
            raise ValueError("No structures remain after filtering; relax the active filtering cutoffs.")

        with (self.output_dir / "filtered_structures.txt").open("w") as handle:
            for pdb_name in rmsa.pdb_name:
                handle.write(f"{pdb_name}\n")
        logging.info("%d structures remain after filtering.", len(rmsa))

        names, _, feature_cache_written, feature_cache_reused = rank_features(
            rmsa, np, plt, self.output_dir, self.config.get("features", {})
        )
        alignment_reference = reference_structure or structures[0]
        center_id = cluster_and_plot(
            rmsa, np, plt, self.output_dir, names, self.config, alignment_reference
        )
        create_simulation_boxes(rmsa, np, af2sim, self.output_dir, names, center_id, self.config)

        report = self.output_dir / "af2rave_report.txt"
        report_reference = reference_structure or structures[0]
        reference_score = plddt_scores.get(report_reference)
        report.write_text(
            "\n".join(
                [
                    "AF2RAVE pipeline complete.",
                    f"Output directory: {self.output_dir}",
                    f"Structures collected: {len(structures)}",
                    f"Reference structure: {report_reference}",
                    f"Reference residue-averaged pLDDT: {reference_score if reference_score is not None else 'not available'}",
                    f"pLDDT score table: {self.output_dir / 'plddt_scores.csv'}",
                    f"Structures after filtering: {len(rmsa)}",
                    f"Cluster centers: {len(center_id)}",
                    f"Log file: {self.log_file}",
                    "",
                ]
            )
        )
        logging.info("AF2RAVE pipeline complete. Report: %s", report)
        return self.build_summary(
            source_dir=source_dirs[0],
            source_dirs=source_dirs,
            structures=structures,
            reference_structure=report_reference,
            reference_score=reference_score,
            filtered_count=len(rmsa),
            cluster_center_count=len(center_id),
            feature_cache_written=feature_cache_written,
            feature_cache_reused=feature_cache_reused,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run AF2RAVE analysis on one or more predicted-structure directories.")
    parser.add_argument("--config", required=True, help="Path to AF2RAVE YAML config.")
    parser.add_argument("--json-out", help="Optional path for a machine-readable JSON run summary")
    args = parser.parse_args()

    config = load_config(args.config)
    json_out = args.json_out or config.get("reporting", {}).get("json_out")
    pipeline = None
    try:
        pipeline = AF2RAVEPipeline(config)
        summary = pipeline.run_all()
        write_json_summary(summary, json_out or pipeline.output_dir / "dockflow_af2rave_summary.json")
    except Exception as exc:
        if pipeline is not None:
            failure = pipeline.build_summary(status="failed", error=str(exc))
            summary_path = json_out or pipeline.output_dir / "dockflow_af2rave_summary.json"
        else:
            failure = {
                "stage": "af2rave",
                "status": "failed",
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "inputs": {"config_path": str(Path(args.config).expanduser().resolve())},
                "reproducibility": {"runtime_versions": runtime_versions()},
                "configuration": config,
                "error": str(exc),
            }
            summary_path = json_out or Path(args.config).expanduser().resolve().with_name(
                "dockflow_af2rave_failure.json"
            )
        write_json_summary(failure, summary_path)
        raise


if __name__ == "__main__":
    main()
