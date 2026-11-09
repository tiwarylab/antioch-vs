#!/usr/bin/env python3
"""
AF2RAVE simulation-box and unbiased-MD pipeline.

This script is intentionally separate from the AF2RAVE feature-selection
pipeline. It starts from selected structure PDBs, maps cached feature atom pairs
into boxed systems, and optionally launches AF2RAVE UnbiasedSimulation runs.
"""

import argparse
import csv
import json
import logging
import os
import shutil
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path

import yaml


def load_config(config_path: str | Path) -> dict:
    with Path(config_path).open("r") as handle:
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


def write_config_snapshot(log_file: Path, config: dict) -> None:
    with log_file.open("a") as handle:
        handle.write("\n" + "=" * 60 + "\n")
        handle.write("--- CONFIG SNAPSHOT ---\n")
        handle.write("=" * 60 + "\n")
        handle.write(yaml.safe_dump(config, sort_keys=False))
        handle.write("\n")


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


def resolve_path(path_value: str | Path | None, base_dir: Path | None = None) -> Path | None:
    if path_value is None:
        return None
    path = Path(path_value).expanduser()
    if not path.is_absolute() and base_dir is not None:
        path = base_dir / path
    return path.resolve()


def natural_sort_key(path: Path) -> list[object]:
    import re

    return [int(text) if text.isdigit() else text.lower() for text in re.split(r"(\d+)", path.name)]


def collect_structures(paths: dict) -> list[Path]:
    structures_dir = resolve_path(paths.get("structures_dir"))
    structures_file = resolve_path(paths.get("structures_file"))

    if structures_file:
        if not structures_file.exists():
            raise FileNotFoundError(f"Structure list not found: {structures_file}")
        structures = [
            resolve_path(line.strip())
            for line in structures_file.read_text().splitlines()
            if line.strip() and not line.strip().startswith("#")
        ]
    elif structures_dir:
        if not structures_dir.exists():
            raise FileNotFoundError(f"Structure directory not found: {structures_dir}")
        structures = sorted(structures_dir.glob("*.pdb"), key=natural_sort_key)
    else:
        raise ValueError("Provide paths.structures_dir or paths.structures_file.")

    structures = [path for path in structures if path is not None]
    missing = [path for path in structures if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing input structure(s): {missing[:5]}")
    if not structures:
        raise FileNotFoundError("No PDB structures found for AF2RAVE MD setup.")
    return structures


def load_feature_pairs(np, feature_cache_npz: Path, n_features: int) -> tuple[list[str], list[tuple[int, int]]]:
    if not feature_cache_npz.exists():
        raise FileNotFoundError(f"Feature cache not found: {feature_cache_npz}")

    with np.load(feature_cache_npz, allow_pickle=True) as cache:
        names = [str(name) for name in cache["names"].tolist()]
        atom_pairs = cache["atom_pairs"].astype(np.int_)

    if n_features < 1:
        raise ValueError("features.n_features must be at least 1.")
    if n_features > len(names):
        raise ValueError(
            f"features.n_features ({n_features}) cannot exceed cached feature count ({len(names)})."
        )

    selected_names = names[:n_features]
    selected_pairs = [tuple(map(int, pair)) for pair in atom_pairs[:n_features]]
    return selected_names, selected_pairs


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


def simulation_box_kwargs(box_config: dict) -> dict:
    create_box_kwargs = {
        "pH": float(box_config.get("pH", 7.0)),
        "padding": float(box_config.get("water_padding", 10.0)),
        "water_model": box_config.get("water_model", "tip3p"),
        "positiveIon": box_config.get("cation", "Na+"),
        "negativeIon": box_config.get("anion", "Cl-"),
        "ionicStrength": float(box_config.get("molar", 0.0)),
    }
    extra_kwargs = box_config.get("create_box_kwargs", {})
    if extra_kwargs:
        if not isinstance(extra_kwargs, dict):
            raise ValueError("simulation_box.create_box_kwargs must be a dictionary.")
        create_box_kwargs.update({key: value for key, value in extra_kwargs.items() if value is not None})
    return create_box_kwargs


def md_kwargs(md_config: dict, run_dir: Path, prefix: str, mapped_pairs: list[tuple[int, int]]) -> dict:
    cv_file = md_config.get("cv_file", f"{prefix}_colvar.dat")
    xtc_file = md_config.get("xtc_file", f"{prefix}.xtc")
    kwargs = {
        "list_of_index": mapped_pairs,
        "temp": float(md_config.get("temp", 310.0)),
        "pressure": md_config.get("pressure", 1.0),
        "dt": float(md_config.get("dt", 0.002)),
        "cutoff": float(md_config.get("cutoff", 10.0)),
        "cv_file": str(run_dir / str(cv_file)) if cv_file is not None else None,
        "cv_freq": md_config.get("cv_freq", 50),
        "xtc_file": str(run_dir / str(xtc_file)) if xtc_file is not None else None,
        "xtc_freq": md_config.get("xtc_freq", 1000),
        "append": bool(md_config.get("append", False)),
        "progress_every": md_config.get("progress_every", 1000),
    }

    if kwargs["pressure"] is not None:
        kwargs["pressure"] = float(kwargs["pressure"])
    for key in ("cv_freq", "xtc_freq", "progress_every"):
        if kwargs[key] is not None:
            kwargs[key] = int(kwargs[key])

    extra_kwargs = md_config.get("unbiased_simulation_kwargs", {})
    if extra_kwargs:
        if not isinstance(extra_kwargs, dict):
            raise ValueError("md.unbiased_simulation_kwargs must be a dictionary.")
        kwargs.update({key: value for key, value in extra_kwargs.items() if value is not None})
    return kwargs


@contextmanager
def pushd(path: Path):
    old_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old_cwd)


def import_af2rave_modules():
    try:
        import af2rave.simulation as af2sim
        import numpy as np
    except ImportError as exc:
        raise ImportError(
            "AF2RAVE MD setup requires numpy and af2rave. Activate the AF2RAVE environment first."
        ) from exc
    return np, af2sim


class AF2RAVEMDPipeline:
    def __init__(self, config: dict):
        self.config = config
        self.paths = config.get("paths", {})
        self.output_dir = ensure_dir(resolve_path(self.paths["output_dir"]))
        self.log_file = self.output_dir / "af2rave_md_pipeline.log"
        configure_logging(self.log_file)
        write_config_snapshot(self.log_file, config)

    def run_all(self) -> dict:
        np, af2sim = import_af2rave_modules()
        structures = collect_structures(self.paths)
        feature_cache = resolve_path(self.paths.get("feature_cache_npz"))
        if feature_cache is None:
            raise ValueError("Provide paths.feature_cache_npz from an AF2RAVE feature-selection run.")
        feature_config = self.config.get("features", {})
        n_features = int(feature_config.get("n_features", 200))
        feature_names, feature_pairs = load_feature_pairs(np, feature_cache, n_features)

        box_config = self.config.get("simulation_box", {})
        md_config = self.config.get("md", {})
        run_md = bool(md_config.get("enabled", False))
        if not box_config.get("enabled", True) and run_md:
            raise ValueError("simulation_box.enabled must be true when md.enabled is true.")

        boxes_dir = ensure_dir(self.output_dir / "simulation_boxes")
        runs_dir = ensure_dir(self.output_dir / "md_runs")
        disulfide_bonds = normalize_disulfide_bonds(box_config.get("disulfide_bonds"))

        manifest_rows = []
        for rank, pdb_path in enumerate(structures, start=1):
            prefix = pdb_path.stem
            run_dir = ensure_dir(runs_dir / prefix)
            boxed_pdb = boxes_dir / f"boxed_{pdb_path.name}"
            mapped_index_npy = run_dir / "feature_atom_indices.npy"
            mapped_index_csv = run_dir / "feature_atom_indices.csv"

            logging.info("[%d/%d] Creating simulation box: %s", rank, len(structures), pdb_path)
            box = af2sim.SimulationBox(str(pdb_path))
            if box_config.get("enabled", True):
                box.create_box(**simulation_box_kwargs(box_config))
                if disulfide_bonds:
                    box.add_disulfide_bond(disulfide_bonds)
                box.save_pdb(str(boxed_pdb))
                mapped_pairs = box.map_atom_index(feature_pairs)
            else:
                shutil.copy2(pdb_path, boxed_pdb)
                mapped_pairs = feature_pairs

            mapped_pairs = [tuple(map(int, pair)) for pair in mapped_pairs]
            np.save(mapped_index_npy, np.array(mapped_pairs, dtype=np.int_))
            with mapped_index_csv.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["rank", "feature_name", "input_atom_i", "input_atom_j", "boxed_atom_i", "boxed_atom_j"])
                for feature_rank, (name, original_pair, mapped_pair) in enumerate(
                    zip(feature_names, feature_pairs, mapped_pairs),
                    start=1,
                ):
                    writer.writerow([feature_rank, name, *original_pair, *mapped_pair])

            final_pdb = None
            checkpoint = None
            if run_md:
                logging.info("[%d/%d] Running AF2RAVE unbiased MD: %s", rank, len(structures), boxed_pdb)
                sim_kwargs = md_kwargs(md_config, run_dir, prefix, mapped_pairs)
                steps = int(md_config.get("steps", 50000000))
                with pushd(run_dir):
                    sim = af2sim.UnbiasedSimulation(str(boxed_pdb), **sim_kwargs)
                    sim.run(steps=steps)
                    final_pdb_name = md_config.get("final_pdb") or f"{prefix}_out.pdb"
                    checkpoint_name = md_config.get("checkpoint") or f"{prefix}_out.chk"
                    final_pdb = run_dir / str(final_pdb_name)
                    checkpoint = run_dir / str(checkpoint_name)
                    sim.save_pdb(str(final_pdb))
                    sim.save_checkpoint(str(checkpoint))

            manifest_rows.append(
                {
                    "rank": rank,
                    "input_pdb": str(pdb_path),
                    "boxed_pdb": str(boxed_pdb),
                    "run_dir": str(run_dir),
                    "feature_atom_indices_npy": str(mapped_index_npy),
                    "feature_atom_indices_csv": str(mapped_index_csv),
                    "md_enabled": run_md,
                    "final_pdb": str(final_pdb) if final_pdb else "",
                    "checkpoint": str(checkpoint) if checkpoint else "",
                }
            )

        manifest = self.output_dir / "af2rave_md_manifest.csv"
        with manifest.open("w", newline="") as handle:
            fieldnames = list(manifest_rows[0])
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(manifest_rows)

        summary = {
            "stage": "af2rave_md",
            "status": "success",
            "timestamp": datetime.now().isoformat(),
            "inputs": {
                "structures": structures,
                "feature_cache_npz": feature_cache,
                "n_features": n_features,
            },
            "configuration": {
                "simulation_box": box_config,
                "md": md_config,
            },
            "outputs": {
                "output_dir": self.output_dir,
                "log_file": self.log_file,
                "manifest": manifest,
                "simulation_boxes_dir": boxes_dir,
                "md_runs_dir": runs_dir,
            },
            "metrics": {
                "structures_processed": len(structures),
                "features_mapped": len(feature_pairs),
                "md_runs_started": len(structures) if run_md else 0,
            },
        }
        write_json_summary(summary, self.output_dir / "af2rave_md_summary.json")
        logging.info("AF2RAVE MD pipeline complete. Manifest: %s", manifest)
        return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Create AF2RAVE simulation boxes and optionally run unbiased MD.")
    parser.add_argument("--config", required=True, help="Path to AF2RAVE MD YAML config.")
    args = parser.parse_args()

    config = load_config(args.config)
    pipeline = AF2RAVEMDPipeline(config)
    pipeline.run_all()


if __name__ == "__main__":
    main()
