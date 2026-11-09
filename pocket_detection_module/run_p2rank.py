#!/usr/bin/env python3
"""Run P2Rank on receptor structures and export every ranked pocket."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
from pathlib import Path

from export_best_pockets import export_best_pockets


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True, help="P2Rank runner YAML configuration.")
    return parser.parse_args()


def load_config(path: Path) -> dict:
    try:
        import yaml
    except ImportError as exc:
        raise ImportError("Running P2Rank from a YAML config requires pyyaml.") from exc
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def required_path(value: object, label: str) -> Path:
    if not value:
        raise ValueError(f"Missing required configuration value: {label}")
    path = Path(str(value)).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    return path


def discover_structures(path: Path) -> list[Path]:
    extensions = (".pdb", ".cif", ".bcif")
    if path.is_file() and path.suffix.lower() in extensions:
        return [path]
    if path.is_dir():
        structures = sorted(
            candidate for candidate in path.rglob("*")
            if candidate.is_file() and candidate.suffix.lower() in extensions
        )
        if structures:
            return structures
    raise FileNotFoundError(f"No PDB, mmCIF, or BinaryCIF receptor files found at: {path}")


def resolve_executable(value: object) -> str:
    if not value:
        raise ValueError("Missing required configuration value: p2rank.executable")
    executable = str(value)
    if "/" in executable:
        path = Path(executable).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"p2rank.executable does not exist: {path}")
        if not path.stat().st_mode & 0o111:
            raise PermissionError(f"p2rank.executable is not executable: {path}")
        return str(path)
    found = shutil.which(executable)
    if not found:
        raise FileNotFoundError(f"p2rank.executable is not on PATH: {executable}")
    return found


def unique_run_name(source: Path, index: int, used: set[str]) -> str:
    name = source.name
    for suffix in (".pdb", ".cif", ".bcif"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    candidate = name
    if candidate in used:
        candidate = f"{index:03d}_{candidate}"
    used.add(candidate)
    return candidate


def find_prediction_csv(run_dir: Path) -> Path:
    matches = sorted(run_dir.rglob("*_predictions.csv"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected exactly one P2Rank '*_predictions.csv' under {run_dir}, found {len(matches)}."
        )
    return matches[0]


def value_from_row(row: dict[str, str], *names: str) -> str:
    normalized = {key.strip().lower(): value for key, value in row.items() if key is not None}
    for name in names:
        if name.lower() in normalized:
            return normalized[name.lower()]
    return ""


def export_ranked_pockets(prediction_rows: list[dict[str, str]], output_csv: Path) -> None:
    columns = [
        "receptor",
        "pocket_rank",
        "p2rank_score",
        "p2rank_probability",
        "center_x",
        "center_y",
        "center_z",
        "size_x",
        "size_y",
        "size_z",
        "source_predictions_csv",
    ]
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(prediction_rows)


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    paths = config.get("paths", {})
    p2rank = config.get("p2rank", {})
    execution = config.get("execution", {})
    docking_box = config.get("docking_box", {})

    input_structures = required_path(paths.get("input_structures"), "paths.input_structures")
    output_dir_value = paths.get("output_dir")
    if not output_dir_value:
        raise ValueError("Missing required configuration value: paths.output_dir")
    output_dir = Path(str(output_dir_value)).expanduser().resolve()
    if output_dir.exists() and any(output_dir.iterdir()) and not execution.get("overwrite", False):
        raise FileExistsError(f"Output directory is not empty: {output_dir}. Set execution.overwrite: true to reuse it.")
    if output_dir.exists() and execution.get("overwrite", False):
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    executable = resolve_executable(p2rank.get("executable"))
    threads = int(p2rank.get("threads", 1))
    if threads < 1:
        raise ValueError("p2rank.threads must be at least 1.")
    profile = p2rank.get("profile", "alphafold")
    if profile is not None and not isinstance(profile, str):
        raise ValueError("p2rank.profile must be a P2Rank config name or null.")
    extra_args = p2rank.get("extra_args", [])
    if not isinstance(extra_args, list) or not all(isinstance(item, str) for item in extra_args):
        raise ValueError("p2rank.extra_args must be a list of command-line strings.")
    size = docking_box.get("size") or [25.0, 25.0, 25.0]
    if not isinstance(size, list) or len(size) != 3:
        raise ValueError("docking_box.size must contain exactly three numbers.")
    box_size = tuple(float(value) for value in size)
    if any(value <= 0 for value in box_size):
        raise ValueError("docking_box.size values must be positive.")

    structures = discover_structures(input_structures)
    raw_dir = output_dir / "p2rank_raw"
    raw_dir.mkdir()
    log_path = output_dir / "p2rank_runner.log"
    manifest_rows = []
    ranked_rows = []
    used_names: set[str] = set()

    with log_path.open("w") as log:
        for index, source in enumerate(structures, start=1):
            receptor = unique_run_name(source, index, used_names)
            receptor_output = raw_dir / receptor
            command = [executable, "predict", "-f", str(source), "-o", str(receptor_output), "-threads", str(threads)]
            if profile:
                command.extend(["-c", profile])
            if not p2rank.get("visualizations", False):
                command.extend(["-visualizations", "0"])
            command.extend(extra_args)
            log.write("$ " + " ".join(command) + "\n")
            log.flush()
            manifest_row = {
                "receptor": receptor,
                "source_structure": str(source),
                "prediction_csv": "",
                "raw_output_dir": str(receptor_output),
                "status": "FAILED",
                "error_message": "",
            }
            try:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
                prediction_csv = find_prediction_csv(receptor_output)
                receptor_rows = []
                with prediction_csv.open(newline="") as handle:
                    reader = csv.DictReader(handle)
                    for row in reader:
                        required = {
                            "pocket rank": value_from_row(row, "rank"),
                            "center_x": value_from_row(row, "center_x"),
                            "center_y": value_from_row(row, "center_y"),
                            "center_z": value_from_row(row, "center_z"),
                        }
                        missing = [name for name, value in required.items() if value == ""]
                        if missing:
                            raise ValueError(
                                f"{prediction_csv} is missing required P2Rank values: {', '.join(missing)}"
                            )
                        receptor_rows.append(
                            {
                                "receptor": receptor,
                                "pocket_rank": required["pocket rank"],
                                "p2rank_score": value_from_row(row, "score"),
                                "p2rank_probability": value_from_row(row, "probability"),
                                "center_x": required["center_x"],
                                "center_y": required["center_y"],
                                "center_z": required["center_z"],
                                "size_x": f"{box_size[0]:g}",
                                "size_y": f"{box_size[1]:g}",
                                "size_z": f"{box_size[2]:g}",
                                "source_predictions_csv": str(prediction_csv),
                            }
                        )
                ranked_rows.extend(receptor_rows)
                manifest_row.update({"prediction_csv": str(prediction_csv), "status": "SUCCESS"})
            except (OSError, ValueError, csv.Error, subprocess.CalledProcessError) as exc:
                manifest_row["error_message"] = f"{type(exc).__name__}: {exc}"
                log.write(f"[FAILED] {receptor}: {manifest_row['error_message']}\n")
                log.flush()
                print(f"WARNING: P2Rank failed for {source.name}; continuing with remaining receptors.")
            manifest_rows.append(manifest_row)

    pockets_dir = output_dir / "pockets"
    pockets_dir.mkdir()
    pocket_csv = pockets_dir / "p2rank_ranked_pockets.csv"
    export_ranked_pockets(ranked_rows, pocket_csv)
    best_pockets_csv, n_best_pockets = export_best_pockets(pocket_csv, "p2rank", output_dir)
    manifest = output_dir / "p2rank_input_manifest.csv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["receptor", "source_structure", "prediction_csv", "raw_output_dir", "status", "error_message"],
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    (output_dir / "p2rank_run.json").write_text(
        json.dumps(
            {
                "input_structures": str(input_structures),
                "n_receptors": len(structures),
                "n_successful_receptors": sum(row["status"] == "SUCCESS" for row in manifest_rows),
                "n_failed_receptors": sum(row["status"] == "FAILED" for row in manifest_rows),
                "n_ranked_pockets": len(ranked_rows),
                "n_best_pockets": n_best_pockets,
                "best_pockets_csv": str(best_pockets_csv),
                "profile": profile,
                "threads": threads,
                "docking_box_size": box_size,
                "executable": executable,
            },
            indent=2,
        )
        + "\n"
    )
    failures = [row for row in manifest_rows if row["status"] == "FAILED"]
    if failures:
        print(f"WARNING: P2Rank failed for {len(failures)} receptor(s); see {manifest} and {log_path}.")
    if len(failures) == len(manifest_rows):
        raise RuntimeError(f"P2Rank failed for every receptor. See {manifest} and {log_path}.")
    print(f"P2Rank raw outputs: {raw_dir}")
    print(f"Ranked pockets: {pocket_csv}")
    print(f"Best docking boxes: {best_pockets_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
