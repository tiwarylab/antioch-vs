"""Shared helpers for external structure-prediction runners."""

from __future__ import annotations

import json
import logging
import re
import subprocess
from datetime import datetime
from pathlib import Path

import yaml


def load_config(config_path: str | Path) -> dict:
    with Path(config_path).open() as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError("Prediction configuration must be a YAML mapping.")
    return config


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def configure_logging(output_dir: Path, runner_name: str, configured_log_file: str | None = None) -> Path:
    if configured_log_file:
        log_file = Path(configured_log_file).expanduser()
        if not log_file.is_absolute():
            log_file = output_dir / log_file
        log_file = log_file.resolve()
    else:
        log_file = output_dir / f"{runner_name}.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    formatter = logging.Formatter("[%(asctime)s] %(message)s", datefmt="%H:%M:%S")
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    for handler in root.handlers[:]:
        root.removeHandler(handler)
        handler.close()

    for handler in (logging.StreamHandler(), logging.FileHandler(log_file)):
        handler.setFormatter(formatter)
        root.addHandler(handler)
    return log_file


def write_config_snapshot(log_file: Path, config: dict) -> None:
    with log_file.open("a") as handle:
        handle.write("\n" + "=" * 60 + "\n--- CONFIG SNAPSHOT ---\n" + "=" * 60 + "\n")
        handle.write(yaml.safe_dump(config, sort_keys=False))
        handle.write("\n")


def run_command(command: list[str], step_name: str, log_file: Path, cwd: Path | None = None) -> None:
    logging.info("STARTING: %s", step_name)
    with log_file.open("a") as handle:
        handle.write("\n" + "=" * 60 + f"\n--- {step_name.upper()} ---\n" + "=" * 60 + "\n")
        handle.write("COMMAND: " + " ".join(command) + "\n\n")
        if cwd is not None:
            handle.write(f"WORKING DIRECTORY: {cwd}\n\n")
        handle.flush()
        result = subprocess.run(command, stdout=handle, stderr=handle, text=True, cwd=cwd)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    logging.info("FINISHED: %s", step_name)


def resolve_output_dir(config: dict) -> Path:
    value = config.get("paths", {}).get("output_dir")
    if not value:
        raise ValueError("paths.output_dir is required.")
    return ensure_dir(Path(value).expanduser().resolve())


def run_specs(config: dict) -> list[dict]:
    runs = config.get("runs") or [{"name": "default"}]
    if not isinstance(runs, list) or not runs:
        raise ValueError("runs must be a non-empty YAML list.")
    names = set()
    normalized = []
    for index, run in enumerate(runs, start=1):
        if not isinstance(run, dict):
            raise ValueError("Each runs entry must be a YAML mapping.")
        name = str(run.get("name", f"run_{index}"))
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name):
            raise ValueError(f"Invalid run name: {name!r}")
        if name in names:
            raise ValueError(f"Duplicate run name: {name}")
        names.add(name)
        normalized.append({**run, "name": name})
    return normalized


def prepare_run_dir(output_dir: Path, run_name: str, overwrite: bool) -> Path:
    run_dir = output_dir / "runs" / run_name
    if run_dir.exists() and any(run_dir.iterdir()) and not overwrite:
        raise FileExistsError(
            f"Run directory is not empty: {run_dir}. Set execution.overwrite: true to reuse it."
        )
    return ensure_dir(run_dir)


def merged_mapping(base: dict | None, override: dict | None) -> dict:
    merged = dict(base or {})
    merged.update(override or {})
    return merged


def write_summary(output_dir: Path, runner: str, rows: list[dict]) -> None:
    summary = {
        "runner": runner,
        "status": "success",
        "timestamp": datetime.now().isoformat(),
        "runs": rows,
    }
    with (output_dir / f"{runner}_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
