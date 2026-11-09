#!/usr/bin/env python3
"""Run native AlphaFold 2 in an Apptainer container from YAML."""

from __future__ import annotations

import argparse
import shlex
from pathlib import Path

from prediction_common import (
    configure_logging,
    load_config,
    merged_mapping,
    prepare_run_dir,
    resolve_output_dir,
    run_command,
    run_specs,
    write_config_snapshot,
    write_summary,
)


def as_flag(name: str, value) -> str:
    return f"--{name}={value}"


def main() -> None:
    parser = argparse.ArgumentParser(description="Run native AlphaFold 2 in Apptainer from YAML.")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    paths = config.get("paths", {})
    environment = config.get("environment", {})
    af2 = config.get("alphafold2", {})
    fasta = Path(paths.get("sequence_file", "")).expanduser().resolve()
    if not fasta.is_file():
        raise FileNotFoundError("paths.sequence_file must point to an existing FASTA file.")

    output_dir = resolve_output_dir(config)
    execution = config.get("execution", {})
    log_file = configure_logging(output_dir, "alphafold2", execution.get("log_file"))
    write_config_snapshot(log_file, config)

    apptainer = str(environment.get("apptainer_executable", "apptainer"))
    image = Path(environment.get("container_image", "")).expanduser().resolve()
    data_dir = Path(environment.get("data_dir", "")).expanduser().resolve()
    if not image.is_file():
        raise FileNotFoundError("environment.container_image must point to an AlphaFold 2 SIF file.")
    if not data_dir.is_dir():
        raise NotADirectoryError("environment.data_dir must point to the AlphaFold 2 database directory.")

    input_mount = str(environment.get("input_mount", "/input"))
    output_mount = str(environment.get("output_mount", "/output"))
    data_mount = str(environment.get("data_mount", "/data"))
    runner_path = str(af2.get("runner_path", "/app/alphafold/run_alphafold.py"))
    defaults = af2.get("settings", {})
    overwrite = bool(execution.get("overwrite", False))
    rows = []

    for run in run_specs(config):
        settings = merged_mapping(defaults, run.get("settings"))
        run_dir = prepare_run_dir(output_dir, run["name"], overwrite)
        inner = [
            "python",
            runner_path,
            as_flag("fasta_paths", f"{input_mount}/{fasta.name}"),
            as_flag("output_dir", output_mount),
            as_flag("data_dir", data_mount),
        ]
        for key, value in settings.items():
            if value is not None:
                inner.append(as_flag(key, value))
        inner.extend(str(value) for value in run.get("extra_args", []))

        command = [apptainer, "exec"]
        if environment.get("use_gpu", True):
            command.append("--nv")
        command.extend(
            [
                "--bind",
                f"{fasta.parent}:{input_mount}",
                "--bind",
                f"{run_dir}:{output_mount}",
                "--bind",
                f"{data_dir}:{data_mount}",
            ]
        )
        for bind in environment.get("extra_binds", []):
            command.extend(["--bind", str(bind)])
        command.extend([str(image), "bash", "-lc", shlex.join(inner)])
        run_command(command, f"AlphaFold 2 {run['name']}", log_file)
        pdb_count = len(list(run_dir.rglob("*.pdb")))
        rows.append({"name": run["name"], "output_dir": str(run_dir), "pdb_files": pdb_count, "settings": settings})

    write_summary(output_dir, "alphafold2", rows)


if __name__ == "__main__":
    main()
