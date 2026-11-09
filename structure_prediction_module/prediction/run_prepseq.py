#!/usr/bin/env python3
"""Run one or more LocalColabFold prediction settings from YAML."""

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


def query_fasta(config: dict, output_dir: Path) -> Path:
    paths = config.get("paths", {})
    if paths.get("sequence_file"):
        fasta = Path(paths["sequence_file"]).expanduser().resolve()
        if not fasta.is_file():
            raise FileNotFoundError(f"FASTA file not found: {fasta}")
        return fasta
    sequence = paths.get("sequence")
    if not sequence:
        raise ValueError("Provide paths.sequence_file or paths.sequence.")
    fasta = output_dir / "query.fasta"
    fasta.write_text(f">{config.get('jobname', 'dockflow_query')}\n{''.join(str(sequence).split())}\n")
    return fasta


def main() -> None:
    parser = argparse.ArgumentParser(description="Run LocalColabFold prediction settings from YAML.")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    output_dir = resolve_output_dir(config)
    execution = config.get("execution", {})
    log_file = configure_logging(output_dir, "localcolabfold", execution.get("log_file"))
    write_config_snapshot(log_file, config)

    executable = config.get("environment", {}).get("localcolabfold_cmd", "/usr/software/localcolabfold/.pixi/envs/default/bin/python3.12 /usr/software/localcolabfold/.pixi/envs/default/bin/colabfold_batch")
    command_prefix = shlex.split(str(executable))
    if not command_prefix:
        raise ValueError("environment.localcolabfold_cmd cannot be empty.")

    fasta = query_fasta(config, output_dir)
    defaults = config.get("localcolabfold", {})
    overwrite = bool(execution.get("overwrite", False))
    rows = []
    for run in run_specs(config):
        settings = merged_mapping(defaults, run.get("settings"))
        run_dir = prepare_run_dir(output_dir, run["name"], overwrite)
        command = command_prefix + [
            str(fasta),
            str(run_dir),
            "--num-seeds",
            str(settings.get("num_seeds", 32)),
            "--num-recycle",
            str(settings.get("num_recycles", 1)),
            "--max-msa",
            str(settings.get("max_msa", "8:16")),
        ]
        if settings.get("use_dropout", True):
            command.append("--use-dropout")
        if settings.get("save_all", False):
            command.append("--save-all")
        if settings.get("msa_mode"):
            command.extend(["--msa-mode", str(settings["msa_mode"])])
        if settings.get("model_type"):
            command.extend(["--model-type", str(settings["model_type"])])
        command.extend(str(value) for value in settings.get("extra_args", []))
        run_command(command, f"LocalColabFold {run['name']}", log_file)
        rows.append({"name": run["name"], "output_dir": str(run_dir), "settings": settings})
    write_summary(output_dir, "localcolabfold", rows)


if __name__ == "__main__":
    main()
