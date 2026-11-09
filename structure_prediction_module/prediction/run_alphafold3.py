#!/usr/bin/env python3
"""Run AlphaFold 3 in an Apptainer container and normalize mmCIF outputs."""

from __future__ import annotations

import argparse
import shlex
import sys
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


def main() -> None:
    parser = argparse.ArgumentParser(description="Run AlphaFold 3 in Apptainer from YAML.")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    paths = config.get("paths", {})
    environment = config.get("environment", {})
    af3 = config.get("alphafold3", {})
    input_json = Path(paths.get("input_json", "")).expanduser().resolve()
    if not input_json.is_file():
        raise FileNotFoundError("paths.input_json must point to an existing AlphaFold 3 JSON input.")

    output_dir = resolve_output_dir(config)
    execution = config.get("execution", {})
    log_file = configure_logging(output_dir, "alphafold3", execution.get("log_file"))
    write_config_snapshot(log_file, config)

    apptainer = str(environment.get("apptainer_executable", "apptainer"))
    image = Path(environment.get("container_image", "")).expanduser().resolve()
    model_dir = Path(environment.get("model_dir", "")).expanduser().resolve()
    database_dir = Path(environment.get("database_dir", "")).expanduser().resolve()
    af2_database_dir = Path(environment.get("af2_database_dir", "")).expanduser().resolve()
    # A path inside the container, so it cannot be checked against this filesystem —
    # only that it was stated. It has no portable default: the value is whatever the
    # image was built with, and silently binding to another site's path would mount
    # the databases where AF3 does not look for them.
    af2_database_mount = str(environment.get("af2_database_mount", "")).strip()
    required_dirs = {"model_dir": model_dir, "database_dir": database_dir, "af2_database_dir": af2_database_dir}
    if not image.is_file():
        raise FileNotFoundError("environment.container_image must point to an AlphaFold 3 SIF file.")
    for name, directory in required_dirs.items():
        if not directory.is_dir():
            raise NotADirectoryError(f"environment.{name} must be an existing directory.")
    if not af2_database_mount:
        raise ValueError(
            "environment.af2_database_mount must be set to the AF2 database path expected "
            "inside the AlphaFold 3 container."
        )

    defaults = af3.get("settings", {})
    overwrite = bool(execution.get("overwrite", False))
    converter = Path(__file__).with_name("convert_mmcif_to_pdb.py")
    rows = []
    for run in run_specs(config):
        settings = merged_mapping(defaults, run.get("settings"))
        run_dir = prepare_run_dir(output_dir, run["name"], overwrite)
        inner = [
            "python",
            str(af3.get("runner_path", "/app/alphafold/run_alphafold.py")),
            f"--json_path=/root/af_input/{input_json.name}",
            "--model_dir=/root/models",
            "--db_dir=/root/public_databases",
            "--output_dir=/root/af_output",
        ]
        for key, value in settings.items():
            if value is not None:
                inner.append(f"--{key}={value}")
        inner.extend(str(value) for value in run.get("extra_args", []))

        command = [apptainer, "exec"]
        if environment.get("use_gpu", True):
            command.append("--nv")
        command.extend(
            [
                "--bind", f"{input_json.parent}:/root/af_input",
                "--bind", f"{run_dir}:/root/af_output",
                "--bind", f"{model_dir}:/root/models",
                "--bind", f"{database_dir}:/root/public_databases",
                "--bind", f"{af2_database_dir}:{af2_database_mount}",
            ]
        )
        for bind in environment.get("extra_binds", []):
            command.extend(["--bind", str(bind)])
        command.extend([str(image), "bash", "-lc", shlex.join(inner)])
        run_command(command, f"AlphaFold 3 {run['name']}", log_file)

        normalized_dir = run_dir / "structures"
        normalize = run.get("normalize_mmcif_to_pdb", af3.get("normalize_mmcif_to_pdb", True))
        if normalize:
            python_executable = str(af3.get("converter_python", sys.executable))
            run_command(
                [python_executable, str(converter), "--input-dir", str(run_dir), "--output-dir", str(normalized_dir)],
                f"Convert AlphaFold 3 mmCIF {run['name']}",
                log_file,
            )
        rows.append(
            {
                "name": run["name"],
                "output_dir": str(run_dir),
                "normalized_structures_dir": str(normalized_dir) if normalize else None,
                "settings": settings,
            }
        )
    write_summary(output_dir, "alphafold3", rows)


if __name__ == "__main__":
    main()
