#!/usr/bin/env python3
"""Run one or more LocalColabFold prediction settings from YAML."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shlex
import shutil
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


FASTA_SUFFIXES = {".fa", ".faa", ".fas", ".fasta", ".fsa"}
VALID_MODES = {"predict", "msa_only", "prediction_only"}


def contains_fasta(directory: Path) -> list[Path]:
    return sorted(
        path for path in directory.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES
    )


def stage_a3m_input(a3m_path: Path, output_dir: Path) -> tuple[Path, Path]:
    """Stage an A3M in an isolated directory before calling LocalColabFold.

    LocalColabFold can treat nearby FASTA files as fresh sequence inputs and
    contact an MSA service. The staged directory deliberately contains only
    the supplied alignment and is also used as the subprocess working directory.
    """
    workspace = output_dir / "localcolabfold_a3m_input"
    workspace.mkdir(parents=True, exist_ok=True)
    fasta_paths = contains_fasta(workspace)
    if fasta_paths:
        names = ", ".join(str(path) for path in fasta_paths)
        raise RuntimeError(
            "Refusing to run LocalColabFold from an A3M workspace containing FASTA files: "
            f"{names}. Remove those files before retrying."
        )

    staged_a3m = workspace / a3m_path.name
    for stale_a3m in workspace.glob("*.a3m"):
        if stale_a3m != staged_a3m:
            stale_a3m.unlink()
            logging.info("Removed stale A3M from controlled LocalColabFold workspace: %s", stale_a3m)
    if a3m_path != staged_a3m.resolve():
        shutil.copy2(a3m_path, staged_a3m)
    logging.info("Staged precomputed A3M in isolated LocalColabFold workspace: %s", staged_a3m)
    return staged_a3m, workspace


def query_input(config: dict, output_dir: Path) -> tuple[Path, Path | None]:
    paths = config.get("paths", {})
    if paths.get("sequence_file"):
        sequence_file = Path(paths["sequence_file"]).expanduser().resolve()
        if not sequence_file.is_file():
            raise FileNotFoundError(f"Sequence input file not found: {sequence_file}")
        if sequence_file.suffix.lower() == ".a3m":
            return stage_a3m_input(sequence_file, output_dir)
        return sequence_file, None
    sequence = paths.get("sequence")
    if not sequence:
        raise ValueError("Provide paths.sequence_file or paths.sequence.")
    fasta = output_dir / "query.fasta"
    fasta.write_text(f">{config.get('jobname', 'dockflow_query')}\n{''.join(str(sequence).split())}\n")
    return fasta, None


def localcolabfold_mode(settings: dict) -> str:
    """Return the configured execution mode, retaining the old boolean alias."""
    mode = settings.get("mode")
    legacy_msa_only = settings.get("msa_only")
    if mode is None:
        mode = "msa_only" if legacy_msa_only else "predict"
    elif legacy_msa_only is not None:
        raise ValueError("Set localcolabfold.mode instead of localcolabfold.msa_only, not both.")
    mode = str(mode).lower()
    if mode not in VALID_MODES:
        choices = ", ".join(sorted(VALID_MODES))
        raise ValueError(f"Unsupported LocalColabFold mode {mode!r}. Choose one of: {choices}.")
    return mode


def rename_fasta_inputs(directory: Path, context: str) -> list[Path]:
    """Disable FASTA files so LocalColabFold cannot treat them as new queries."""
    renamed = []
    for fasta_path in contains_fasta(directory):
        candidate = fasta_path.with_name(f"{fasta_path.name}.disabled")
        suffix = 1
        while candidate.exists():
            candidate = fasta_path.with_name(f"{fasta_path.name}.disabled.{suffix}")
            suffix += 1
        fasta_path.rename(candidate)
        logging.warning(
            "%s: renamed FASTA input %s to %s so LocalColabFold uses the A3M.",
            context,
            fasta_path,
            candidate,
        )
        renamed.append(candidate)
    return renamed


def prediction_only_input(config: dict, output_dir: Path) -> tuple[Path, Path]:
    """Stage a user-supplied A3M for prediction without MSA retrieval."""
    paths = config.get("paths", {})
    if paths.get("sequence"):
        raise ValueError(
            "localcolabfold.mode: prediction_only accepts paths.sequence_file pointing to an A3M, "
            "not paths.sequence."
        )
    a3m_path_value = paths.get("sequence_file")
    if not a3m_path_value:
        raise ValueError(
            "localcolabfold.mode: prediction_only requires paths.sequence_file pointing to a deep .a3m file."
        )
    a3m_path = Path(a3m_path_value).expanduser().resolve()
    if a3m_path.suffix.lower() != ".a3m":
        raise ValueError(
            f"prediction_only requires an .a3m input file, not: {a3m_path}"
        )
    if not a3m_path.is_file():
        raise FileNotFoundError(f"A3M input file not found: {a3m_path}")
    return stage_a3m_input(a3m_path, output_dir)


def validate_msa_only_output(output_dir: Path, overwrite: bool) -> None:
    """Avoid silently replacing or reusing a previously generated deep MSA."""
    existing_a3ms = sorted(path for path in output_dir.iterdir() if path.is_file() and path.suffix.lower() == ".a3m")
    if existing_a3ms and not overwrite:
        names = ", ".join(str(path) for path in existing_a3ms)
        raise FileExistsError(
            "MSA-only output directory already contains A3M file(s): "
            f"{names}. Choose a new paths.output_dir or set execution.overwrite: true."
        )


def generated_a3m(output_dir: Path) -> Path:
    """Return the single deep A3M written by an MSA-only ColabFold phase."""
    a3m_paths = sorted(path for path in output_dir.iterdir() if path.is_file() and path.suffix.lower() == ".a3m")
    if len(a3m_paths) != 1:
        found = ", ".join(str(path) for path in a3m_paths) or "none"
        raise RuntimeError(
            "The MSA-only phase must produce exactly one deep A3M for the prediction handoff; "
            f"found: {found}"
        )
    return a3m_paths[0]


def write_a3m_provenance(run_dir: Path, a3m_path: Path) -> dict:
    """Record that a prediction used a supplied alignment without MSA retrieval."""
    digest = hashlib.sha256(a3m_path.read_bytes()).hexdigest()
    provenance = {
        "input_mode": "precomputed_a3m",
        "msa_server_queried": False,
        "a3m_input": str(a3m_path),
        "a3m_sha256": digest,
    }
    output_path = run_dir / "dockflow_a3m_provenance.json"
    output_path.write_text(json.dumps(provenance, indent=2) + "\n")
    return provenance


def build_localcolabfold_command(
    command_prefix: list[str], input_file: Path, run_dir: Path, settings: dict, mode: str, overwrite: bool
) -> list[str]:
    """Build a LocalColabFold command for one configured run."""
    command = command_prefix + [
        str(input_file),
        str(run_dir),
    ]
    if mode == "msa_only":
        command.append("--msa-only")
    else:
        command.extend([
            "--max-msa",
            str(settings.get("max_msa", "8:16")),
            "--num-seeds",
            str(settings.get("num_seeds", 32)),
            "--num-recycle",
            str(settings.get("num_recycles", 1)),
        ])
        if settings.get("use_dropout", True):
            command.append("--use-dropout")
        if settings.get("save_all", False):
            command.append("--save-all")
        if settings.get("model_type"):
            command.extend(["--model-type", str(settings["model_type"])])
    if settings.get("msa_mode") and mode != "prediction_only":
        command.extend(["--msa-mode", str(settings["msa_mode"])])
    if overwrite:
        command.append("--overwrite-existing-results")
    command.extend(str(value) for value in settings.get("extra_args", []))
    return command


def prepare_prediction_run(output_dir: Path, run_name: str, overwrite: bool) -> Path:
    """Prepare one prediction output directory without allowing a FASTA query."""
    existing_run_dir = output_dir / "runs" / run_name
    if existing_run_dir.is_dir():
        rename_fasta_inputs(existing_run_dir, f"prediction run {run_name}")
    run_dir = prepare_run_dir(output_dir, run_name, overwrite)
    rename_fasta_inputs(run_dir, f"prediction run {run_name}")
    return run_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Run LocalColabFold prediction settings from YAML.")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    output_dir = resolve_output_dir(config)
    execution = config.get("execution", {})
    log_file = configure_logging(output_dir, "localcolabfold", execution.get("log_file"))
    write_config_snapshot(log_file, config)

    executable = config.get("environment", {}).get("localcolabfold_cmd", "colabfold_batch")
    command_prefix = shlex.split(str(executable))
    if not command_prefix:
        raise ValueError("environment.localcolabfold_cmd cannot be empty.")

    defaults = config.get("localcolabfold", {})
    overwrite = bool(execution.get("overwrite", False))
    runs = run_specs(config)
    run_settings = [(run, merged_mapping(defaults, run.get("settings"))) for run in runs]
    modes = {localcolabfold_mode(settings) for _, settings in run_settings}
    if len(modes) != 1:
        raise ValueError(
            "Use one localcolabfold.mode per YAML configuration. Run msa_only, predict, and "
            "prediction_only as separate workflows."
        )
    mode = modes.pop()
    if mode == "msa_only" and len(runs) != 1:
        raise ValueError(
            "localcolabfold.mode: msa_only generates one deep A3M directly in paths.output_dir. "
            "Use exactly one runs entry, then use prediction_only with the full reduced-MSA run list."
        )

    rows = []
    if mode in {"msa_only", "predict"}:
        validate_msa_only_output(output_dir, overwrite)
        msa_input, msa_workspace = query_input(config, output_dir)
        if msa_input.suffix.lower() == ".a3m":
            raise ValueError(
                f"localcolabfold.mode: {mode} requires paths.sequence or a FASTA paths.sequence_file. "
                "Use prediction_only to predict from an existing A3M."
            )
        msa_settings = run_settings[0][1]
        msa_command = build_localcolabfold_command(
            command_prefix, msa_input, output_dir, msa_settings, "msa_only", overwrite
        )
        run_command(command=msa_command, step_name="LocalColabFold deep MSA", log_file=log_file, cwd=msa_workspace)
        rows.append({"name": "deep_msa", "mode": "msa_only", "output_dir": str(output_dir), "settings": msa_settings})
        if mode == "msa_only":
            write_summary(output_dir, "localcolabfold", rows)
            return
        prediction_input, a3m_workspace = stage_a3m_input(generated_a3m(output_dir), output_dir)
    else:
        prediction_input, a3m_workspace = prediction_only_input(config, output_dir)

    for run, settings in run_settings:
        run_dir = prepare_prediction_run(output_dir, run["name"], overwrite)
        command = build_localcolabfold_command(
            command_prefix, prediction_input, run_dir, settings, "prediction_only", overwrite
        )
        run_command(
            command=command,
            step_name=f"LocalColabFold prediction-only {run['name']}",
            log_file=log_file,
            cwd=a3m_workspace,
        )
        provenance = write_a3m_provenance(run_dir, prediction_input)
        rows.append({
            "name": run["name"],
            "mode": "prediction_only",
            "output_dir": str(run_dir),
            "settings": settings,
            "a3m_provenance": provenance,
        })
    write_summary(output_dir, "localcolabfold", rows)


if __name__ == "__main__":
    main()
