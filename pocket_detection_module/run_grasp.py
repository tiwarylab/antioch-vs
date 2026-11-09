#!/usr/bin/env python3
"""Run GrASP production inference on receptor PDBs and export all pockets."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from export_best_pockets import export_best_pockets


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True, help="GrASP runner YAML configuration.")
    return parser.parse_args()


def load_config(path: Path) -> dict:
    try:
        import yaml
    except ImportError as exc:
        raise ImportError("Running GrASP from a YAML config requires pyyaml.") from exc
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def required_path(value, label: str) -> Path:
    if not value:
        raise ValueError(f"Missing required configuration value: {label}")
    path = Path(value).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    return path


def discover_pdbs(path: Path) -> list[Path]:
    if path.is_file() and path.suffix.lower() == ".pdb":
        return [path]
    if path.is_dir():
        pdbs = sorted(candidate for candidate in path.rglob("*.pdb") if candidate.is_file())
        if pdbs:
            return pdbs
    raise FileNotFoundError(f"No receptor PDB files found at: {path}")


def symlink_grasp_source(grasp_root: Path, workspace: Path) -> None:
    """Expose GrASP code and pretrained models without sharing run directories."""
    excluded = {"benchmark_data_dir", "test_metrics", "__pycache__"}
    for source in grasp_root.iterdir():
        if source.name in excluded:
            continue
        (workspace / source.name).symlink_to(source)


def stage_inputs(pdbs: list[Path], input_dir: Path) -> list[dict[str, str]]:
    input_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    seen = set()
    for index, source in enumerate(pdbs, start=1):
        name = source.name
        if name in seen:
            name = f"{index:03d}_{name}"
        seen.add(name)
        staged = input_dir / name
        shutil.copy2(source, staged)
        rows.append({"source_pdb": str(source), "staged_pdb": str(staged)})
    return rows


def run_command(command: list[str], cwd: Path, log_path: Path) -> None:
    with log_path.open("a") as handle:
        handle.write("$ " + " ".join(command) + "\n")
        handle.flush()
        subprocess.run(command, cwd=cwd, check=True, stdout=handle, stderr=subprocess.STDOUT)


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    paths = config.get("paths", {})
    grasp = config.get("grasp", {})
    execution = config.get("execution", {})
    clustering = config.get("pocket_clustering", {})
    docking_box = config.get("docking_box", {})

    input_structures = required_path(paths.get("input_structures"), "paths.input_structures")
    grasp_root = required_path(grasp.get("root"), "grasp.root")
    grasp_python = required_path(grasp.get("python"), "grasp.python")
    color_python = required_path(grasp.get("color_python", grasp.get("python")), "grasp.color_python")
    output_dir_value = paths.get("output_dir")
    if not output_dir_value:
        raise ValueError("Missing required configuration value: paths.output_dir")
    output_dir = Path(output_dir_value).expanduser().resolve()
    if output_dir.exists() and any(output_dir.iterdir()) and not execution.get("overwrite", False):
        raise FileExistsError(f"Output directory is not empty: {output_dir}. Set execution.overwrite: true to reuse it.")

    required_scripts = [grasp_root / name for name in ("parse_files.py", "infer_test_set.py", "color_pdb.py")]
    missing_scripts = [str(path) for path in required_scripts if not path.is_file()]
    if missing_scripts:
        raise FileNotFoundError("GrASP root is missing required scripts:\n" + "\n".join(missing_scripts))

    output_dir.mkdir(parents=True, exist_ok=True)
    workspace_root = output_dir / "grasp_workspaces"
    raw_colors = output_dir / "grasp_colored_pdbs"
    pockets_dir = output_dir / "pockets"
    docking_boxes_dir = output_dir / "docking_boxes"
    if execution.get("overwrite", False):
        for generated_path in (workspace_root, raw_colors, pockets_dir, docking_boxes_dir):
            if generated_path.exists():
                shutil.rmtree(generated_path)
    elif workspace_root.exists():
        raise FileExistsError(f"GrASP workspace already exists: {workspace_root}")
    workspace_root.mkdir()

    staged_rows = stage_inputs(
        discover_pdbs(input_structures),
        workspace_root / "staged_inputs",
    )
    log_path = output_dir / "grasp_runner.log"
    log_path.write_text("")
    model_path = str(grasp.get("model_path", "train_full/trained_model_s_train_full_ag_multi_1680643832.8660116/cv_0/epoch_49"))
    n_tasks = int(grasp.get("n_tasks", 4))
    if n_tasks < 1:
        raise ValueError("grasp.n_tasks must be at least 1.")

    environment = dict(os.environ)
    # GrASP's pyximport defaults to ~/.pyxbld. Keep its compiled extension
    # private to this run so it cannot reuse an incompatible build cache from
    # another Python version or a failed concurrent compilation.
    environment["HOME"] = str(workspace_root)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    raw_colors.mkdir(exist_ok=True)
    colored_pdbs = []
    missing_outputs = []
    for row in staged_rows:
        staged_input = Path(row["staged_pdb"])
        receptor_workspace = workspace_root / "runs" / staged_input.stem
        receptor_workspace.mkdir(parents=True)
        symlink_grasp_source(grasp_root, receptor_workspace)
        receptor_input = receptor_workspace / "benchmark_data_dir" / "production" / "unprocessed_inputs" / staged_input.name
        receptor_input.parent.mkdir(parents=True)
        shutil.copy2(staged_input, receptor_input)
        row["staged_pdb"] = str(receptor_input)
        try:
            commands = (
                [
                    str(grasp_python),
                    "-c",
                    "import pyximport; pyximport.install(); import fast_distance_computation; "
                    "print('GrASP Cython distance extension is ready.')",
                ],
                [str(grasp_python), "parse_files.py", "production"],
                [str(grasp_python), "infer_test_set.py", model_path, "--infer_set", "production", "--n_tasks", str(n_tasks)],
                [str(color_python), "color_pdb.py", model_path],
            )
            for command in commands:
                with log_path.open("a") as handle:
                    handle.write(f"[{staged_input.name}] $ " + " ".join(command) + "\n")
                    handle.flush()
                    subprocess.run(
                        command,
                        cwd=receptor_workspace,
                        check=True,
                        stdout=handle,
                        stderr=subprocess.STDOUT,
                        env=environment,
                    )
            colored_dir = receptor_workspace / "test_metrics" / "production" / "colors" / model_path
            matches = sorted(colored_dir.rglob(f"{staged_input.stem}_probs.pdb")) if colored_dir.is_dir() else []
            if len(matches) != 1:
                raise RuntimeError(f"Expected one colored PDB for {staged_input.name}, found {len(matches)}.")
            colored_pdb = matches[0]
            shutil.copy2(colored_pdb, raw_colors / colored_pdb.name)
            colored_pdbs.append(colored_pdb)
            row.update({"colored_pdb": str(colored_pdb), "color_status": "FOUND", "failure_reason": ""})
        except (OSError, RuntimeError, subprocess.CalledProcessError) as exc:
            row.update(
                {
                    "colored_pdb": "",
                    "color_status": "MISSING",
                    "failure_reason": f"{type(exc).__name__}: {exc}",
                }
            )
            missing_outputs.append(row["source_pdb"])
            print(f"WARNING: GrASP failed for {staged_input.name}; continuing with remaining receptors.")

    manifest_path = output_dir / "grasp_input_manifest.csv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["source_pdb", "staged_pdb", "colored_pdb", "color_status", "failure_reason"],
        )
        writer.writeheader()
        writer.writerows(staged_rows)

    if not colored_pdbs:
        raise RuntimeError(
            "GrASP did not produce any colored PDBs. "
            f"See {manifest_path} and {log_path} for per-receptor failure details."
        )
    if missing_outputs:
        print(
            f"WARNING: GrASP did not produce colored PDBs for {len(missing_outputs)} receptor(s); "
            f"continuing with successful receptors. See {manifest_path} and {log_path}."
        )

    postprocessor = Path(__file__).with_name("postprocess_grasp_pockets.py")
    pocket_command = [
        sys.executable,
        str(postprocessor),
        "--scored-pdbs",
        str(raw_colors),
        "--out-dir",
        str(output_dir / "pockets"),
        "--score-threshold",
        str(clustering.get("score_threshold", 0.5)),
        "--cluster-distance",
        str(clustering.get("cluster_distance", 6.0)),
        "--min-atoms",
        str(clustering.get("min_atoms", 3)),
        "--box-padding",
        str(docking_box.get("padding", 5.0)),
        "--min-box-size",
        str(docking_box.get("min_size", 15.0)),
    ]
    fixed_box_size = docking_box.get("size")
    if fixed_box_size is not None:
        if not isinstance(fixed_box_size, list) or len(fixed_box_size) != 3:
            raise ValueError("docking_box.size must contain exactly three numbers or be null.")
        pocket_command.extend(["--fixed-box-size", *(str(float(value)) for value in fixed_box_size)])
    run_command(pocket_command, output_dir, log_path)
    best_pockets_csv, n_best_pockets = export_best_pockets(
        output_dir / "pockets" / "grasp_ranked_pockets.csv", "grasp", output_dir
    )
    (output_dir / "grasp_run.json").write_text(
        json.dumps(
            {
                "input_structures": str(input_structures),
                "n_receptors": len(staged_rows),
                "n_successful_receptors": len(staged_rows) - len(missing_outputs),
                "n_failed_receptors": len(missing_outputs),
                "n_best_pockets": n_best_pockets,
                "input_manifest": str(manifest_path),
                "ranked_pockets_csv": str(output_dir / "pockets" / "grasp_ranked_pockets.csv"),
                "best_pockets_csv": str(best_pockets_csv),
            },
            indent=2,
        )
        + "\n"
    )
    print(f"GrASP-colored PDBs: {raw_colors}")
    print(f"Ranked pockets: {output_dir / 'pockets' / 'grasp_ranked_pockets.csv'}")
    print(f"Best docking boxes: {best_pockets_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
