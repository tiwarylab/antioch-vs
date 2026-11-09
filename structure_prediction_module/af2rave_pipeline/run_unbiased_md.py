#!/usr/bin/env python3
"""
run_unbiased_md.py — build AMBER unbiased-MD systems for AF2RAVE cluster centers.

Wraps ctmd_module/md_submodule/setup_MD_from_PDB.sh (the shared soluble/membrane system
builder) in a YAML-driven, logged, manifested runner, mirroring the docking module's
one-config/one-runner pattern. Build-only: it builds system.prmtop/inpcrd per cluster
center and emits a launcher to run the MD; it does not run MD itself.

A single membrane flag selects the engine path (soluble octahedral TIP3P vs POPC bilayer),
exactly as setup_MD_from_PDB.sh does. Target-specific prep (pH, protected residues, required
disulfides, truncation) is passed through as DOCKFLOW_* environment variables.

Usage:
    python3 run_unbiased_md.py --config unbiased_md_config.yaml
"""

import argparse
import json
import logging
import os
import subprocess
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
    structures_dir = resolve_path(paths.get("cluster_centers_dir"))
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
            raise FileNotFoundError(f"Cluster-center directory not found: {structures_dir}")
        structures = sorted(structures_dir.glob("*.pdb"), key=natural_sort_key)
    else:
        raise ValueError("Provide paths.cluster_centers_dir or paths.structures_file.")

    structures = [p for p in structures if p is not None]
    missing = [p for p in structures if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing input structure(s): {missing[:5]}")
    if not structures:
        raise FileNotFoundError("No PDB cluster centers found for unbiased-MD setup.")
    return structures


def _as_list(value) -> list[str]:
    if value in (None, False, ""):
        return []
    if isinstance(value, (list, tuple)):
        return [str(v) for v in value]
    return [str(value)]


def build_prep_env(config: dict) -> dict:
    """Map the YAML prep/membrane/environment sections onto setup_MD_from_PDB.sh's DOCKFLOW_* vars."""
    prep = config.get("prep", {}) or {}
    membrane = bool((config.get("membrane", {}) or {}).get("enabled", False))
    env_cfg = config.get("environment", {}) or {}

    env = {
        "DOCKFLOW_MEMBRANE": "1" if membrane else "0",
        "DOCKFLOW_PH": str(prep.get("ph", 7.0)),
        "DOCKFLOW_PKA_CUTOFF": str(prep.get("pka_cutoff", 2.0)),
        "DOCKFLOW_SALT_CONC": str(prep.get("salt_conc", 0.15)),
        "DOCKFLOW_DISULFIDE_SG_CUTOFF": str(prep.get("disulfide_sg_cutoff", 2.3)),
        "DOCKFLOW_PROTECT_RES": " ".join(_as_list(prep.get("protect_res"))),
        "DOCKFLOW_REQUIRED_DISULFIDES": " ".join(_as_list(prep.get("required_disulfides"))),
        "DOCKFLOW_TRUNCATE_RANGES": " ".join(_as_list(prep.get("truncate_ranges"))),
    }
    if prep.get("keep_sugars"):
        env["DOCKFLOW_KEEP_SUGARS"] = "1"
    if prep.get("truncate_chain"):
        env["DOCKFLOW_TRUNCATE_CHAIN"] = str(prep["truncate_chain"])
    dockflow_root = env_cfg.get("dockflow_root")
    if dockflow_root:
        env["DOCKFLOW_ROOT"] = str(resolve_path(dockflow_root))
    if env_cfg.get("env_setup"):
        env["DOCKFLOW_ENV_SETUP"] = str(resolve_path(env_cfg["env_setup"]))
    return env


def setup_script_path(config: dict) -> Path:
    """Locate setup_MD_from_PDB.sh: under configured dockflow_root, else relative to this file."""
    env_cfg = config.get("environment", {}) or {}
    root = env_cfg.get("dockflow_root")
    if root:
        candidate = resolve_path(root) / "ctmd_module" / "md_submodule" / "setup_MD_from_PDB.sh"
        if candidate.exists():
            return candidate
    # fall back to this checkout: .../structure_prediction_module/af2rave_pipeline/ -> repo root
    repo_root = Path(__file__).resolve().parents[2]
    return repo_root / "ctmd_module" / "md_submodule" / "setup_MD_from_PDB.sh"


def emit_launcher(output_dir: Path, built_dirs: list[Path], launcher_cfg: dict) -> Path | None:
    """Write submit_md.sh that runs production for every built system (build-only: never auto-runs)."""
    if not launcher_cfg.get("emit", True):
        return None
    prod_ns = launcher_cfg.get("production_ns", 100)
    temp = launcher_cfg.get("temperature", 303.0)
    gpus = int(launcher_cfg.get("gpus_per_job", 2))
    md_cmd = launcher_cfg.get("md_run_command")

    rels = [os.path.relpath(d, output_dir) for d in built_dirs]
    lines = [
        "#!/bin/bash",
        "# submit_md.sh — run unbiased-MD production for every built system (generated by run_unbiased_md.py).",
        "# Build-only workflow: review, then submit to your scheduler. Edit MD_CMD / GPUS_PER_JOB as needed.",
        "set -uo pipefail",
        f'PROD_NS="{prod_ns}"',
        f'TEMP="{temp}"',
        f'GPUS_PER_JOB="{gpus}"',
        "# Per-system MD command: receives <system_dir> <production_ns>. Wire this to your MD image/runner",
        "# (01_Min on CPU, 02_Min2..prod on GPU; membrane barostat for the membrane case).",
    ]
    if md_cmd:
        lines.append(f'MD_CMD="{md_cmd}"')
    else:
        lines += [
            'MD_CMD="${MD_CMD:-}"',
            'if [ -z "$MD_CMD" ]; then',
            '  echo "Set MD_CMD to your per-system MD runner (receives <system_dir> <prod_ns>)."; exit 2',
            "fi",
        ]
    lines.append("")
    lines.append("SYSTEMS=(")
    lines += [f'  "{r}"' for r in rels]
    lines.append(")")
    lines += [
        "",
        'echo "Built systems: ${#SYSTEMS[@]}  (pair ${GPUS_PER_JOB}/job one-per-GPU at your scheduler)"',
        'for d in "${SYSTEMS[@]}"; do',
        '  echo "[MD] $d"',
        '  $MD_CMD "$d" "$PROD_NS"',
        "done",
    ]
    launcher = output_dir / "submit_md.sh"
    launcher.write_text("\n".join(lines) + "\n")
    launcher.chmod(0o755)
    return launcher


class UnbiasedMDPipeline:
    def __init__(self, config: dict):
        self.config = config
        self.paths = config.get("paths", {})
        self.output_dir = ensure_dir(resolve_path(self.paths["output_dir"]))
        self.log_file = self.output_dir / "unbiased_md_pipeline.log"
        configure_logging(self.log_file)
        write_config_snapshot(self.log_file, config)

    def run_all(self) -> dict:
        structures = collect_structures(self.paths)
        prep_env = build_prep_env(self.config)
        setup_script = setup_script_path(self.config)
        if not setup_script.exists():
            raise FileNotFoundError(f"setup_MD_from_PDB.sh not found at {setup_script}")
        membrane = prep_env["DOCKFLOW_MEMBRANE"] == "1"
        logging.info("Unbiased-MD build: %d centers, membrane=%s, engine=%s",
                     len(structures), membrane, setup_script)

        run_env = {**os.environ, **prep_env}
        rows = []
        built_dirs = []
        n_built = n_skipped = n_failed = 0
        for rank, pdb in enumerate(structures, start=1):
            folder = pdb.stem
            sys_dir = self.output_dir / folder
            logging.info("[%d/%d] Building system: %s", rank, len(structures), folder)
            proc = subprocess.run(
                ["bash", str(setup_script), str(pdb), folder],
                cwd=str(self.output_dir), env=run_env,
                capture_output=True, text=True,
            )
            (sys_dir).mkdir(parents=True, exist_ok=True)
            (sys_dir / "setup_MD_from_PDB.out").write_text(proc.stdout + "\n" + proc.stderr)

            prmtop = sys_dir / "system.prmtop"
            inpcrd = sys_dir / "system.inpcrd"
            if (sys_dir / "SKIPPED_INELIGIBLE").exists():
                status = "skipped_ineligible"; n_skipped += 1
            elif prmtop.exists() and inpcrd.exists():
                status = "built"; n_built += 1; built_dirs.append(sys_dir)
            else:
                status = "failed"; n_failed += 1
                logging.warning("  build failed (rc=%d) — see %s/setup_MD_from_PDB.out", proc.returncode, folder)
            rows.append({
                "rank": rank, "center": str(pdb), "system_dir": str(sys_dir),
                "status": status,
                "prmtop": str(prmtop) if prmtop.exists() else "",
                "inpcrd": str(inpcrd) if inpcrd.exists() else "",
            })

        manifest = self.output_dir / "unbiased_md_manifest.csv"
        import csv
        with manifest.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

        launcher = emit_launcher(self.output_dir, built_dirs, self.config.get("launcher", {}) or {})
        if launcher and (self.config.get("launcher", {}) or {}).get("submit"):
            logging.info("launcher.submit is true — running %s", launcher)
            subprocess.run(["bash", str(launcher)], cwd=str(self.output_dir), check=False)

        summary = {
            "stage": "unbiased_md",
            "status": "success" if n_failed == 0 else "completed_with_failures",
            "timestamp": datetime.now().isoformat(),
            "inputs": {"structures": structures, "membrane": membrane, "engine": setup_script},
            "outputs": {
                "output_dir": self.output_dir, "log_file": self.log_file,
                "manifest": manifest, "launcher": launcher,
            },
            "metrics": {
                "centers": len(structures), "built": n_built,
                "skipped_ineligible": n_skipped, "failed": n_failed,
            },
        }
        write_json_summary(summary, self.output_dir / "unbiased_md_summary.json")
        logging.info("Unbiased-MD build complete: %d built, %d skipped, %d failed. Manifest: %s",
                     n_built, n_skipped, n_failed, manifest)
        return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Build AF2RAVE unbiased-MD systems (soluble or membrane).")
    parser.add_argument("--config", required=True, help="Path to unbiased-MD YAML config.")
    args = parser.parse_args()
    config = load_config(args.config)
    UnbiasedMDPipeline(config).run_all()


if __name__ == "__main__":
    main()
