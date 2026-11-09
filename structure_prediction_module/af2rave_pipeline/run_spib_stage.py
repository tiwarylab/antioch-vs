#!/usr/bin/env python3
"""
run_spib_stage.py — orchestrate the SPIB stage for AF2RAVE unbiased-MD trajectories.

Wraps the COLVAR -> SPIB -> analysis chain (the logic in the image's `runspib`) into a
YAML-driven, logged, manifested runner, mirroring run_unbiased_md.py and the docking runners.
Beyond the stock path it handles:
  - membrane trajectories: strip lipids (POPC) as well as water when staging (stock
    extract_mdout.sh strips water only), and avoids its hardcoded absolute path;
  - truncated receptors: take the top-N features after dropping the cut ranges, and relabel
    each only_protein.pdb back to the original numbering (+ CYX->CYS) so the feature names map.

Steps (per the config): stage trajectories -> [relabel] -> select top-N features ->
generate_colvar_plumed -> run_spib (sweep) -> SPIB analysis.

Usage:
    python3 run_spib_stage.py --config spib_stage_config.yaml
"""

import argparse
import csv
import json
import logging
import os
import re
import shlex
import subprocess
from datetime import datetime
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent


def load_config(config_path):
    with Path(config_path).open("r") as handle:
        return yaml.safe_load(handle)


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def configure_logging(log_file: Path) -> None:
    formatter = logging.Formatter("[%(asctime)s] %(message)s", datefmt="%H:%M:%S")
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    for h in root.handlers[:]:
        root.removeHandler(h); h.close()
    ch = logging.StreamHandler(); ch.setFormatter(formatter); root.addHandler(ch)
    fh = logging.FileHandler(log_file); fh.setFormatter(formatter); root.addHandler(fh)


def write_config_snapshot(log_file: Path, config: dict) -> None:
    with log_file.open("a") as handle:
        handle.write("\n" + "=" * 60 + "\n--- CONFIG SNAPSHOT ---\n" + "=" * 60 + "\n")
        handle.write(yaml.safe_dump(config, sort_keys=False) + "\n")


def _json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    return value


def write_json_summary(summary: dict, json_out) -> None:
    p = Path(json_out); p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w") as handle:
        json.dump(_json_safe(summary), handle, indent=2, sort_keys=False); handle.write("\n")


def resolve_path(path_value, base_dir=None):
    if path_value is None:
        return None
    path = Path(path_value).expanduser()
    if not path.is_absolute() and base_dir is not None:
        path = base_dir / path
    return path.resolve()


def run_cmd(cmd, cwd=None, env=None, log_to=None):
    """Run a subprocess, stream rc, append combined output to log_to if given."""
    logging.info("  $ %s", " ".join(str(c) for c in cmd))
    proc = subprocess.run([str(c) for c in cmd], cwd=str(cwd) if cwd else None,
                          env=env, capture_output=True, text=True)
    if log_to is not None:
        with Path(log_to).open("a") as fh:
            fh.write(f"\n$ {' '.join(str(c) for c in cmd)}\n{proc.stdout}\n{proc.stderr}\n")
    if proc.returncode != 0:
        logging.warning("  rc=%d (see %s)", proc.returncode, log_to or "stderr")
    return proc


def parse_ranges(values):
    out = []
    for item in values or []:
        m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", str(item))
        if not m:
            raise ValueError(f"Bad range '{item}' (expected START-END).")
        out.append((int(m.group(1)), int(m.group(2))))
    return out


def feature_resids(name):
    # resid = digits after a 3-letter resname token; avoids atom-name digits (CD1, NH1, ...)
    return [int(m.group(2)) for m in re.finditer(r"([A-Z]{3})(\d+)([A-Za-z]+)", name)]


def filter_features(in_csv: Path, out_csv: Path, n_features: int, drop_ranges):
    rows = list(csv.DictReader(in_csv.open()))
    if not rows:
        raise ValueError(f"No rows in {in_csv}")
    ranges = parse_ranges(drop_ranges)
    kept = [r for r in rows
            if not any(lo <= x <= hi for x in feature_resids(r["feature_name"]) for lo, hi in ranges)]
    top = kept[:n_features]
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(top)
    logging.info("features: %d -> %d kept (dropped %s) -> top %d -> %s",
                 len(rows), len(kept), drop_ranges or "none", len(top), out_csv.name)
    return out_csv


def auto_relabel_params(relabel_cfg, drop_ranges):
    """For a single truncation range [S,E]: split = S+2 caps, offset = (E-S+1) - 2 caps."""
    split = relabel_cfg.get("split")
    offset = relabel_cfg.get("offset")
    if split is not None and offset is not None:
        return int(split), int(offset)
    ranges = parse_ranges(drop_ranges)
    if len(ranges) != 1:
        raise ValueError("relabel.split/offset must be given explicitly unless there is exactly one drop_range.")
    s, e = ranges[0]
    return s + 2, (e - s + 1) - 2


def md_folders(base: Path, trajectory_name: str):
    skip = {"spib_results", "graph_outputs", "amber_MD_template", "amber_mem_MD_template",
            "protein_outputs_direct", "protein_outputs_MD"}
    out = []
    for d in sorted(base.iterdir()):
        if not d.is_dir() or d.name in skip or d.name.startswith("."):
            continue
        if (d / trajectory_name).exists() or (d / "system.prmtop").exists():
            out.append(d)
    return out


class SPIBStagePipeline:
    def __init__(self, config: dict):
        self.config = config
        self.paths = config.get("paths", {}) or {}
        self.md_dir = resolve_path(self.paths.get("md_output_dir"))
        if self.md_dir is None or not self.md_dir.exists():
            raise FileNotFoundError(f"paths.md_output_dir not found: {self.md_dir}")
        self.work_dir = resolve_path(self.paths.get("output_dir")) or self.md_dir
        ensure_dir(self.work_dir)
        self.log_file = self.work_dir / "spib_stage.log"
        configure_logging(self.log_file)
        write_config_snapshot(self.log_file, config)
        env = config.get("environment", {}) or {}
        self.cpptraj = env.get("cpptraj", "cpptraj")
        self.colvar_python = shlex.split(env.get("colvar_python", "python3"))
        self.spib_python = shlex.split(env.get("spib_python", "python3"))
        root = env.get("dockflow_root")
        self.af2_dir = HERE if root is None else (resolve_path(root) / "structure_prediction_module" / "af2rave_pipeline")

    def stage_trajectory(self, folder: Path, membrane: bool, trajectory_name: str) -> bool:
        parm = folder / "system.prmtop"
        nc = folder / trajectory_name
        if not parm.exists() or not nc.exists():
            logging.warning("  stage skip %s (missing system.prmtop or %s)", folder.name, trajectory_name)
            return False
        mask = ":WAT,Na+,Cl-,K+" + (",PC,PA,OL" if membrane else "")
        cin = folder / ".stage.in"
        cin.write_text(
            f"parm {parm}\ntrajin {nc}\nautoimage\n"
            f"strip {mask} outprefix {folder}/stripped\n"
            f"trajout {folder}/06p5_Prod.xtc xtc\n"
            f"trajout {folder}/only_protein.pdb pdb onlyframes 1\nrun\n"
        )
        proc = run_cmd([self.cpptraj, "-i", str(cin)], cwd=folder, log_to=folder / "stage.log")
        return (folder / "06p5_Prod.xtc").exists() and (folder / "only_protein.pdb").exists() and proc.returncode == 0

    def relabel(self, folder: Path, split: int, offset: int, cyx_to_cys: bool):
        pdb = folder / "only_protein.pdb"
        cmd = self.colvar_python + [str(self.af2_dir / "relabel_only_protein.py"), str(pdb),
                                    "--split", str(split), "--offset", str(offset)]
        if not cyx_to_cys:
            cmd.append("--keep-cyx")
        run_cmd(cmd, cwd=folder, log_to=folder / "stage.log")

    def run_all(self) -> dict:
        cfg = self.config
        membrane = bool((cfg.get("membrane", {}) or {}).get("enabled", False))
        feat_cfg = cfg.get("features", {}) or {}
        relabel_cfg = cfg.get("relabel", {}) or {}
        traj_name = self.paths.get("trajectory_name", "06p5_Prod.nc")

        features_csv = resolve_path(self.paths.get("features_csv"))
        if features_csv is None or not features_csv.exists():
            raise FileNotFoundError(f"paths.features_csv not found: {features_csv}")

        folders = md_folders(self.md_dir, traj_name)
        if not folders:
            raise FileNotFoundError(f"No MD folders with {traj_name}/system.prmtop under {self.md_dir}")
        logging.info("SPIB stage: %d folders, membrane=%s, work_dir=%s", len(folders), membrane, self.work_dir)

        # 1. stage (+ optional relabel)
        staged, stage_failed = [], []
        do_relabel = bool(relabel_cfg.get("enabled", False))
        if do_relabel:
            split, offset = auto_relabel_params(relabel_cfg, feat_cfg.get("drop_ranges"))
            logging.info("relabel enabled: split=%d offset=+%d cyx_to_cys=%s",
                         split, offset, relabel_cfg.get("cyx_to_cys", True))
        for fol in folders:
            logging.info("staging %s", fol.name)
            if self.stage_trajectory(fol, membrane, traj_name):
                if do_relabel:
                    self.relabel(fol, split, offset, bool(relabel_cfg.get("cyx_to_cys", True)))
                staged.append(fol)
            else:
                stage_failed.append(fol)

        # 2. top-N features (after dropping truncated ranges)
        feat_out = self.work_dir / "features_selected.csv"
        filter_features(features_csv, feat_out, int(feat_cfg.get("n_features", 200)),
                        feat_cfg.get("drop_ranges"))

        # 3. COLVAR (needs MDAnalysis + PLUMED). generate_colvar_plumed discovers */COLVAR in cwd.
        gcp = self.af2_dir / "generate_colvar_plumed.py"
        colvar_log = self.work_dir / "colvar_generation.log"
        run_cmd(self.colvar_python + [str(gcp), "--features-csv", str(feat_out)],
                cwd=self.md_dir, env=os.environ.copy(), log_to=colvar_log)
        af2_structures = resolve_path(self.paths.get("af2_structures_dir"))
        if af2_structures and af2_structures.exists():
            for run_fol in sorted(p for p in af2_structures.iterdir() if p.is_dir()):
                run_cmd(self.colvar_python + [str(gcp), "--features-csv", str(feat_out),
                                              "--pdb-dir", str(run_fol)],
                        cwd=self.md_dir, env=os.environ.copy(), log_to=colvar_log)

        # 4. SPIB sweep
        spib_cfg = resolve_path((cfg.get("spib", {}) or {}).get("config")) or (self.af2_dir / "spib_config.yaml")
        run_cmd(self.spib_python + [str(self.af2_dir / "run_spib.py"), "--config", str(spib_cfg)],
                cwd=self.md_dir, env=os.environ.copy(), log_to=self.work_dir / "spib_run.log")

        # 5. analysis
        analysis_cfg = cfg.get("analysis", {}) or {}
        if analysis_cfg.get("enabled", True):
            an_env = os.environ.copy()
            an_env["SPIB_BASE"] = str(self.md_dir)
            an_env["SPIB_OUT"] = str(self.work_dir)
            an_env["SPIB_NUM_FINAL"] = str(analysis_cfg.get("num_final_structures", 5))
            if af2_structures and af2_structures.exists():
                an_env["AF2_STRUCTURE_ROOT"] = str(af2_structures)
            run_cmd(self.spib_python + [str(self.af2_dir / "spib_analysis.py")],
                    cwd=self.md_dir, env=an_env, log_to=self.work_dir / "spib_analysis.log")

        summary = {
            "stage": "spib_stage",
            "status": "success" if not stage_failed else "completed_with_failures",
            "timestamp": datetime.now().isoformat(),
            "inputs": {"md_output_dir": self.md_dir, "features_csv": features_csv,
                       "membrane": membrane, "n_features": int(feat_cfg.get("n_features", 200)),
                       "drop_ranges": feat_cfg.get("drop_ranges") or []},
            "outputs": {"work_dir": self.work_dir, "features_selected": feat_out,
                        "spib_results": self.work_dir / "spib_results",
                        "graph_outputs": self.work_dir / "graph_outputs", "log_file": self.log_file},
            "metrics": {"folders": len(folders), "staged": len(staged), "stage_failed": len(stage_failed)},
        }
        write_json_summary(summary, self.work_dir / "spib_stage_summary.json")
        logging.info("SPIB stage complete: %d staged, %d failed.", len(staged), len(stage_failed))
        return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the AF2RAVE SPIB stage (COLVAR -> SPIB -> analysis).")
    parser.add_argument("--config", required=True, help="Path to SPIB-stage YAML config.")
    args = parser.parse_args()
    SPIBStagePipeline(load_config(args.config)).run_all()


if __name__ == "__main__":
    main()
