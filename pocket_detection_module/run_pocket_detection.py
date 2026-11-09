#!/usr/bin/env python3
"""Dispatch a unified pocket-detection YAML to GrASP or P2Rank."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True, help="Unified pocket-detection YAML configuration.")
    return parser.parse_args()


def load_method(config_path: Path) -> str:
    try:
        import yaml
    except ImportError as exc:
        raise ImportError("Running a pocket-detection YAML requires pyyaml.") from exc
    with config_path.open() as handle:
        config = yaml.safe_load(handle) or {}
    method = config.get("method")
    if method not in {"grasp", "p2rank"}:
        raise ValueError("method must be either 'grasp' or 'p2rank'.")
    return method


def main() -> int:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    if not config_path.is_file():
        raise FileNotFoundError(f"Configuration file does not exist: {config_path}")
    method = load_method(config_path)
    runner = Path(__file__).with_name(f"run_{method}.py")
    subprocess.run([sys.executable, str(runner), "--config", str(config_path)], check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
