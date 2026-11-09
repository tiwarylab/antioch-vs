#!/usr/bin/env python3
"""Plot DUD-E-style ROC enrichment curves for multiple docking targets.

The x-axis is the percentage of decoys recovered and the y-axis is the
percentage of actives recovered. Each legend entry includes the adjusted
DUD-E LogAUC written by calculate_enrichment.py.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_target(value: str) -> tuple[str, Path]:
    """Parse a required NAME=ENRICHMENT_OUTPUT_DIR target argument."""
    if "=" not in value:
        raise argparse.ArgumentTypeError(
            "Targets must use NAME=ENRICHMENT_OUTPUT_DIR, for example JAK2=/path/to/output."
        )
    name, path_text = value.split("=", 1)
    name = name.strip()
    path = Path(path_text).expanduser()
    if not name or not path_text.strip():
        raise argparse.ArgumentTypeError("Both target name and output directory are required.")
    return name, path


def load_target(name: str, output_dir: Path, prefix: str) -> tuple[pd.DataFrame, float]:
    curve_path = output_dir / f"{prefix}_roc_curve.csv"
    summary_path = output_dir / f"{prefix}_summary.json"
    if not curve_path.is_file():
        raise FileNotFoundError(f"{name}: ROC curve not found: {curve_path}")
    if not summary_path.is_file():
        raise FileNotFoundError(f"{name}: summary JSON not found: {summary_path}")

    curve = pd.read_csv(curve_path)
    required_columns = {"Decoy_Recall", "Active_Recall"}
    missing = required_columns.difference(curve.columns)
    if missing:
        raise ValueError(f"{name}: ROC curve is missing columns: {', '.join(sorted(missing))}")

    summary = json.loads(summary_path.read_text())
    logauc = summary.get("dude_logauc_percent")
    if logauc is None:
        raise ValueError(f"{name}: dude_logauc_percent is missing from {summary_path}")
    return curve.sort_values("Decoy_Recall"), float(logauc)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target",
        action="append",
        required=True,
        type=parse_target,
        metavar="NAME=OUTPUT_DIR",
        help="Target label and calculate_enrichment.py output directory. Repeat for each target.",
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Output PNG or PDF path.",
    )
    parser.add_argument(
        "--prefix",
        default="enrichment",
        help="Output filename prefix from calculate_enrichment.py (default: enrichment).",
    )
    parser.add_argument(
        "--title",
        default="Docking Enrichment Curves",
        help="Plot title.",
    )
    parser.add_argument("--dpi", type=int, default=300, help="Raster output resolution (default: 300).")
    args = parser.parse_args()

    names = [name for name, _ in args.target]
    if len(names) != len(set(names)):
        parser.error("Each --target name must be unique.")

    fig, ax = plt.subplots(figsize=(8, 6))
    for name, output_dir in args.target:
        curve, logauc = load_target(name, output_dir, args.prefix)
        ax.plot(
            curve["Decoy_Recall"] * 100.0,
            curve["Active_Recall"] * 100.0,
            linewidth=2,
            label=f"{name} (adjusted LogAUC = {logauc:.2f})",
        )

    ax.plot([0, 100], [0, 100], color="0.55", linestyle="--", linewidth=1, label="Random")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Decoys recovered (%)")
    ax.set_ylabel("Actives recovered (%)")
    ax.set_title(args.title)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0, frameon=False)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
