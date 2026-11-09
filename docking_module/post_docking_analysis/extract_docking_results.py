#!/usr/bin/env python3
"""
Dispatch post-docking result extraction to the Vina or DOCK6 extractor.

This wrapper keeps the user-facing command simple:

    python docking_module/post_docking_analysis/extract_docking_results.py \
      --software dock6 \
      -i /path/to/dock6_output \
      --out_dir /path/to/summary \
      --out_csv top_hits.csv

All arguments after --software are forwarded to the selected extractor.
"""

import argparse
import runpy
import sys
from pathlib import Path

EXTRACTORS = {
    "vina": "extract_vina_results.py",
    "vina-cpu": "extract_vina_results.py",
    "vina-gpu": "extract_vina_results.py",
    "dock6": "extract_dock6_results.py",
}


def _get_option(args, names, default=None):
    for i, value in enumerate(args):
        for name in names:
            prefix = f"{name}="
            if value.startswith(prefix):
                return value[len(prefix):]
            if value == name and i + 1 < len(args):
                return args[i + 1]
    return default


def _default_json_out(forwarded_args):
    out_dir = Path(_get_option(forwarded_args, ["--out_dir"], ".")).resolve()
    return out_dir / "dockflow_post_docking_summary.json"


def main():
    parser = argparse.ArgumentParser(
        description="Run the correct ANTIOCH-VS post-docking extractor.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  DOCK6:
    python docking_module/post_docking_analysis/extract_docking_results.py \\
      --software dock6 \\
      -i /path/to/dock6_output \\
      --out_dir /path/to/dock6_summary \\
      --out_csv top_ligands.csv \\
      --best_per_ligand \\
      --export_mol2

  Vina:
    python docking_module/post_docking_analysis/extract_docking_results.py \\
      --software vina \\
      -i /path/to/vina_output \\
      --out_dir /path/to/vina_summary \\
      --out_csv top_hits.csv
""",
    )
    parser.add_argument(
        "--software",
        required=True,
        choices=sorted(EXTRACTORS),
        help="Docking software used to generate the outputs.",
    )
    parser.add_argument(
        "--json_out",
        default=None,
        help=("Path for an agent-readable JSON summary. If omitted, writes "
              "dockflow_post_docking_summary.json inside --out_dir."),
    )

    args, forwarded_args = parser.parse_known_args()
    extractor = Path(__file__).resolve().parent / EXTRACTORS[args.software]

    if not forwarded_args:
        parser.error("Provide extractor arguments after --software, for example: -i OUTPUT_DIR --out_dir SUMMARY_DIR")

    json_out = Path(args.json_out).resolve() if args.json_out else _default_json_out(forwarded_args)
    forwarded_args = [*forwarded_args, "--json_out", str(json_out)]

    original_argv = sys.argv[:]
    sys.argv = [str(extractor), *forwarded_args]
    try:
        runpy.run_path(str(extractor), run_name="__main__")
    finally:
        sys.argv = original_argv


if __name__ == "__main__":
    main()
