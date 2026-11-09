#!/usr/bin/env python3
"""
project_colvar.py
=================
Project a COLVAR trajectory into the SPIB latent space using a trained
SPIB model (.pkl file).

Usage:
    # Project a single COLVAR:
    python3 project_colvar.py --model spib_results/spib_dt_532.pkl \
                              --colvar jak2_taf2_ifdc_rank1/COLVAR

    # Save with a custom output name:
    python3 project_colvar.py --model spib_results/spib_dt_532.pkl \
                              --colvar jak2_taf2_ifdc_rank1/COLVAR \
                              --output my_projection.dat

    # Project all COLVARs in the current directory tree:
    python3 project_colvar.py --model spib_results/spib_dt_532.pkl --all

Requires (set in environment before running):
    source /scratch/software/micromamba_root/sourceme.sh
    micromamba run -n chem_env python3 project_colvar.py ...

Author: OpenClaw
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from af2rave import Colvar
import af2rave.spib as af2spib


def discover_colvars(base: Path) -> list[Path]:
    """Find all COLVAR files in immediate subdirectories."""
    colvars = sorted(base.glob("*/COLVAR"))
    return [c for c in colvars if c.is_file()]


def project_and_save(model_path: str, colvar_path: Path,
                     output_path: Path | None = None) -> Path:
    """Load a SPIB model, project a COLVAR, and save the result.

    The output file contains: time, CV1, CV2, ..., CVk, state
    where 'state' is the predicted state number (0-indexed).
    NaN is written for points outside the convex hull of latent space.

    Returns the path to the saved file.
    """
    # Load model
    result = af2spib.SPIBResult.from_file(model_path)
    dt = result.dt
    n_states = int(result.n_converged_states)
    print(f"  Model: {Path(model_path).name}  (dt={dt}, states={n_states})")

    # Load COLVAR
    cv = Colvar.from_file(str(colvar_path))
    print(f"  Input: {colvar_path.name}  ({cv.shape[0]} features x "
          f"{cv.shape[1]} frames)")

    # Project into latent space
    projected = result.project_colvar(cv)
    n_cvs = projected.shape[0]
    n_frames = projected.shape[1]

    # Predict state labels for each frame
    # project_state_label expects shape (2, n_points)
    latent_X = projected.data  # (n_cvs, n_frames)
    state_labels = result.project_state_label(latent_X)  # (n_frames,)

    # Fix headers and add 'state' column
    cv_headers = [f"CV{i + 1}" for i in range(n_cvs)]
    all_headers = ["time"] + cv_headers + ["state"]

    # Determine output path
    if output_path is None:
        out_name = f"spib_cv_dt{dt}.dat"
        output_path = colvar_path.parent / out_name

    # Write file manually (Colvar.write doesn't support extra columns)
    with open(output_path, "w") as f:
        f.write("#! FIELDS " + " ".join(all_headers) + "\n")
        time_arr = projected.time.flatten() if projected.time is not None else np.arange(n_frames)
        for frame_idx in range(n_frames):
            t = float(time_arr[frame_idx])
            cv_vals = projected.data[:, frame_idx]
            state = float(state_labels[frame_idx])
            # Format: time CV1 CV2 ... CVk state
            parts = [f"{t:.6f}"] + [f"{float(v):.6f}" for v in cv_vals]
            # State: integer if valid, NaN if outside convex hull
            if np.isnan(state):
                parts.append("nan")
            else:
                parts.append(f"{int(state)}")
            f.write(" ".join(parts) + "\n")

    print(f"  Saved: {output_path.name}  ({n_cvs} CVs + state x {n_frames} frames)")

    return output_path


def main():
    parser = argparse.ArgumentParser(
        description="Project a COLVAR trajectory into SPIB latent space",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Project one COLVAR file:
  python3 %(prog)s --model spib_results/spib_dt_532.pkl \\
                    --colvar jak2_taf2_ifdc_rank1/COLVAR

  # Custom output filename:
  python3 %(prog)s --model spib_results/spib_dt_532.pkl \\
                    --colvar jak2_taf2_ifdc_rank1/COLVAR \\
                    --output projection.dat

  # Project all COLVARs found in subdirectories:
  python3 %(prog)s --model spib_results/spib_dt_532.pkl --all

  # Only specific folders:
  python3 %(prog)s --model spib_results/spib_dt_532.pkl \\
                    --folders jak2_taf2_ifdc_rank1 rmsa04_r326_m2_s025
        """,
    )
    parser.add_argument(
        "--model", type=str, required=True,
        help="Path to the trained SPIB model (.pkl file)",
    )
    parser.add_argument(
        "--colvar", type=str, default=None,
        help="Path to a single COLVAR file to project",
    )
    parser.add_argument(
        "--all", action="store_true", dest="project_all",
        help="Project all COLVAR files found in subdirectories",
    )
    parser.add_argument(
        "--output", type=str, default=None,
        help="Output filename (default: spib_cv_dt<N>.dat next to the COLVAR)",
    )
    parser.add_argument(
        "--folders", nargs="+", default=None,
        help="Process only these folder names (used with --all)",
    )

    args = parser.parse_args()

    # Validate model file
    model_path = Path(args.model)
    if not model_path.exists():
        print(f"[ERROR] Model file not found: {model_path}")
        sys.exit(1)

    # Determine which COLVARs to process
    base_dir = Path(__file__).resolve().parent

    if args.project_all:
        colvars = discover_colvars(base_dir)
        if not colvars:
            print("[ERROR] No COLVAR files found in subdirectories.")
            sys.exit(1)
        if args.folders:
            colvars = [c for c in colvars if c.parent.name in args.folders]
            missing = set(args.folders) - {c.parent.name for c in colvars}
            if missing:
                print(f"[WARN] Folders without COLVAR: {missing}")
        print(f"Found {len(colvars)} COLVAR files to project\n")

        outputs = []
        for cv_path in colvars:
            try:
                out = project_and_save(str(model_path), cv_path)
                outputs.append(out)
            except Exception as e:
                print(f"  [FAIL] {cv_path}: {e}")

        print(f"\nDone. Projected {len(outputs)}/{len(colvars)} trajectories.")

    elif args.colvar:
        colvar_path = Path(args.colvar)
        if not colvar_path.exists():
            print(f"[ERROR] COLVAR file not found: {colvar_path}")
            sys.exit(1)

        output_path = Path(args.output) if args.output else None
        project_and_save(str(model_path), colvar_path, output_path)

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
