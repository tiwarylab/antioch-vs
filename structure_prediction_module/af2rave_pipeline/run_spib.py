#!/usr/bin/env python3
"""
run_spib.py — Run SPIB on AF2RAVE COLVAR trajectories
======================================================
Run SPIB (State Predicted Information Bottleneck) on COLVAR files generated
by generate_colvar_plumed.py.

Each COLVAR file is treated as an independent trajectory. The af2rave
Colvar reader automatically strips the time column, so all 200 pairwise
distance features are passed to SPIB as input.

Two modes of operation:

  1) Single time lag (legacy):
       python3 run_spib.py --time-lag 3000

  2) Sweep over log-spaced time lags (recommended):
       python3 run_spib.py --sweep 50 5000 12
       # Runs SPIB at 12 log-spaced points between 50 and 5000 frames.

Each time-lag result is saved as a separate .pkl file:
    spib_results/spib_dt_<lag>.pkl

Usage:
    # Sweep mode (default, recommended):
    python3 run_spib.py --sweep 50 5000 12

    # Single time lag:
    python3 run_spib.py --time-lag 3000

    # Use tICA initialization with 4 clusters:
    python3 run_spib.py --init tica:4 --sweep 50 5000 12

    # Tune beta (KL weight) and learning rate:
    python3 run_spib.py --beta 0.1 --lr 5e-4 --sweep 100 5000 8

    # Full control over hyperparameters:
    python3 run_spib.py --beta 0.02 --lr 1e-3 --batch-size 256 \\
                        --refinements 12 --patience 5 --rc-dim 3 \\
                        --sweep 50 5000 10

    # List available COLVAR files only:
    python3 run_spib.py --list

Requires (set in environment before running):
    micromamba run -n chem_env python3 run_spib.py

Author: OpenClaw
"""

import argparse
import sys,os
from pathlib import Path

import numpy as np
import yaml,pickle

from af2rave import Colvar
import af2rave.spib as af2spib


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(".") #Path(__file__).resolve().parent
COLVAR_PATTERN = "*/COLVAR"  # One COLVAR per simulation subfolder
CONFIG_DEFAULT = Path(__file__).resolve().parent / "spib_config.yaml"


def load_config(config_path: str | None) -> dict:
    """Load YAML config file; return empty dict if not found."""
    if config_path is None:
        # Try default location next to the script
        if CONFIG_DEFAULT.exists():
            config_path = str(CONFIG_DEFAULT)
        else:
            return {}

    path = Path(config_path)
    if not path.exists():
        print(f"[WARN] Config file not found: {path}")
        return {}

    with open(path) as f:
        cfg = yaml.safe_load(f) or {}
    print(f"  Loaded config: {path}")
    return cfg


def _first_not_none(*vals):
    """Return the first non-None value (CLI > YAML > default precedence)."""
    for v in vals:
        if v is not None:
            return v
    return None


def discover_colvars(base: Path) -> list[Path]:
    """Find all COLVAR files in immediate subdirectories."""
    colvars = sorted(base.glob(COLVAR_PATTERN))
    return [c for c in colvars if c.is_file()]


def log_spaced_lags(lo: int, hi: int, n: int) -> list[int]:
    """Return n integers log-spaced between lo and hi (inclusive)."""
    if lo < 1:
        raise ValueError(f"--sweep min must be >= 1, got {lo}")
    if hi < lo:
        raise ValueError(f"--sweep max ({hi}) < min ({lo})")
    if n < 2:
        raise ValueError(f"--sweep points must be >= 2, got {n}")

    lags = np.round(
        np.exp(np.linspace(np.log(lo), np.log(hi), n))
    ).astype(int)
    # Ensure unique and sorted; deduplicate rounding collisions
    lags = sorted(set(lags))
    return list(lags)


def run_spib_for_lag(colvar_paths: list[str], time_lag: int, args,
                     output_dir: Path) -> dict:
    """Run SPIB for a single time lag and save the result.

    Returns a summary dict with key metrics.
    """
    print(f"\n  {'─' * 50}")
    print(f"  Time lag = {int(time_lag)} frames "
          f"(~{time_lag * 10 / 1000:.1f} ns)")

    result_path = output_dir / f"spib_dt_{time_lag}.pkl"
    print("\tExpected result path:",str(result_path))
    if os.path.exists(str(result_path)):
        print("\tFound SPIB output. Skipping")
        resobj=pickle.load(open(str(result_path),"rb"))
        return {
            "dt": time_lag,
            "states": resobj.n_converged_states,
            "bias": np.sum(resobj.bias),
            "file": str(result_path),
        }

    # Build the kwargs dict for the SPIB kernel from CLI args
    spib_kwargs = {
        "beta": args.beta,
        "learning_rate": args.lr,
        "batch_size": args.batch_size,
        "refinements": args.refinements,
        "patience": args.patience,
        "threshold": args.threshold,
        "RC_dim": args.rc_dim,
        "seed": args.seed,
        "neuron_num1": args.neuron_num1,
        "neuron_num2": args.neuron_num2,
        "lr_scheduler_step_size": args.lr_scheduler_step,
        "lr_scheduler_gamma": args.lr_scheduler_gamma,
    }

    spib = af2spib.SPIBProcess(
        colvar_paths,
        init=args.init,
    )

    # Pass beta and other hyperparams through run() kwargs → spib kernel
    result = spib.run(time_lag=time_lag, **spib_kwargs)

    result.to_file(str(result_path))

    states = int(result.n_converged_states)
    # apparent_bias can be a scalar, numpy scalar, or ndarray
    bias_raw = result.apparent_bias
    if hasattr(bias_raw, "mean"):
        bias_val = float(np.mean(bias_raw))
    else:
        bias_val = float(bias_raw)

    print(f"  → Saved: {result_path.name}")
    print(f"    Converged states: {states}  |  "
          f"Apparent bias: {bias_val:.4f}")

    return {
        "dt": time_lag,
        "states": states,
        "bias": bias_val,
        "file": str(result_path),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Run SPIB on AF2RAVE COLVAR trajectories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:

  # Sweep 12 log-spaced time lags from 50 to 5000 frames:
  python3 %(prog)s --sweep 50 5000 12

  # Single time lag only:
  python3 %(prog)s --time-lag 3000

  # Sweep with tICA initialization (4 clusters):
  python3 %(prog)s --init tica:4 --sweep 50 5000 15

  # Tune beta (KL weight) and learning rate:
  python3 %(prog)s --beta 0.1 --lr 5e-4 --sweep 100 5000 8

  # Full control over hyperparameters:
  python3 %(prog)s --beta 0.02 --lr 1e-3 --batch-size 256 \\
                    --refinements 12 --patience 5 --rc-dim 3 \\
                    --sweep 50 5000 10

  # Only process specific folders:
  python3 %(prog)s --folders jak2_taf2_ifdc_rank1 rmsa04_r326_m2_s025 --sweep 50 5000 10

  # Load options from YAML config (CLI overrides YAML):
  python3 %(prog)s --config spib_config.yaml
  python3 %(prog)s --config my_run.yaml --beta 0.1
        """,
    )

    parser.add_argument(
        "--config", type=str, default=None,
        help="Path to YAML config file (default: spib_config.yaml next to script). "
             "CLI args override YAML values.",
    )

    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument(
        "--time-lag", type=int, default=None,
        help="Single SPIB time lag in COLVAR frames. "
             "Each frame = 10 ps (Amber ntwx=5000 × dt=2fs).",
    )
    group.add_argument(
        "--sweep", type=int, nargs=3, metavar=("MIN", "MAX", "POINTS"),
        help="Run SPIB for log-spaced time lags between MIN and MAX "
             "with POINTS number of samples. "
             "Example: --sweep 50 5000 12",
    )

    parser.add_argument(
        "--init", type=str, default="split",
        help="Label initialization: 'split' (default) or 'tica:N' "
             "(e.g. 'tica:4' for 4 clusters)",
    )

    # --- SPIB Hyperparameters ---
    hp = parser.add_argument_group("SPIB Hyperparameters")
    hp.add_argument(
        "--beta", type=float, default=0.05,
        help="KL divergence weight in loss (L = CE + beta * KL). "
             "Higher = stronger regularization. Default: 0.05",
    )
    hp.add_argument(
        "--lr", type=float, default=1e-4,
        help="Adam learning rate. Default: 1e-4",
    )
    hp.add_argument(
        "--batch-size", type=int, default=512,
        help="Training batch size. Default: 512",
    )
    hp.add_argument(
        "--refinements", type=int, default=8,
        help="Max EM-style label refinement rounds. Default: 8",
    )
    hp.add_argument(
        "--patience", type=int, default=2,
        help="Epochs of no change before stopping refinements. Default: 2",
    )
    hp.add_argument(
        "--threshold", type=float, default=0.01,
        help="Convergence threshold for state population change. Default: 0.01",
    )
    hp.add_argument(
        "--rc-dim", type=int, default=2,
        help="Number of reaction coordinates (latent dim). Default: 2",
    )
    hp.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility. Default: 42",
    )
    hp.add_argument(
        "--neuron-num1", type=int, default=32,
        help="First hidden layer neurons (encoder). Default: 32",
    )
    hp.add_argument(
        "--neuron-num2", type=int, default=128,
        help="Second hidden layer neurons (encoder). Default: 128",
    )
    hp.add_argument(
        "--lr-scheduler-step", type=int, default=5,
        help="LR scheduler step size. Default: 5",
    )
    hp.add_argument(
        "--lr-scheduler-gamma", type=float, default=0.8,
        help="LR scheduler decay factor. Default: 0.8",
    )

    parser.add_argument(
        "--output", type=str, default="spib_results",
        help="Output directory for SPIB results (default: spib_results)",
    )
    parser.add_argument(
        "--list", action="store_true",
        help="List available COLVAR files and exit",
    )
    parser.add_argument(
        "--folders", nargs="+", default=None,
        help="Process only these folder names (default: all folders with COLVAR)",
    )

    args = parser.parse_args()

    # --- Load YAML config (CLI > YAML > defaults) ---
    cfg = load_config(args.config)
    spib_cfg = cfg.get("spib", {})

    def _val(cli, yaml_key, default=None):
        """Resolve precedence: CLI value > YAML value > default."""
        if cli is not None:
            return cli
        return _first_not_none(spib_cfg.get(yaml_key), cfg.get(yaml_key), default)

    # Resolve each option through the precedence chain
    time_lag   = _val(args.time_lag, "time_lag", None)
    sweep      = _val(args.sweep, "sweep", None)
    init       = _val(args.init, "init", "split")
    beta       = _val(args.beta, "beta", 0.05)
    lr         = _val(args.lr, "lr", 1e-4)
    batch_size = _val(args.batch_size, "batch_size", 512)
    refinements = _val(args.refinements, "refinements", 8)
    patience   = _val(args.patience, "patience", 2)
    threshold  = _val(args.threshold, "threshold", 0.01)
    rc_dim     = _val(args.rc_dim, "rc_dim", 2)
    seed       = _val(args.seed, "seed", 42)
    neuron_num1 = _val(args.neuron_num1, "neuron_num1", 32)
    neuron_num2 = _val(args.neuron_num2, "neuron_num2", 128)
    lr_scheduler_step = _val(args.lr_scheduler_step, "lr_scheduler_step", 5)
    lr_scheduler_gamma = _val(args.lr_scheduler_gamma, "lr_scheduler_gamma", 0.8)
    output     = _val(args.output, "output", "spib_results")
    do_list    = _val(None if args.list else None, "list", False)  # bool flag: treat absent as None
    folders    = _val(args.folders, "folders", None)

    # Validate mutually exclusive mode
    if time_lag is None and sweep is None:
        print("[ERROR] Must specify --time-lag, --sweep, or set one in config.")
        sys.exit(1)
    if time_lag is not None and sweep is not None:
        print("[ERROR] Cannot use both --time-lag and --sweep (or both in config).")
        sys.exit(1)

    # Build a namespace that mirrors the original `args` object
    class Args:
        pass

    resolved = Args()
    resolved.time_lag = time_lag
    resolved.sweep = sweep
    resolved.init = init
    resolved.beta = beta
    resolved.lr = lr
    resolved.batch_size = batch_size
    resolved.refinements = refinements
    resolved.patience = patience
    resolved.threshold = threshold
    resolved.rc_dim = rc_dim
    resolved.seed = seed
    resolved.neuron_num1 = neuron_num1
    resolved.neuron_num2 = neuron_num2
    resolved.lr_scheduler_step = lr_scheduler_step
    resolved.lr_scheduler_gamma = lr_scheduler_gamma
    resolved.output = output
    resolved.list = do_list
    resolved.folders = folders

    args = resolved

    # --- Discover COLVAR files ---
    all_colvars = discover_colvars(BASE_DIR)

    if not all_colvars:
        print(f"[ERROR] No COLVAR files found in {BASE_DIR}")
        print("       Run generate_colvar_plumed.py first.")
        sys.exit(1)

    # Filter by folder if requested
    if args.folders:
        colvars = [c for c in all_colvars if c.parent.name in args.folders]
        missing = set(args.folders) - {c.parent.name for c in colvars}
        if missing:
            print(f"[WARN] Folders without COLVAR: {missing}")
    else:
        colvars = all_colvars

    if args.list:
        print(f"Found {len(colvars)} COLVAR files:")
        for c in colvars:
            cv = Colvar.from_file(str(c))
            print(f"  {c.parent.name}/COLVAR  ->  "
                  f"{cv.shape[0]} features x {cv.shape[1]} frames")
        sys.exit(0)

    # --- Inspect each file ---
    total_frames = 0
    min_frames = float("inf")
    for c in colvars:
        cv = Colvar.from_file(str(c))
        total_frames += cv.shape[1]
        min_frames = min(min_frames, cv.shape[1])
        print(f"  {c.parent.name}: {cv.shape[0]} features x {cv.shape[1]} frames")

    print(f"\nTotal frames across all trajectories: {total_frames}")
    print(f"Minimum frames (shortest trajectory): {min_frames}")

    # --- Determine time lags ---
    if args.sweep:
        lo, hi, n = args.sweep
        lags = log_spaced_lags(lo, hi, n)
        print(f"\nSweep mode: {len(lags)} log-spaced time lags "
              f"from {lags[0]} to {lags[-1]}")
        print(f"  Lags: {lags}")

        # Warn if any lag exceeds shortest trajectory
        for lag in lags:
            if lag >= min_frames:
                print(f"  [WARN] dt={lag} >= shortest trajectory "
                      f"({min_frames} frames) — may cause issues")
    else:
        lags = [args.time_lag]

    # --- Prepare output directory ---
    output_dir = BASE_DIR / args.output
    output_dir.mkdir(exist_ok=True)

    colvar_paths = [str(c) for c in colvars]

    print(f"\n{'=' * 60}")
    print(f"SPIB Run Summary")
    print(f"{'=' * 60}")
    print(f"  Init         : {args.init}")
    print(f"  Trajectories : {len(colvar_paths)}")
    print(f"  Time lags    : {len(lags)}")
    print(f"  Beta         : {args.beta}")
    print(f"  LR           : {args.lr}")
    print(f"  Batch size   : {args.batch_size}")
    print(f"  RC dim       : {args.rc_dim}")
    print(f"  Refinements  : {args.refinements}")
    print(f"  Patience     : {args.patience}")
    print(f"  Seed         : {args.seed}")
    print(f"  Output dir   : {output_dir}")
    print(f"{'=' * 60}")

    # --- Run SPIB for each time lag ---
    summaries = []
    for lag in lags:
        try:
            summary = run_spib_for_lag(colvar_paths, lag, args, output_dir)
            summaries.append(summary)
        except Exception as e:
            print(f"  [FAIL] dt={lag}: {e}")
            summaries.append({
                "dt": lag,
                "states": None,
                "bias": None,
                "file": None,
                "error": str(e),
            })

    # --- Final summary table ---
    print(f"\n{'=' * 60}")
    print("All Results")
    print(f"{'=' * 60}")
    print(f"  {'dt':>6}  {'States':>6}  {'Bias':>8}  File")
    print(f"  {'────':>6}  {'──────':>6}  {'──────':>8}  ──────────")
    for s in summaries:
        if s.get("error"):
            print(f"  {int(s['dt']):>6}  {'FAIL':>6}  {'—':>8}  "
                  f"[{s['error'][:40]}]")
        else:
            fname = Path(s["file"]).name
            print(f"  {int(s['dt']):>6}  {int(s['states']):>6}  "
                  f"{float(s['bias']):>8.4f}  {fname}")

    print(f"\nDone. Results in {output_dir}/")


if __name__ == "__main__":
    main()
