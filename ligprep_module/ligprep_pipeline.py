#!/usr/bin/env python3
"""
ligprep_pipeline.py
===================
Converts SMILES or UMF files into docking-ready PDBQT / Mol2 formats
via the UMDFP toolkit (smiles2umd → umfgen → umfdump).

Usage:
    python ligprep_pipeline.py [config.yaml]

If no config is given, defaults to ligprep_config.yaml in the same directory.
The three case configs (config_case1.yaml .. 3) are minimal — missing keys
are filled with the defaults defined below.

Important umfdump quirks:
  - Output paths are resolved relative to the UMF file's directory.
  - With -s (split), files are named {output_prefix}_{N}.{ext} in that dir.
  - Without -s, a single combined file is written at the output path.
  - umfgen writes diagnostic text to stdout and creates "test.umf" in CWD.
"""

import os
import sys
import yaml
import shutil
import subprocess
import argparse
from datetime import datetime

# ── Default configuration (matches ligprep_config.yaml schema) ───────
config=None
DEFAULTS = {
    "umdfp": {
        "smiles2umd": "smiles2umd",
        "umfgen":     "umfgen",
        "umfdump":    "umfdump",
        "python":     "python3",
    },
    "source": {
        "type": "smiles",
        "path": None,
    },
    "smiles_conversion": {
        "cpus": 6,
        "batch_size": 4096,
        "charge_method": "gasteiger",
    },
    "umf_intermediate": {
        "dir": ".",
        "suffix": "",
    },
    "processing": {
        "target_format": "pdbqt",
        "output_dir": "prepared_ligands",
        "num_blocks": 1,
        "split_files": True,
        "scratch_dir": "./"
    },
    "logging": {
        "log_file": "ligprep_execution.log",
        "verbose": True,
    },
}


# ── Helpers ─────────────────────────────────────────────────────────

def deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge *override* into a copy of *base*."""
    merged = base.copy()
    for key, value in override.items():
        if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_config(config_path: str) -> dict:
    """Load YAML config and merge with defaults."""
    with open(config_path, "r") as f:
        user_cfg = yaml.safe_load(f)
    return deep_merge(DEFAULTS, user_cfg)


def resolve_binary(config: dict, name: str, script_dir: str) -> str:
    """
    Resolve a UMDFP binary path.
    Order: absolute path → relative to antioch-vs/UMDFP/bin/ sibling dir → PATH.
    """
    raw = config.get("umdfp", {}).get(name, name)
    if os.path.isabs(raw):
        return raw
    candidates = [
        os.path.join(script_dir, "..", "..", "UMDFP", "bin", name),
        raw,  # falls back to PATH lookup via subprocess
    ]
    for c in candidates:
        resolved = os.path.normpath(c)
        if os.path.isfile(resolved):
            return resolved
    return raw


def run_command(cmd: str, log_file: str = None, verbose: bool = True) -> tuple:
    """Run a shell command; log output; raise on failure.

    Returns (stdout_tail, stderr_tail) on success.
    NOTE: umfdump writes massive stdout; redirect to temp files to avoid SIGPIPE.
    """
    prc=config["processing"]
    my_scratch_dir=config["processing"]["scratch_dir"]
    if verbose:
        print(f"Using scratch directory: {my_scratch_dir}")
        print(f"  $ {cmd}")
    import tempfile
    tmp_out = tempfile.NamedTemporaryFile(
        suffix='.out', prefix='ligprep_', delete=False,
        dir=my_scratch_dir
    )
    tmp_err = tempfile.NamedTemporaryFile(
        suffix='.err', prefix='ligprep_', delete=False,
        dir=my_scratch_dir
    )
    try:
        result = subprocess.run(
            cmd, shell=True, check=False,
            stdout=tmp_out, stderr=tmp_err
        )
        tmp_out.close()
        tmp_err.close()
        with open(tmp_out.name, 'r') as f:
            out = f.read()
        with open(tmp_err.name, 'r') as f:
            err = f.read()
        if log_file:
            with open(log_file, "a") as f:
                f.write(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] CMD: {cmd}\n")
                # Write only last 200 lines of stdout to keep log manageable
                out_lines = out.strip().split('\n') if out.strip() else []
                if len(out_lines) > 200:
                    f.write(f"STDOUT (last 200 of {len(out_lines)} lines):\n")
                    f.write('\n'.join(out_lines[-200:]) + '\n')
                elif out_lines:
                    f.write(f"STDOUT:\n{out}\n")
                if err.strip():
                    f.write(f"STDERR:\n{err}\n")
                f.write("-" * 60 + "\n")
        if result.returncode != 0:
            error_msg = (
                f"Command failed (exit {result.returncode}):\n"
                f"STDOUT (last 50 lines):\n"
                f"{chr(10).join(out.strip().split(chr(10))[-50:]) if out.strip() else '(empty)'}\n"
                f"STDERR:\n{err}"
            )
            if log_file:
                with open(log_file, "a") as f:
                    f.write(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] ERROR: {error_msg}\n")
            print(error_msg, file=sys.stderr)
            raise subprocess.CalledProcessError(result.returncode, cmd, out, err)
        return out, err
    finally:
        for tmp in (tmp_out.name, tmp_err.name):
            try:
                os.unlink(tmp)
            except OSError:
                pass


def count_smiles_lines(path: str) -> int:
    """Count non-empty, non-comment lines in a SMILES file."""
    count = 0
    with open(path, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                count += 1
    return count


def count_files_in_dir(directory: str) -> int:
    """Count regular files in a directory (non-recursive)."""
    if not os.path.isdir(directory):
        return 0
    return sum(1 for f in os.listdir(directory)
               if os.path.isfile(os.path.join(directory, f)))


def clean_dir(path: str):
    """Remove and recreate a directory."""
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)


# ── Main pipeline ───────────────────────────────────────────────────

def prepare_ligands(config_path: str):
    global config
    """Run the full ligand preparation pipeline."""

    # Load and merge config
    config = load_config(config_path)
    script_dir = os.path.dirname(os.path.abspath(__file__))

    source_cfg      = config["source"]
    smiles_cfg      = config["smiles_conversion"]
    umf_int_cfg     = config["umf_intermediate"]
    proc_cfg        = config["processing"]
    log_cfg         = config["logging"]

    source_type   = source_cfg["type"]
    source_path   = source_cfg["path"]
    target_format = proc_cfg["target_format"]
    num_blocks    = proc_cfg.get("num_blocks", 1)
    output_dir    = proc_cfg["output_dir"]
    split_files   = proc_cfg.get("split_files", True)
    log_file      = log_cfg.get("log_file")
    verbose       = log_cfg.get("verbose", True)

    # Resolve binary paths
    smiles2umd_bin = resolve_binary(config, "smiles2umd", script_dir)
    umfgen_bin     = resolve_binary(config, "umfgen", script_dir)
    umfdump_bin    = resolve_binary(config, "umfdump", script_dir)
    python_bin     = resolve_binary(config, "python", script_dir)

    # Intermediate directory for UMD/UMF files
    int_dir  = os.path.abspath(proc_cfg["scratch_dir"]) #os.path.abspath(umf_int_cfg.get("dir", "."))
    suffix   = umf_int_cfg.get("suffix", "")
    os.makedirs(int_dir, exist_ok=True)

    # Clear old log
    if log_file and os.path.exists(log_file):
        os.remove(log_file)

    if verbose:
        print("=" * 60)
        print("Ligand Preparation Pipeline")
        print(f"Config   : {config_path}")
        print(f"Source   : {source_path} ({source_type})")
        print(f"Format   : {target_format}")
        print(f"Blocks   : {num_blocks}")
        print(f"Split    : {split_files}")
        print(f"ConvDump : {proc_cfg['scratch_dir']}")
        print("=" * 60)
        print()
        print("=" * 60)
        print("Binary paths for UMDFP")
        print(f"smiles2umd   : {smiles2umd_bin}")
        print(f"umfgen       : {umfgen_bin}")
        print(f"umfdump      : {umfdump_bin}")
        print(f"python       : {python_bin}")
        print("=" * 60)
        

    # ── Step 1: SMILES → UMD → UMF (if needed) ────────────────────
    umf_path = None
    total_molecules = None

    if source_type == "smiles":
        if not source_path or not os.path.exists(source_path):
            raise FileNotFoundError(f"SMILES file not found: {source_path}")

        total_molecules = count_smiles_lines(source_path)
        if verbose:
            print(f"\nFound {total_molecules} molecules in SMILES file.")

        # Derive intermediate filenames
        base_name = os.path.splitext(os.path.basename(source_path))[0]
        umd_file  = os.path.join(os.path.abspath("./"), f"{base_name}{suffix}.umd") # UMD file needs to stay in the working directory
        umf_path  = os.path.join("/mnt", f"{base_name}{suffix}.umf")

        # Clean up old intermediates
        for ext in (".umd", ".umf", ".umfp"):
            old = os.path.join(int_dir, f"{base_name}{suffix}{ext}")
            if os.path.exists(old):
                os.remove(old)
        # Also clean up per-CPU .umd files from previous runs
        for old in list(os.listdir(int_dir)):
            if old.startswith(base_name + suffix) and old.endswith(".umd"):
                os.remove(os.path.join(int_dir, old))

        cpus       = smiles_cfg.get("cpus", 12)
        batch_size = smiles_cfg.get("batch_size", 4096)

        # SMILES → UMD  (smiles2umd is a Python script; run via python_bin)
        if verbose:
            print(f"\nStep 1/2: Converting SMILES → UMD ({cpus} CPUs, batch={batch_size}) using {source_path} ...")
        smiles2umd_cmd = (
            f"{python_bin} {smiles2umd_bin} "
            f"-c {cpus} -b {batch_size} -o {umd_file} -w {proc_cfg['scratch_dir']} {source_path}"
        )
        run_command(smiles2umd_cmd, log_file=log_file, verbose=verbose)

        # UMD → UMF
        # NOTE: umfgen writes the output to hardcoded "test.umf" + "test.umfp" in CWD
        # and writes diagnostic text to stdout.  Do NOT redirect stdout.
        if verbose:
            print("Step 2/2: Converting UMD → UMF...")
        umfgen_cmd = (
            f"cd {int_dir} && "
            f"{umfgen_bin} {umd_file} {umf_path}"
        )
        run_command(umfgen_cmd, log_file=log_file, verbose=verbose)

        # Also move the pointer file (.umfp) alongside the UMF
        umf_ptr = umf_path + "p"  # e.g. bulk.umf → bulk.umfp
        test_ptr = os.path.join(int_dir, "test.umfp")
        if os.path.exists(test_ptr):
            os.rename(test_ptr, umf_ptr)

    else:
        # Source is already a UMF file
        umf_path = source_path
        if not umf_path or not os.path.exists(umf_path):
            raise FileNotFoundError(f"UMF file not found: {umf_path}")
        if verbose:
            print(f"\nUsing existing UMF file: {umf_path}")

    if not os.path.exists(umf_path):
        raise FileNotFoundError(f"UMF file was not created: {umf_path}")

    umf_dir = os.path.dirname(os.path.abspath(umf_path))
    umf_basename = os.path.basename(umf_path)

    # If total_molecules is still None (UMF source), we can't count easily.
    if total_molecules is None:
        if verbose:
            print("UMF source: molecule count unknown; using large limit.")
        total_molecules = 10_000_000

    # ── Step 2: UMF → target format via umfdump ────────────────────
    #
    # CRITICAL: umfdump resolves all output paths relative to the UMF file's
    # directory (umf_dir), NOT the current working directory.
    #
    # Strategy: create staging directories INSIDE umf_dir, run umfdump with
    # relative paths, then move results to the user-specified output location.

    if verbose:
        print(f"\nDumping molecules ({target_format}, "
              f"{'split' if split_files else 'combined'} files)...")

    # Create a staging area inside umf_dir so umfdump can find it
    stage_root = os.path.join(umf_dir, "__ligprep_staging__")
    clean_dir(stage_root)

    base_output = os.path.abspath(output_dir)

    if num_blocks == 1:
        # ── Single block ──────────────────────────────────────────
        clean_dir(base_output)

        if split_files:
            stage_dir = os.path.join(stage_root, "single")
            os.makedirs(stage_dir, exist_ok=True)

            # The output prefix for umfdump (relative to umf_dir)
            rel_prefix = os.path.relpath(
                os.path.join(stage_dir, "mol"), umf_dir
            )

            print(os.listdir("/mnt"))
            cmd = (
                f"cd {umf_dir} && "
                f"{umfdump_bin} /mnt/{umf_basename} {rel_prefix} "
                f"-f {target_format} -s -n"
            )
        else:
            # Single combined file — umfdump writes to the path directly
            stage_file = os.path.join(stage_root, f"combined.{target_format}")
            rel_file = os.path.relpath(stage_file, umf_dir)

            cmd = (
                f"cd {umf_dir} && "
                f"{umfdump_bin} ./{umf_basename} {rel_file} "
                f"-f {target_format}"
            )

        if verbose:
            print(f"\nDumping all molecules...")
        run_command(cmd, log_file=log_file, verbose=verbose)

        # Move results from staging to final location
        if split_files:
            for f in os.listdir(stage_dir):
                src = os.path.join(stage_dir, f)
                dst = os.path.join(base_output, f)
                shutil.move(src, dst)
        else:
            dst = os.path.join(base_output, f"combined.{target_format}")
            if os.path.exists(stage_file):
                shutil.move(stage_file, dst)

    else:
        # ── Multi-block ───────────────────────────────────────────
        items_per_block = (total_molecules + num_blocks - 1) // num_blocks

        if verbose:
            print(f"\nDividing {total_molecules} molecules into "
                  f"{num_blocks} blocks (~{items_per_block} each)...")

        for i in range(num_blocks):
            current_skip  = i * items_per_block
            current_limit = items_per_block

            block_dir = os.path.abspath(f"{base_output}_block_{i + 1}")
            clean_dir(block_dir)

            if split_files:
                stage_block = os.path.join(stage_root, f"block{i}")
                os.makedirs(stage_block, exist_ok=True)
                rel_prefix = os.path.relpath(
                    os.path.join(stage_block, "mol"), umf_dir
                )

                cmd = (
                    f"cd {umf_dir} && "
                    f"{umfdump_bin} ./{umf_basename} {rel_prefix} "
                    f"-f {target_format} -s "
                    f"-k {current_skip} -l {current_limit}"
                )
            else:
                stage_file = os.path.join(
                    stage_root, f"block{i}.{target_format}"
                )
                rel_file = os.path.relpath(stage_file, umf_dir)

                cmd = (
                    f"cd {umf_dir} && "
                    f"{umfdump_bin} ./{umf_basename} {rel_file} "
                    f"-f {target_format} "
                    f"-k {current_skip} -l {current_limit}"
                )

            if verbose:
                print(f"  Block {i + 1}/{num_blocks}: "
                      f"skip={current_skip}, limit={current_limit}")
            run_command(cmd, log_file=log_file, verbose=verbose)

            # Move results from staging to final location
            if split_files:
                for f in os.listdir(stage_block):
                    src = os.path.join(stage_block, f)
                    dst = os.path.join(block_dir, f)
                    os.rename(src, dst)
                os.rmdir(stage_block)
            else:
                dst = os.path.join(block_dir, f"combined.{target_format}")
                if os.path.exists(stage_file):
                    os.rename(stage_file, dst)

    # Clean up staging area
    try:
        shutil.rmtree(stage_root)
    except OSError:
        pass

    # ── Summary ────────────────────────────────────────────────────
    if verbose:
        print(f"\n✓ Ligand preparation complete.")
        print(f"  Output : {base_output}")
        if log_file and os.path.exists(log_file):
            print(f"  Log    : {log_file}")

        # Count output files
        if num_blocks == 1:
            n = count_files_in_dir(base_output)
            print(f"  Files  : {n} in {base_output}")
        else:
            for i in range(num_blocks):
                bd = f"{base_output}_block_{i + 1}"
                n = count_files_in_dir(bd)
                print(f"  Block {i+1}: {n} files in {bd}")


# ── Entry point ─────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Ligand preparation pipeline — SMILES/UMF → PDBQT/Mol2",
        epilog="Example: python ligprep_pipeline.py config_case1.yaml",
    )
    parser.add_argument(
        "config",
        nargs="?",
        default=None,
        help="Path to YAML config (default: ligprep_config.yaml)",
    )
    args = parser.parse_args()

    if args.config:
        config_path = args.config
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(script_dir, "ligprep_config.yaml")
        if not os.path.exists(config_path):
            print(f"No config specified and default not found at {config_path}",
                  file=sys.stderr)
            sys.exit(1)

    if not os.path.exists(config_path):
        print(f"Config file not found: {config_path}", file=sys.stderr)
        sys.exit(1)

    prepare_ligands(config_path)


if __name__ == "__main__":
    main()
