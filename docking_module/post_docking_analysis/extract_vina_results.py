#!/usr/bin/env python3
"""
Extract Vina docking results from multiple directories.

Recursively searches for PDBQT output files, extracts the top binding affinity
(Vina score), maps it to the absolute file path, and optionally converts
the top pose (Model 1) to SDF format using Meeko's mk_export.py.

Optional batch mode:
    Input directories can be split into batches either by immediate subfolders
    or by recursively discovered leaf folders that directly contain matching
    PDBQT files. Each batch is processed independently into a batch parquet
    first, then all batches are merged into one master parquet at the end.
"""

import argparse
import json
import multiprocessing as mp
import os
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **_kwargs):
        return iterable

# Catches both decimal scores, e.g. -8.7, and integer-looking scores, e.g. -8
VINA_REGEX = re.compile(r"^REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)")
pd = None


def _json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _json_safe(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def write_json_summary(summary, json_out):
    json_path = Path(json_out)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    with json_path.open("w") as handle:
        json.dump(_json_safe(summary), handle, indent=2, sort_keys=False)
        handle.write("\n")


def build_json_summary(args, status, total_files=0, total_skipped=0, df=None, csv_path=None,
                       sdf_attempted=0, sdf_success=0, error=None):
    out_dir = Path(args.out_dir).resolve()
    parquet_path = out_dir / args.out_parquet
    csv_target = csv_path or (out_dir / args.out_csv if args.out_csv else None)

    metrics = {
        "files_scanned": int(total_files),
        "scored_rows": int(len(df)) if df is not None else 0,
        "skipped_files": int(total_skipped),
        "unique_ligand_count": int(df["Ligand"].nunique()) if df is not None and "Ligand" in df.columns else 0,
        "scores_parquet_exists": parquet_path.exists(),
        "top_hits_csv_exists": csv_target.exists() if csv_target else False,
        "sdf_export_attempted": int(sdf_attempted),
        "sdf_export_success": int(sdf_success),
    }

    if df is not None and not df.empty and "Vina_Score" in df.columns:
        metrics.update({
            "score_column": "Vina_Score",
            "best_score": float(df["Vina_Score"].min()),
            "median_score": float(df["Vina_Score"].median()),
        })

    top_hits_preview = []
    if df is not None and not df.empty:
        top_hits_preview = df.head(min(args.top_n, 10)).to_dict(orient="records")

    return {
        "stage": "post_docking_analysis",
        "extractor": "vina",
        "status": status,
        "timestamp": datetime.now().isoformat(),
        "inputs": {
            "input_dirs": [str(Path(path).resolve()) for path in args.input_dirs],
            "input_file": str(Path(args.input_file).resolve()) if args.input_file else None,
            "pattern": args.pattern,
        },
        "configuration": {
            "cpus": args.cpus,
            "batch_by_subdirs": args.batch_by_subdirs,
            "batch_discovery": args.batch_discovery,
            "export_sdf": args.export_sdf,
            "top_n": args.top_n,
            "csv_n": args.csv_n,
            "sdf_n": args.sdf_n,
        },
        "outputs": {
            "out_dir": str(out_dir),
            "scores_parquet": str(parquet_path),
            "top_hits_csv": str(csv_target) if csv_target else None,
            "sdf_dir": str(out_dir / args.sdf_dir) if args.export_sdf else None,
            "json_summary": str(Path(args.json_out).resolve()),
        },
        "metrics": metrics,
        "top_hits_preview": top_hits_preview,
        "warnings": [],
        "error": error,
    }


def _extract_data_worker(file_path):
    try:
        path_obj = Path(file_path)
        # Keep ligand names exactly as originally defined by the file stem.
        ligand_name = path_obj.stem.replace('_out', '')

        with open(file_path) as f:
            for line in f:
                match = VINA_REGEX.match(line)
                if match:
                    score = float(match.group(1))
                    return ligand_name, score, str(path_obj.resolve())
    except Exception:
        return None
    return None


def check_mk_export_available():
    """
    Check whether Meeko's mk_export.py is available on PATH before attempting
    SDF export.
    """
    if shutil.which("mk_export.py"):
        return True

    print("\nError: 'mk_export.py' was not found on your PATH.")
    print("To export SDF files, activate the docking_pipeline environment first, e.g.:")
    print("    conda activate docking_pipeline")
    print("or:")
    print("    micromamba activate docking_pipeline")
    print("If Meeko is not installed in that environment, install it, e.g.:")
    print("    pip install meeko")
    return False


def extract_and_convert_top_pose(original_pdbqt, output_dir, ligand_name):
    """
    Extracts only MODEL 1 from the Vina output and converts it to SDF.
    """
    temp_pdbqt = output_dir / f"{ligand_name}_pose1.pdbqt"

    try:
        with open(original_pdbqt) as f_in, open(temp_pdbqt, 'w') as f_out:
            for line in f_in:
                f_out.write(line)
                if line.startswith("ENDMDL"):
                    break
    except Exception as e:
        print(f"Error reading {ligand_name}: {e}")
        return False

    try:
        subprocess.run(
            ["mk_export.py", temp_pdbqt.name],
            cwd=output_dir,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        temp_pdbqt.unlink()
        return True

    except subprocess.CalledProcessError:
        print(f"Error: mk_export.py failed for {ligand_name}. Is Meeko installed correctly?")
        if temp_pdbqt.exists():
            temp_pdbqt.unlink()
        return False
    except FileNotFoundError:
        print("Error: 'mk_export.py' command not found. Activate docking_pipeline or install Meeko.")
        if temp_pdbqt.exists():
            temp_pdbqt.unlink()
        return False


def gather_files(scan_dirs, pattern):
    """
    Recursively gather files matching pattern from all scan directories.
    """
    all_files = []
    for d in scan_dirs:
        path_obj = Path(d)
        if path_obj.is_dir():
            all_files.extend([str(f) for f in path_obj.rglob(pattern)])
        else:
            print(f"Warning: Directory '{d}' does not exist or is not a folder. Skipping.")
    return all_files


def parse_vina_files(files, cpus=None, label="Parsing PDBQTs"):
    """
    Parse a list of PDBQT files and return a sorted dataframe plus counts.
    """
    if not files:
        return pd.DataFrame(columns=['Ligand', 'Vina_Score', 'File_Path']), 0, 0

    print(f"Found {len(files)} files. Parsing for Vina scores...")
    results = []
    n_procs = cpus or mp.cpu_count()

    with mp.Pool(processes=n_procs) as pool:
        iterator = pool.imap_unordered(_extract_data_worker, files)
        for res in tqdm(iterator, total=len(files), desc=label, unit="file"):
            if res:
                results.append(res)

    skipped = len(files) - len(results)
    print(f"[+] Parsed {len(results)} valid scored files out of {len(files)} total files.")
    print(f"[+] Skipped {skipped} files with no valid Vina score or unreadable content.")

    if not results:
        return pd.DataFrame(columns=['Ligand', 'Vina_Score', 'File_Path']), len(files), skipped

    df = pd.DataFrame(results, columns=['Ligand', 'Vina_Score', 'File_Path'])
    df = df.sort_values(by='Vina_Score', ascending=True).reset_index(drop=True)
    return df, len(files), skipped


def _relative_batch_name(input_dir, batch_dir):
    """
    Make a readable, unique batch name from the batch path relative to the
    original input directory.
    """
    try:
        rel = batch_dir.relative_to(input_dir)
        if str(rel) != ".":
            return str(rel).replace(os.sep, "__")
    except ValueError:
        pass
    return batch_dir.name


def get_batches_from_input_dirs(scan_dirs, pattern, discovery="immediate"):
    """
    Build batch definitions from input directories.

    discovery='immediate':
      - If an input directory has immediate subdirectories, process each
        immediate subdirectory as a separate batch.
      - If it has no immediate subdirectories, process the input directory
        itself as one batch.
      - Files are found recursively inside each batch directory.

    discovery='recursive_leaf':
      - Recursively find directories that directly contain at least one file
        matching the pattern.
      - Process each such directory as its own batch.
      - Files are gathered non-recursively within that batch directory so files
        are not double-counted through both parent and child directories.
    """
    batches = []

    if discovery not in {"immediate", "recursive_leaf"}:
        raise ValueError("discovery must be either 'immediate' or 'recursive_leaf'")

    for d in scan_dirs:
        input_dir = Path(d)
        if not input_dir.is_dir():
            print(f"Warning: Directory '{d}' does not exist or is not a folder. Skipping.")
            continue

        if discovery == "recursive_leaf":
            # A batch is any directory, at any depth, that directly contains
            # matching output files. This is useful for layouts such as:
            # main/group/chunk_0000/*out.pdbqt
            # main/group/chunk_0001/*out.pdbqt
            for batch_dir in sorted([p for p in input_dir.rglob("*") if p.is_dir()] + [input_dir]):
                files = [str(f) for f in sorted(batch_dir.glob(pattern)) if f.is_file()]
                if files:
                    batch_name = _relative_batch_name(input_dir, batch_dir)
                    batches.append((batch_name, batch_dir, files))

            if not any(str(Path(f).resolve()).startswith(str(input_dir.resolve())) for _, _, files in batches for f in files):
                print(f"Warning: No files matching '{pattern}' found under '{input_dir}'. Skipping.")
            continue

        # Default: immediate subdirectories as batches.
        subdirs = sorted([p for p in input_dir.iterdir() if p.is_dir()])

        if subdirs:
            for subdir in subdirs:
                files = [str(f) for f in subdir.rglob(pattern)]
                if files:
                    batches.append((subdir.name, subdir, files))
                else:
                    print(f"Warning: No files matching '{pattern}' found in batch '{subdir}'. Skipping.")
        else:
            files = [str(f) for f in input_dir.rglob(pattern)]
            if files:
                batches.append((input_dir.name, input_dir, files))
            else:
                print(f"Warning: No files matching '{pattern}' found in '{input_dir}'. Skipping.")

    return batches


def safe_parquet_name(batch_name):
    """
    Create a filesystem-safe parquet filename for a batch.
    """
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", batch_name).strip("_")
    return f"{cleaned or 'batch'}.parquet"


def process_batches(scan_dirs, pattern, out_dir, cpus, batch_parquet_dir, keep_batch_parquets, batch_discovery):
    """
    Process subfolders as independent batches, write batch parquets,
    merge them into one master dataframe, and optionally delete batch parquets.
    """
    batches = get_batches_from_input_dirs(scan_dirs, pattern, discovery=batch_discovery)

    if not batches:
        print(f"No files matching '{pattern}' found in any batch.")
        return pd.DataFrame(columns=['Ligand', 'Vina_Score', 'File_Path']), 0, 0, []

    batch_dir = out_dir / batch_parquet_dir
    batch_dir.mkdir(parents=True, exist_ok=True)

    batch_dfs = []
    batch_parquet_paths = []
    total_files = 0
    total_skipped = 0

    print(f"--- Processing {len(batches)} batch(es) ---")

    for batch_name, batch_path, files in batches:
        print(f"\n=== Batch: {batch_name} ===")
        print(f"Source: {batch_path}")

        batch_df, n_files, n_skipped = parse_vina_files(
            files,
            cpus=cpus,
            label=f"Parsing {batch_name}"
        )

        total_files += n_files
        total_skipped += n_skipped

        if batch_df.empty:
            print(f"Warning: No valid Vina scores extracted for batch '{batch_name}'.")
            continue

        batch_df.insert(0, 'Batch', batch_name)
        batch_df.insert(1, 'Batch_Dir', str(batch_path.resolve()))

        batch_parquet_path = batch_dir / safe_parquet_name(batch_name)
        batch_df.to_parquet(batch_parquet_path, index=False)
        batch_parquet_paths.append(batch_parquet_path)
        batch_dfs.append(batch_df)
        print(f"[+] Saved batch parquet ({len(batch_df)} ligands): {batch_parquet_path}")

    if not batch_dfs:
        return pd.DataFrame(columns=['Batch', 'Batch_Dir', 'Ligand', 'Vina_Score', 'File_Path']), total_files, total_skipped, batch_parquet_paths

    merged_df = pd.concat(batch_dfs, ignore_index=True)
    merged_df = merged_df.sort_values(by='Vina_Score', ascending=True).reset_index(drop=True)

    if not keep_batch_parquets:
        for p in batch_parquet_paths:
            try:
                p.unlink()
            except OSError as e:
                print(f"Warning: Could not remove batch parquet '{p}': {e}")
        try:
            batch_dir.rmdir()
        except OSError:
            pass
        print("[+] Removed temporary batch parquets after creating the master parquet.")
    else:
        print(f"[+] Kept {len(batch_parquet_paths)} batch parquet file(s) in: {batch_dir}")

    return merged_df, total_files, total_skipped, batch_parquet_paths


def main():
    global pd
    parser = argparse.ArgumentParser(description="Extract Vina docking scores to Parquet/CSV and export SDFs.")

    # Inputs
    parser.add_argument('-i', '--input_dirs', nargs='+', default=[],
                        help="One or more directories to scan for docking outputs.")
    parser.add_argument('-f', '--input_file', type=str, default=None,
                        help="Text file containing directory paths to scan (one path per line).")
    parser.add_argument('-p', '--pattern', type=str, default="*out.pdbqt",
                        help="File pattern to match Vina outputs (default: '*out.pdbqt').")

    # Outputs
    parser.add_argument('--out_dir', type=str, default=".",
                        help="Directory to save output files (default: current directory).")
    parser.add_argument('--out_parquet', type=str, default="docking_summary.parquet",
                        help="Filename for the full/master Parquet database.")
    parser.add_argument('--out_csv', type=str, default=None,
                        help="Filename for the lightweight top-hits CSV (e.g., 'top_hits.csv').")
    parser.add_argument('--json_out', type=str, default=None,
                        help=("Path for an agent-readable JSON summary. If omitted, writes "
                              "dockflow_vina_extraction_summary.json inside --out_dir."))

    # Independent Export Limits
    parser.add_argument('-t', '--top_n', type=int, default=10,
                        help="Number of ligands to display in the terminal (default: 10).")
    parser.add_argument('--csv_n', type=int, default=1000,
                        help="Number of ligands to save in the CSV (default: 1000).")
    parser.add_argument('--sdf_n', type=int, default=10,
                        help="Number of ligands to export to individual SDFs (default: 10).")

    # Processing
    parser.add_argument('-c', '--cpus', type=int, default=None,
                        help="Number of CPU cores to use.")
    parser.add_argument('--batch_by_subdirs', action='store_true',
                        help=("Process subdirectories as separate batches, save batch parquet(s), "
                              "then merge into the master parquet."))
    parser.add_argument('--batch_discovery', choices=['immediate', 'recursive_leaf'], default='immediate',
                        help=("How to define batches when --batch_by_subdirs is used. "
                              "'immediate' processes each immediate subdirectory as a batch and scans "
                              "recursively inside it. 'recursive_leaf' recursively finds directories that "
                              "directly contain matching files and processes each as its own batch. "
                              "Default: immediate."))
    parser.add_argument('--batch_parquet_dir', type=str, default="batch_parquets",
                        help="Subdirectory inside --out_dir for batch parquet files (default: 'batch_parquets').")
    parser.add_argument('--keep_batch_parquets', action='store_true',
                        help="Keep per-batch parquet files after the master parquet is created. Default: remove them.")

    # SDF Export Arguments
    parser.add_argument('--export_sdf', action='store_true',
                        help="Export the top hits to individual SDF files (best pose only).")
    parser.add_argument('--sdf_dir', type=str, default="top_sdfs",
                        help="Subdirectory name inside --out_dir to save the SDF files (default: 'top_sdfs').")

    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.json_out is None:
        args.json_out = str(out_dir / "dockflow_vina_extraction_summary.json")

    try:
        import pandas as pandas_module
    except ImportError as exc:
        summary = build_json_summary(args, status="failed", error=str(exc))
        write_json_summary(summary, args.json_out)
        raise
    pd = pandas_module

    # Validate inputs
    if not args.input_dirs and not args.input_file:
        parser.error("You must provide directories to scan using either -i/--input_dirs or -f/--input_file.")

    scan_dirs = list(args.input_dirs)
    if args.input_file:
        try:
            with open(args.input_file) as f:
                for line in f:
                    clean_path = line.strip()
                    if clean_path and not clean_path.startswith('#'):
                        scan_dirs.append(clean_path)
        except FileNotFoundError:
            print(f"Error: The input file '{args.input_file}' was not found.")
            summary = build_json_summary(args, status="failed", error=f"Input file not found: {args.input_file}")
            write_json_summary(summary, args.json_out)
            return

    parquet_path = out_dir / args.out_parquet
    csv_path = out_dir / args.out_csv if args.out_csv else None

    # Check SDF exporter once, before doing expensive parsing/conversion work.
    if args.export_sdf and not check_mk_export_available():
        print("SDF export is disabled for this run. Score extraction will continue.")
        args.export_sdf = False

    if args.batch_by_subdirs:
        df, total_files, total_skipped, _ = process_batches(
            scan_dirs=scan_dirs,
            pattern=args.pattern,
            out_dir=out_dir,
            cpus=args.cpus,
            batch_parquet_dir=args.batch_parquet_dir,
            keep_batch_parquets=args.keep_batch_parquets,
            batch_discovery=args.batch_discovery
        )
    else:
        # 1. Gather files
        print("--- Scanning directories ---")
        all_files = gather_files(scan_dirs, args.pattern)

        if not all_files:
            print(f"No files matching '{args.pattern}' found.")
            summary = build_json_summary(args, status="failed", error=f"No files matching '{args.pattern}' found.")
            write_json_summary(summary, args.json_out)
            return

        # 2. Extract Data
        df, total_files, total_skipped = parse_vina_files(
            all_files,
            cpus=args.cpus,
            label="Parsing PDBQTs"
        )

    if df.empty:
        print("No valid Vina scores extracted.")
        summary = build_json_summary(
            args,
            status="failed",
            total_files=total_files,
            total_skipped=total_skipped,
            df=df,
            error="No valid Vina scores extracted.",
        )
        write_json_summary(summary, args.json_out)
        return

    # 3. Save Full/Master Dataset
    df = df.sort_values(by='Vina_Score', ascending=True).reset_index(drop=True)
    try:
        df.to_parquet(parquet_path, index=False)
    except Exception as exc:
        summary = build_json_summary(
            args,
            status="failed",
            total_files=total_files,
            total_skipped=total_skipped,
            df=df,
            error=str(exc),
        )
        write_json_summary(summary, args.json_out)
        raise
    print(f"\n[+] Saved master database ({len(df)} ligands): {parquet_path}")
    print(f"[+] Final parse summary: {len(df)} scored files, {total_skipped} skipped files, {total_files} total files scanned.")

    # 4. Save Optional CSV
    if args.out_csv:
        df.head(args.csv_n).to_csv(csv_path, index=False)
        print(f"[+] Saved top {min(len(df), args.csv_n)} ligands to CSV: {csv_path}")

    # 5. Display Clean Terminal Output
    print(f"\n--- Top {args.top_n} Ligands ---")
    pd.set_option('display.max_colwidth', 80)
    print(df.head(args.top_n).to_string(index=False))

    # 6. Optional: Export top poses to SDF
    if args.export_sdf:
        sdf_target_dir = out_dir / args.sdf_dir
        sdf_target_dir.mkdir(parents=True, exist_ok=True)

        # Grab only the requested number for SDF export
        sdf_df = df.head(args.sdf_n)

        print(f"\n--- Exporting top {len(sdf_df)} poses to SDF ---")
        print(f"Destination: {sdf_target_dir}")

        success_count = 0
        for _, row in tqdm(sdf_df.iterrows(), total=len(sdf_df), desc="Converting", unit="mol"):
            if extract_and_convert_top_pose(row['File_Path'], sdf_target_dir, row['Ligand']):
                success_count += 1

        print(f"[+] Successfully converted {success_count}/{len(sdf_df)} ligands to SDF.")

    else:
        sdf_df = []
        success_count = 0

    summary = build_json_summary(
        args,
        status="success",
        total_files=total_files,
        total_skipped=total_skipped,
        df=df,
        csv_path=csv_path,
        sdf_attempted=len(sdf_df),
        sdf_success=success_count,
    )
    write_json_summary(summary, args.json_out)
    print(f"[+] Saved agent-readable JSON summary: {args.json_out}")


if __name__ == "__main__":
    main()
