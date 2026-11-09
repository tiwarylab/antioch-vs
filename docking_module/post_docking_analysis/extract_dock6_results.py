#!/usr/bin/env python3
"""
Extract DOCK6 screening results from scored Mol2 files.

Recursively searches for DOCK6 scored Mol2 output files, extracts the metadata
written before each Mol2 block, stores all scored poses in a Parquet table, and
optionally exports top hits to CSV and SDF.

Optional batch mode mirrors the Vina extractor:
    Input directories can be split into batches either by immediate subfolders
    or by recursively discovered leaf folders that directly contain matching
    scored Mol2 files. Each batch is processed independently into a batch
    parquet first, then all batches are merged into one master parquet.
"""

import argparse
import json
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **_kwargs):
        return iterable


DOCK6_META_REGEX = re.compile(r"^##########\s*(?P<key>[^:]+):\s*(?P<value>.*?)\s*$")
MOLECULE_MARKER = "@<TRIPOS>MOLECULE"
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


def build_json_summary(args, status, total_files=0, total_skipped=0, df=None, rank_df=None,
                       csv_path=None, export_attempted=0, export_success=0, error=None):
    out_dir = Path(args.out_dir).resolve()
    parquet_path = out_dir / args.out_parquet
    csv_target = csv_path or (out_dir / args.out_csv if args.out_csv else None)
    ranking_df = rank_df if rank_df is not None else df

    metrics = {
        "files_scanned": int(total_files),
        "scored_rows": int(len(df)) if df is not None else 0,
        "ranked_rows": int(len(ranking_df)) if ranking_df is not None else 0,
        "skipped_files": int(total_skipped),
        "unique_ligand_count": int(df["Ligand"].nunique()) if df is not None and "Ligand" in df.columns else 0,
        "scores_parquet_exists": parquet_path.exists(),
        "top_hits_csv_exists": csv_target.exists() if csv_target else False,
        "pose_export_attempted": int(export_attempted),
        "pose_export_success": int(export_success),
    }

    if df is not None and not df.empty:
        for column in ("Grid_Score", "Total_Score"):
            if column in df.columns:
                metrics.update({
                    "score_column": column,
                    "best_score": float(df[column].min()),
                    "median_score": float(df[column].median()),
                })
                break

    top_hits_preview = []
    if ranking_df is not None and not ranking_df.empty:
        top_hits_preview = ranking_df.head(min(args.top_n, 10)).to_dict(orient="records")

    return {
        "stage": "post_docking_analysis",
        "extractor": "dock6",
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
            "best_per_ligand": args.best_per_ligand,
            "export_sdf": args.export_sdf,
            "export_mol2": args.export_mol2,
            "top_n": args.top_n,
            "csv_n": args.csv_n,
            "sdf_n": args.sdf_n,
        },
        "outputs": {
            "out_dir": str(out_dir),
            "scores_parquet": str(parquet_path),
            "top_hits_csv": str(csv_target) if csv_target else None,
            "pose_export_dir": str(out_dir / args.sdf_dir) if args.export_sdf or args.export_mol2 else None,
            "json_summary": str(Path(args.json_out).resolve()),
        },
        "metrics": metrics,
        "top_hits_preview": top_hits_preview,
        "warnings": [],
        "error": error,
    }


def normalize_key(key):
    return re.sub(r"[^A-Za-z0-9_]+", "_", key.strip()).strip("_")


def parse_value(value):
    value = value.strip()
    if value == "":
        return value

    try:
        number = float(value)
    except ValueError:
        return value

    if number.is_integer():
        return int(number)
    return number


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
    return sorted(all_files)


def parse_dock6_file(file_path):
    """
    Parse one DOCK6 scored Mol2 file into one row per scored pose.
    """
    path_obj = Path(file_path)
    rows = []
    metadata = {}
    pose_index = 0
    current_start_line = None
    current_mol_name = None
    previous_pose = None
    pending_name_line = False

    try:
        with path_obj.open("r", errors="replace") as handle:
            for line_number, line in enumerate(handle, start=1):
                meta_match = DOCK6_META_REGEX.match(line)
                if meta_match:
                    key = normalize_key(meta_match.group("key"))
                    value = parse_value(meta_match.group("value"))
                    metadata[key] = value
                    continue

                if line.startswith(MOLECULE_MARKER):
                    if previous_pose is not None:
                        previous_pose["Mol2_End_Line"] = line_number - 1
                        rows.append(previous_pose)

                    pose_index += 1
                    current_start_line = line_number
                    current_mol_name = None
                    pending_name_line = True
                    previous_pose = None
                    continue

                if pending_name_line:
                    current_mol_name = line.strip()
                    ligand = str(metadata.get("Name") or current_mol_name or path_obj.stem)
                    row = {
                        "Ligand": ligand,
                        "Pose_Index": pose_index,
                        "Source_File": str(path_obj.resolve()),
                        "Mol2_Start_Line": current_start_line,
                        "Mol2_End_Line": None,
                    }
                    row.update(metadata)
                    if "Name" not in row and current_mol_name:
                        row["Name"] = current_mol_name
                    previous_pose = row
                    metadata = {}
                    pending_name_line = False

            if previous_pose is not None:
                previous_pose["Mol2_End_Line"] = line_number
                rows.append(previous_pose)
    except Exception:
        return []

    return rows


def _extract_data_worker(file_path):
    return parse_dock6_file(file_path)


def parse_dock6_files(files, cpus=None, label="Parsing DOCK6 Mol2 files"):
    """
    Parse a list of DOCK6 scored Mol2 files and return a sorted dataframe.
    """
    if not files:
        return pd.DataFrame(), 0, 0

    print(f"Found {len(files)} files. Parsing for DOCK6 scores...")
    results = []
    n_procs = cpus or mp.cpu_count()

    with mp.Pool(processes=n_procs) as pool:
        iterator = pool.imap_unordered(_extract_data_worker, files)
        for file_rows in tqdm(iterator, total=len(files), desc=label, unit="file"):
            results.extend(file_rows)

    skipped = len(files) - sum(1 for row_group in results_by_file(files, results) if row_group)
    print(f"[+] Extracted {len(results)} scored poses from {len(files)} total files.")
    print(f"[+] Skipped {skipped} files with no valid DOCK6 scored Mol2 poses or unreadable content.")

    if not results:
        return pd.DataFrame(), len(files), skipped

    df = pd.DataFrame(results)
    df = sort_dock6_dataframe(df)
    return df, len(files), skipped


def results_by_file(files, rows):
    found = {row["Source_File"] for row in rows}
    for file_path in files:
        resolved = str(Path(file_path).resolve())
        yield resolved in found


def sort_dock6_dataframe(df, score_column="Grid_Score"):
    """
    Sort from best to worst score. For DOCK6 grid scores, lower is better.
    """
    if df.empty:
        return df

    if score_column in df.columns:
        df = df.sort_values(by=score_column, ascending=True, na_position="last")
    elif "Total_Score" in df.columns:
        df = df.sort_values(by="Total_Score", ascending=True, na_position="last")
    return df.reset_index(drop=True)


def best_pose_per_ligand(df):
    """
    Return the best-scoring pose for each ligand from an already sorted table.
    """
    if df.empty or "Ligand" not in df.columns:
        return df
    return df.drop_duplicates(subset=["Ligand"], keep="first").reset_index(drop=True)


def _relative_batch_name(input_dir, batch_dir):
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
            for batch_dir in sorted([p for p in input_dir.rglob("*") if p.is_dir()] + [input_dir]):
                files = [str(f) for f in sorted(batch_dir.glob(pattern)) if f.is_file()]
                if files:
                    batch_name = _relative_batch_name(input_dir, batch_dir)
                    batches.append((batch_name, batch_dir, files))
            continue

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
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", batch_name).strip("_")
    return f"{cleaned or 'batch'}.parquet"


def process_batches(scan_dirs, pattern, out_dir, cpus, batch_parquet_dir, keep_batch_parquets, batch_discovery):
    batches = get_batches_from_input_dirs(scan_dirs, pattern, discovery=batch_discovery)

    if not batches:
        print(f"No files matching '{pattern}' found in any batch.")
        return pd.DataFrame(), 0, 0, []

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

        batch_df, n_files, n_skipped = parse_dock6_files(
            files,
            cpus=cpus,
            label=f"Parsing {batch_name}"
        )

        total_files += n_files
        total_skipped += n_skipped

        if batch_df.empty:
            print(f"Warning: No valid DOCK6 scores extracted for batch '{batch_name}'.")
            continue

        batch_df.insert(0, "Batch", batch_name)
        batch_df.insert(1, "Batch_Dir", str(batch_path.resolve()))

        batch_parquet_path = batch_dir / safe_parquet_name(batch_name)
        try:
            write_parquet_or_raise(batch_df, batch_parquet_path)
        except RuntimeError as exc:
            print(f"Error: {exc}")
            raise
        batch_parquet_paths.append(batch_parquet_path)
        batch_dfs.append(batch_df)
        print(f"[+] Saved batch parquet ({len(batch_df)} poses): {batch_parquet_path}")

    if not batch_dfs:
        return pd.DataFrame(), total_files, total_skipped, batch_parquet_paths

    merged_df = pd.concat(batch_dfs, ignore_index=True)
    merged_df = sort_dock6_dataframe(merged_df)

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


def extract_mol2_pose(row, output_path):
    """
    Extract one Mol2 pose block from a scored DOCK6 multi-mol2 file.
    """
    source_file = Path(row["Source_File"])
    start_line = int(row["Mol2_Start_Line"])
    end_line = int(row["Mol2_End_Line"])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with source_file.open("r", errors="replace") as handle, output_path.open("w") as out:
        for line_number, line in enumerate(handle, start=1):
            if line_number < start_line:
                continue
            if line_number > end_line:
                break
            out.write(line)


def safe_pose_name(row):
    ligand = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(row["Ligand"])).strip("_") or "ligand"
    source = re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        Path(row["Source_File"]).stem,
    ).strip("_") or "source"
    pose_index = int(row["Pose_Index"])
    return f"{ligand}_{source}_pose{pose_index:06d}"


def write_parquet_or_raise(df, parquet_path):
    try:
        df.to_parquet(parquet_path, index=False)
    except ImportError as exc:
        raise RuntimeError(
            "Could not write Parquet because pandas cannot find a Parquet engine. "
            "Install pyarrow or fastparquet in the active environment, for example: "
            "pip install pyarrow"
        ) from exc


def check_obabel_available(command):
    if shutil.which(command):
        return True

    print(f"\nError: '{command}' was not found on your PATH.")
    print("To export Dock6 Mol2 poses to SDF, activate an environment with Open Babel, e.g.:")
    print("    conda activate obabel_env")
    print("or install Open Babel in the active environment.")
    return False


def convert_mol2_to_sdf(input_mol2, output_sdf, obabel_command="obabel"):
    try:
        subprocess.run(
            [obabel_command, str(input_mol2), "-O", str(output_sdf)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return output_sdf.exists()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def export_top_poses(df, output_dir, n_poses, export_sdf=False, keep_mol2=False, obabel_command="obabel"):
    export_df = df.head(n_poses)
    output_dir.mkdir(parents=True, exist_ok=True)

    success = 0
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        for _, row in tqdm(export_df.iterrows(), total=len(export_df), desc="Exporting", unit="pose"):
            pose_name = safe_pose_name(row)
            mol2_target = output_dir / f"{pose_name}.mol2" if keep_mol2 else tmpdir / f"{pose_name}.mol2"

            try:
                extract_mol2_pose(row, mol2_target)
            except Exception as exc:
                print(f"Warning: Could not extract Mol2 pose {pose_name}: {exc}")
                continue

            if export_sdf:
                sdf_target = output_dir / f"{pose_name}.sdf"
                if convert_mol2_to_sdf(mol2_target, sdf_target, obabel_command=obabel_command):
                    success += 1
                else:
                    print(f"Warning: SDF conversion failed for {pose_name}")
            else:
                success += 1

    return success, len(export_df)


def main():
    global pd
    parser = argparse.ArgumentParser(description="Extract DOCK6 docking scores to Parquet/CSV and export top poses.")

    parser.add_argument("-i", "--input_dirs", nargs="+", default=[],
                        help="One or more directories to scan for DOCK6 scored Mol2 outputs.")
    parser.add_argument("-f", "--input_file", type=str, default=None,
                        help="Text file containing directory paths to scan (one path per line).")
    parser.add_argument("-p", "--pattern", type=str, default="*_flex_scored*.mol2",
                        help="File pattern to match DOCK6 outputs (default: '*_flex_scored*.mol2').")

    parser.add_argument("--out_dir", type=str, default=".",
                        help="Directory to save output files (default: current directory).")
    parser.add_argument("--out_parquet", type=str, default="dock6_summary.parquet",
                        help="Filename for the full/master Parquet database.")
    parser.add_argument("--out_csv", type=str, default=None,
                        help="Filename for the lightweight top-hits CSV (e.g., 'top_hits.csv').")
    parser.add_argument("--json_out", type=str, default=None,
                        help=("Path for an agent-readable JSON summary. If omitted, writes "
                              "dockflow_dock6_extraction_summary.json inside --out_dir."))

    parser.add_argument("-t", "--top_n", type=int, default=10,
                        help="Number of poses to display in the terminal (default: 10).")
    parser.add_argument("--csv_n", type=int, default=1000,
                        help="Number of poses to save in the CSV (default: 1000).")
    parser.add_argument("--sdf_n", type=int, default=10,
                        help="Number of poses to export to SDF/Mol2 (default: 10).")
    parser.add_argument("--best_per_ligand", action="store_true",
                        help=("Use only the best-scoring pose per ligand for terminal display, CSV, "
                              "and SDF/Mol2 export. The master parquet still keeps all scored poses."))

    parser.add_argument("-c", "--cpus", type=int, default=None,
                        help="Number of CPU cores to use while parsing.")
    parser.add_argument("--batch_by_subdirs", action="store_true",
                        help=("Process subdirectories as separate batches, save batch parquet(s), "
                              "then merge into the master parquet."))
    parser.add_argument("--batch_discovery", choices=["immediate", "recursive_leaf"], default="immediate",
                        help=("How to define batches when --batch_by_subdirs is used. "
                              "'immediate' processes each immediate subdirectory as a batch and scans "
                              "recursively inside it. 'recursive_leaf' recursively finds directories that "
                              "directly contain matching files and processes each as its own batch. "
                              "Default: immediate."))
    parser.add_argument("--batch_parquet_dir", type=str, default="batch_parquets",
                        help="Subdirectory inside --out_dir for batch parquet files (default: 'batch_parquets').")
    parser.add_argument("--keep_batch_parquets", action="store_true",
                        help="Keep per-batch parquet files after the master parquet is created. Default: remove them.")

    parser.add_argument("--export_sdf", action="store_true",
                        help="Export top poses to individual SDF files.")
    parser.add_argument("--export_mol2", action="store_true",
                        help="Export top poses to individual Mol2 files.")
    parser.add_argument("--sdf_dir", type=str, default="top_sdfs",
                        help="Subdirectory name inside --out_dir to save exported poses (default: 'top_sdfs').")
    parser.add_argument("--obabel_command", type=str, default="obabel",
                        help="Open Babel command for Mol2 to SDF conversion (default: obabel).")

    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.json_out is None:
        args.json_out = str(out_dir / "dockflow_dock6_extraction_summary.json")

    try:
        import pandas as pandas_module
    except ImportError as exc:
        summary = build_json_summary(args, status="failed", error=str(exc))
        write_json_summary(summary, args.json_out)
        raise
    pd = pandas_module

    if not args.input_dirs and not args.input_file:
        parser.error("You must provide directories to scan using either -i/--input_dirs or -f/--input_file.")

    scan_dirs = list(args.input_dirs)
    if args.input_file:
        try:
            with open(args.input_file) as handle:
                for line in handle:
                    clean_path = line.strip()
                    if clean_path and not clean_path.startswith("#"):
                        scan_dirs.append(clean_path)
        except FileNotFoundError:
            print(f"Error: The input file '{args.input_file}' was not found.")
            summary = build_json_summary(args, status="failed", error=f"Input file not found: {args.input_file}")
            write_json_summary(summary, args.json_out)
            return

    parquet_path = out_dir / args.out_parquet
    csv_path = out_dir / args.out_csv if args.out_csv else None

    if args.export_sdf and not check_obabel_available(args.obabel_command):
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
            batch_discovery=args.batch_discovery,
        )
    else:
        print("--- Scanning directories ---")
        all_files = gather_files(scan_dirs, args.pattern)
        if not all_files:
            print(f"No files matching '{args.pattern}' found.")
            summary = build_json_summary(args, status="failed", error=f"No files matching '{args.pattern}' found.")
            write_json_summary(summary, args.json_out)
            return

        df, total_files, total_skipped = parse_dock6_files(
            all_files,
            cpus=args.cpus,
            label="Parsing DOCK6 Mol2 files",
        )

    if df.empty:
        print("No valid DOCK6 scores extracted.")
        summary = build_json_summary(
            args,
            status="failed",
            total_files=total_files,
            total_skipped=total_skipped,
            df=df,
            error="No valid DOCK6 scores extracted.",
        )
        write_json_summary(summary, args.json_out)
        return

    df = sort_dock6_dataframe(df)
    try:
        write_parquet_or_raise(df, parquet_path)
    except RuntimeError as exc:
        print(f"Error: {exc}")
        summary = build_json_summary(
            args,
            status="failed",
            total_files=total_files,
            total_skipped=total_skipped,
            df=df,
            error=str(exc),
        )
        write_json_summary(summary, args.json_out)
        return
    print(f"\n[+] Saved master database ({len(df)} poses): {parquet_path}")
    print(f"[+] Final parse summary: {len(df)} scored poses, {total_skipped} skipped files, {total_files} total files scanned.")

    rank_df = best_pose_per_ligand(df) if args.best_per_ligand else df
    rank_label = "ligands" if args.best_per_ligand else "poses"
    if args.best_per_ligand:
        print(f"[+] Collapsed to best pose per ligand for top-hit outputs: {len(rank_df)} ligands.")

    if args.out_csv:
        rank_df.head(args.csv_n).to_csv(csv_path, index=False)
        print(f"[+] Saved top {min(len(rank_df), args.csv_n)} {rank_label} to CSV: {csv_path}")

    print(f"\n--- Top {args.top_n} Dock6 {rank_label.title()} ---")
    pd.set_option("display.max_colwidth", 80)
    display_columns = [
        col for col in [
            "Ligand",
            "Grid_Score",
            "Grid_vdw_energy",
            "Grid_es_energy",
            "HA_RMSDs",
            "Pose_Index",
            "Source_File",
        ]
        if col in rank_df.columns
    ]
    print(rank_df.head(args.top_n)[display_columns].to_string(index=False))

    if args.export_sdf or args.export_mol2:
        pose_target_dir = out_dir / args.sdf_dir
        print(f"\n--- Exporting top {args.sdf_n} poses ---")
        print(f"Destination: {pose_target_dir}")

        success_count, attempted = export_top_poses(
            df=rank_df,
            output_dir=pose_target_dir,
            n_poses=args.sdf_n,
            export_sdf=args.export_sdf,
            keep_mol2=args.export_mol2,
            obabel_command=args.obabel_command,
        )

        export_label = "SDF" if args.export_sdf else "Mol2"
        print(f"[+] Successfully exported {success_count}/{attempted} poses to {export_label}.")
    else:
        success_count = 0
        attempted = 0

    summary = build_json_summary(
        args,
        status="success",
        total_files=total_files,
        total_skipped=total_skipped,
        df=df,
        rank_df=rank_df,
        csv_path=csv_path,
        export_attempted=attempted,
        export_success=success_count,
    )
    write_json_summary(summary, args.json_out)
    print(f"[+] Saved agent-readable JSON summary: {args.json_out}")


if __name__ == "__main__":
    main()
