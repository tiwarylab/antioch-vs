#!/usr/bin/env python3
"""Collect the original docking pose files for the top fraction of a consensus ranking.

Pipeline context
----------------
This is the traceback step that runs AFTER ``build_unweighted_rrf_consensus.py``.
It answers: "take the top X% of the fused ranking, and for each of those ligands
pull back the single best docking pose across the N receptors, then put those
pose files in one place."

Labels are optional. A benchmark consensus carries ``Class`` and the manifest
reports it; an unlabeled screen has no ``Class`` column at all, and the manifest
simply omits it rather than carrying a column of blanks.

What it does
------------
1. Reads ``consensus_all.parquet`` and selects the top fraction by ``Consensus_Rank``
   (the fused ranking over every ligand, actives and decoys alike where labelled).
2. For each selected ligand, looks across every receptor's ``enrichment_ranked.parquet``
   and picks the *best docking pose*:
     - default ``--best_by score``: the receptor/pose with the lowest raw docking
       score (Grid_Score / Vina_Score; use --higher_is_better to flip). This is the
       physically best-scoring pose and is generally what "best pose" means.
     - ``--best_by rank``: the receptor where the ligand achieved the best within-
       receptor rank (this is what contributed most to the consensus, and can differ
       from the best raw score when receptors have different score distributions).
3. Materializes each best pose into ``--out_dir/poses/``:
     - If the receptor table has per-pose line offsets (Mol2_Start_Line/Mol2_End_Line),
       the pose is EXTRACTED by line range from its shared Source_File .mol2.
     - For Vina multi-model PDBQT files, the best-scoring MODEL is written by
       default. This is the same model used by ``extract_vina_results.py``
       (Vina writes its best pose first).
4. Writes ``collection_manifest.csv`` with full provenance for every collected pose.

The traceback tables are the same ``enrichment_ranked.parquet`` files the consensus
was built from, so the selection is guaranteed consistent with the consensus.

Assumptions worth checking once against your data
-------------------------------------------------
* Mol2 line offsets are 1-based inclusive (lines[start-1:end]). Override with
  --mol2_line_base 0 if your extractor wrote 0-based offsets.
* Ligand IDs in consensus_all match the ``--ligand_column`` in enrichment_ranked
  (both are the (optionally -again stripped) names the consensus used).
* This script does not create AMBER inputs itself. Use
  ``stage_consensus_ctmd_pairs.py`` after collection to make one CTMD-ready
  project per receptor-ligand pair. The receptor template supplied here must
  therefore point to the PDB structure used for the chosen receptor.

Example (SLURM-friendly, one engine)
------------------------------------
    python collect_top_consensus_poses.py \
        --consensus_all  dock6_selected_af2_unweighted_rrf_consensus/consensus_all.parquet \
        --input rmsa03_r605_m2_s080=dock6/rmsa03_r605_m2_s080/dude_enrichment \
        --input rmsa08_r251_m1_s080=dock6/rmsa08_r251_m1_s080/dude_enrichment \
        --input rmsa08_r336_m2_s122=dock6/rmsa08_r336_m2_s122/dude_enrichment \
        --input rmsa08_r403_m1_s094=dock6/rmsa08_r403_m1_s094/dude_enrichment \
        --input rmsa08_r499_m1_s070=dock6/rmsa08_r499_m1_s070/dude_enrichment \
        --out_dir  dock6_top1pct_poses \
        --top_fraction 0.01
"""
from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

SCORE_CANDIDATES = ("Grid_Score", "Vina_Score", "Total_Score", "Score", "score")
# Columns that point at the original pose file, in order of preference.
FILE_COLUMN_CANDIDATES = ("Source_File", "File_Path", "Pose_File", "Path", "Source_Parquet")
START_LINE_CANDIDATES = ("Mol2_Start_Line", "Start_Line")
END_LINE_CANDIDATES = ("Mol2_End_Line", "End_Line")


# --------------------------------------------------------------------------- #
# argument parsing
# --------------------------------------------------------------------------- #
def parse_named_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("Each --input must use NAME=PATH.")
    name, path_text = value.split("=", 1)
    name, path_text = name.strip(), path_text.strip()
    if not name or not path_text:
        raise argparse.ArgumentTypeError("Both NAME and PATH are required.")
    return name, Path(path_text).expanduser()


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.\-]+", "_", str(value)).strip("_") or "ligand"


def resolve_ranked_parquet(path: Path) -> Path:
    """Accept an enrichment output dir or a parquet file directly."""
    if path.is_file():
        if path.suffix != ".parquet":
            raise ValueError(f"Input is not a Parquet file: {path}")
        return path
    for candidate in ("enrichment_ranked.parquet", "consensus_all.parquet"):
        hit = path / candidate
        if hit.is_file():
            return hit
    raise FileNotFoundError(f"No enrichment_ranked.parquet found under: {path}")


def choose_column(frame: pd.DataFrame, requested: str | None,
                  candidates: tuple[str, ...], kind: str) -> str:
    if requested:
        if requested not in frame.columns:
            raise ValueError(f"{kind} column {requested!r} not found. Have: {list(frame.columns)}")
        return requested
    for column in candidates:
        if column in frame.columns:
            return column
    raise ValueError(f"Could not auto-detect a {kind} column. Have: {list(frame.columns)}")


def first_present(frame: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    for column in candidates:
        if column in frame.columns:
            return column
    return None


# --------------------------------------------------------------------------- #
# core selection
# --------------------------------------------------------------------------- #
def load_top_set(consensus_all: Path, top_fraction: float,
                 ligand_column: str) -> pd.DataFrame:
    con = pd.read_parquet(consensus_all)
    if ligand_column not in con.columns:
        # consensus_all stores the ligand as a column named "Ligand" after reset_index
        raise ValueError(f"{ligand_column!r} not in consensus_all columns: {list(con.columns)}")
    if "Consensus_Rank" not in con.columns:
        if "Consensus_RRF_Score" not in con.columns:
            raise ValueError("consensus_all needs Consensus_Rank or Consensus_RRF_Score.")
        con = con.sort_values("Consensus_RRF_Score", ascending=False).reset_index(drop=True)
        con["Consensus_Rank"] = con.index + 1

    total = len(con)
    cutoff = max(1, math.ceil(total * top_fraction))
    top = con.sort_values("Consensus_Rank").head(cutoff).copy()
    keep = [c for c in (ligand_column, "Class", "Consensus_RRF_Score", "Consensus_Rank",
                        "Best_Receptor_Rank", "Top_Support_Votes", "Supporting_Receptors")
            if c in top.columns]
    return top[keep], total, cutoff


def load_receptor_tables(inputs: list[tuple[str, Path]], ligand_column: str,
                         score_column: str | None) -> tuple[pd.DataFrame, str]:
    frames = []
    resolved_score = score_column
    for name, path in inputs:
        parquet = resolve_ranked_parquet(path)
        frame = pd.read_parquet(parquet)
        if ligand_column not in frame.columns:
            raise ValueError(f"{name}: ligand column {ligand_column!r} missing in {parquet}")
        sc = choose_column(frame, score_column, SCORE_CANDIDATES, "score")
        if resolved_score is None:
            resolved_score = sc
        elif sc != resolved_score:
            raise ValueError(
                f"{name}: detected score column {sc!r} but earlier receptor used "
                f"{resolved_score!r}. Pass --score_column explicitly to force one.")
        frame = frame.copy()
        frame["__receptor__"] = name
        frame["__score__"] = pd.to_numeric(frame[sc], errors="coerce")
        frame["__rank__"] = pd.to_numeric(frame.get("Rank"), errors="coerce")
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True)
    return combined, resolved_score


def pick_best_pose(rows: pd.DataFrame, best_by: str, lower_is_better: bool) -> pd.Series:
    """From all (receptor, pose) rows for one ligand, return the best one."""
    valid = rows.dropna(subset=["__score__"])
    if valid.empty:
        valid = rows
    if best_by == "rank":
        ranked = valid.dropna(subset=["__rank__"])
        pool = ranked if not ranked.empty else valid
        return pool.loc[pool["__rank__"].idxmin()]
    # best_by == "score"
    idx = valid["__score__"].idxmin() if lower_is_better else valid["__score__"].idxmax()
    return valid.loc[idx]


# --------------------------------------------------------------------------- #
# materialization
# --------------------------------------------------------------------------- #
def _is_meta_line(line: str) -> bool:
    """True for blank lines and DOCK '##########' metadata-comment lines."""
    stripped = line.strip()
    return stripped == "" or stripped.startswith("##########")


def collect_pose_lines(all_lines: list[str], start: int, end: int,
                       line_base: int, header_mode: str) -> list[str]:
    """Return the lines for one pose from an in-memory mol2.

    In DOCK6 mol2 output each molecule's '##########' metadata block sits ABOVE
    its '@<TRIPOS>MOLECULE' line, while Mol2_Start_Line/End_Line span from
    '@<TRIPOS>MOLECULE' to just before the next molecule. So the raw [start, end]
    slice loses this molecule's own header and picks up the NEXT molecule's header
    at the tail. header_mode fixes that:

      * 'reattach' (default): drop the trailing comment block (next molecule's) and
        prepend the '##########' block that sits directly above start (this
        molecule's own header).
      * 'strip': just drop the trailing comment block; do not add a header.
      * 'raw': exact [start, end] slice (original behavior).
    """
    n = len(all_lines)
    s = (start - 1) if line_base == 1 else start
    e = (end - 1) if line_base == 1 else end
    s = max(0, min(s, n - 1))
    e = max(0, min(e, n - 1))
    if header_mode == "raw":
        return all_lines[s:e + 1]

    # trim trailing blank/comment lines: they are the next molecule's header
    body_end = e
    while body_end >= s and _is_meta_line(all_lines[body_end]):
        body_end -= 1
    body = all_lines[s:body_end + 1]
    if header_mode == "strip":
        return body

    # reattach: walk up from start collecting this molecule's own header block
    idx = s - 1
    header: list[str] = []
    while idx >= 0 and _is_meta_line(all_lines[idx]):
        header.append(all_lines[idx])
        idx -= 1
    header.reverse()
    while header and header[0].strip() == "":  # drop leading blanks
        header.pop(0)
    return header + body


def extract_mol2_ranges(source: Path, jobs: list[dict], line_base: int,
                        header_mode: str, dry_run: bool) -> None:
    """Read one shared .mol2 once and write each requested pose (header-aware)."""
    if dry_run:
        return
    with source.open("r", errors="replace") as handle:
        all_lines = handle.readlines()
    for job in jobs:
        lines = collect_pose_lines(all_lines, job["start"], job["end"], line_base, header_mode)
        text = "".join(lines)
        if text and not text.endswith("\n"):
            text += "\n"
        Path(job["dest"]).write_text(text)


VINA_RESULT_RE = re.compile(r"VINA RESULT:\s*(-?\d+(?:\.\d+)?)")


def split_pdbqt_models(lines: list[str]) -> list[tuple[float | None, list[str]]]:
    """Split a (possibly multi-model) pdbqt into (affinity, block) pairs.

    Each block runs from its MODEL line through its ENDMDL line. Affinity is the
    first number on that model's 'REMARK VINA RESULT:' line (None if absent).
    """
    models: list[tuple[float | None, list[str]]] = []
    current: list[str] | None = None
    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith("MODEL"):
            current = [line]
        elif current is not None:
            current.append(line)
            if stripped.startswith("ENDMDL"):
                affinity = None
                for candidate in current:
                    match = VINA_RESULT_RE.search(candidate)
                    if match:
                        affinity = float(match.group(1))
                        break
                models.append((affinity, current))
                current = None
    return models


def extract_best_vina_model(lines: list[str], lower_is_better: bool) -> tuple[float | None, list[str]] | None:
    """Return (affinity, block) for the best-scoring MODEL, or None if no MODEL records."""
    models = split_pdbqt_models(lines)
    if not models:
        return None
    scored = [m for m in models if m[0] is not None]
    pool = scored if scored else models
    best = min(pool, key=lambda m: m[0]) if lower_is_better else max(pool, key=lambda m: m[0])
    affinity, block = best
    # normalize the MODEL header to 1 so the single-pose file is unambiguous
    if block and block[0].lstrip().startswith("MODEL"):
        block = ["MODEL 1\n"] + block[1:]
    return affinity, block


def run_mk_export(mk_cmd: str, pdbqt: Path, sdf: Path) -> tuple[bool, str]:
    """Run Meeko's mk_export.py to write an SDF for a single-pose pdbqt."""
    try:
        result = subprocess.run([mk_cmd, str(pdbqt), "-s", str(sdf)],
                                capture_output=True, text=True)
    except FileNotFoundError:
        return False, f"{mk_cmd} not found on PATH"
    if result.returncode != 0:
        msg = (result.stderr or result.stdout).strip().splitlines()
        return False, (msg[-1] if msg else f"exit {result.returncode}")
    if not sdf.exists():
        return False, "command succeeded but no SDF was written"
    return True, "ok"


def place_file(src: Path, dest: Path, mode: str, dry_run: bool) -> None:
    """Copy or symlink one file (receptor structure) to dest."""
    if dry_run or mode == "none":
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    if mode == "symlink":
        if dest.exists() or dest.is_symlink():
            dest.unlink()
        dest.symlink_to(src.resolve())
    else:
        shutil.copy2(src, dest)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--consensus_all", required=True, type=Path,
                        help="consensus_all.parquet produced by the consensus script.")
    parser.add_argument("--input", action="append", required=True, type=parse_named_path,
                        metavar="RECEPTOR=ENRICHMENT_DIR_OR_PARQUET",
                        help="Repeat once per receptor (the same inputs given to the consensus).")
    parser.add_argument("--out_dir", required=True, type=Path)
    parser.add_argument("--top_fraction", type=float, default=0.01,
                        help="Top fraction of the consensus ranking to collect (default 0.01 = 1%%).")
    parser.add_argument("--ligand_column", default="Ligand")
    parser.add_argument("--score_column", default=None,
                        help="Raw docking score column. Auto-detected (Grid_Score/Vina_Score) by default.")
    parser.add_argument("--higher_is_better", action="store_true",
                        help="Set if larger raw docking scores are better. Default: lower is better.")
    parser.add_argument("--best_by", choices=("score", "rank"), default="score",
                        help="Define best pose by lowest raw score (default) or best within-receptor rank.")
    parser.add_argument("--file_column", default=None,
                        help="Column holding the pose file path. Auto-detected by default.")
    parser.add_argument("--mol2_line_base", type=int, choices=(0, 1), default=1,
                        help="Are Mol2_Start/End_Line 1-based (default) or 0-based?")
    parser.add_argument("--mol2_header", choices=("reattach", "strip", "raw"), default="reattach",
                        help="DOCK6 mol2 records carry a '##########' metadata block ABOVE "
                             "@<TRIPOS>MOLECULE. 'reattach' (default) prepends this molecule's own "
                             "header and drops the next molecule's trailing block; 'strip' only drops "
                             "the trailing block; 'raw' keeps the exact Start/End slice.")
    parser.add_argument("--copy_whole_source", action="store_true",
                        help="Copy the entire Source_File instead of extracting the pose line range.")
    parser.add_argument("--pose_selection", choices=("best_model", "whole_file"), default="best_model",
                        help="For whole-file pose sources (vina multi-model pdbqt): 'best_model' "
                             "(default) writes only the single best-scoring MODEL; 'whole_file' copies "
                             "every pose. Files with no MODEL records are always copied whole.")
    parser.add_argument("--export_sdf", action="store_true",
                        help="After writing each pdbqt pose, run Meeko's mk_export.py to also write an "
                             "SDF of the ligand (ligand only; receptor not exported). Vina/pdbqt poses only.")
    parser.add_argument("--mk_export_cmd", default="mk_export.py",
                        help="Executable for the SDF export (default: mk_export.py on PATH).")
    parser.add_argument("--receptor_template", default=None,
                        help="Path template for the receptor structure of the best receptor, with a "
                             "'{receptor}' placeholder. Example (vina): "
                             "'/.../minimized_receptors/{receptor}.pdb'. The template should "
                             "normally resolve to a PDB when the collection will be staged for CTMD. "
                             "When given, each hit becomes a receptor+ligand pair.")
    parser.add_argument("--layout", choices=("pairs", "flat"), default=None,
                        help="'pairs' = one folder per hit holding the ligand pose and its receptor "
                             "(default when --receptor_template is set); 'flat' = all poses in poses/ "
                             "(default otherwise).")
    parser.add_argument("--receptor_mode", choices=("copy", "symlink", "none"), default="copy",
                        help="How to place the receptor next to each pose: copy (default, portable), "
                             "symlink (no duplication; stays on this filesystem), or none (manifest only).")
    parser.add_argument("--dry_run", action="store_true",
                        help="Select and write the manifest, but do not copy/extract any pose files.")
    args = parser.parse_args()

    if not 0 < args.top_fraction <= 1:
        parser.error("--top_fraction must be in (0, 1].")
    names = [n for n, _ in args.input]
    if len(names) != len(set(names)):
        parser.error("Each receptor name in --input must be unique.")

    lower_is_better = not args.higher_is_better

    # 1. top fraction of the fused ranking
    top, total, cutoff = load_top_set(args.consensus_all, args.top_fraction, args.ligand_column)
    # An unlabeled screen has no Class column. Reporting "? active, ? decoy" would
    # read as "we could not determine the labels" when the truth is that there are
    # none to determine.
    labels_present = "Class" in top.columns
    if labels_present:
        composition = (f"{(top['Class'] == 'active').sum()} active, "
                       f"{(top['Class'] == 'decoy').sum()} decoy")
    else:
        composition = "unlabeled run, no active/decoy split"
    print(f"[i] consensus ligands: {total} | top {args.top_fraction:g} -> {cutoff} ligands "
          f"({composition})")

    # 2. traceback tables + best pose per ligand
    combined, score_column = load_receptor_tables(args.input, args.ligand_column, args.score_column)
    print(f"[i] score column: {score_column} ({'lower' if lower_is_better else 'higher'} is better) "
          f"| best pose by: {args.best_by}")

    file_column = args.file_column or first_present(combined, FILE_COLUMN_CANDIDATES)
    start_col = first_present(combined, START_LINE_CANDIDATES)
    end_col = first_present(combined, END_LINE_CANDIDATES)
    extract_mode = (not args.copy_whole_source) and bool(start_col and end_col)
    if not file_column:
        print(f"[!] No pose-file column found (looked for {FILE_COLUMN_CANDIDATES}). "
              f"Manifest will still be written.", file=sys.stderr)
    print(f"[i] pose files via: {file_column} | "
          f"{'extract line range' if extract_mode else 'copy whole file'}")

    top_ligs = set(top[args.ligand_column])
    pool = combined[combined[args.ligand_column].isin(top_ligs)]
    by_ligand = pool.groupby(args.ligand_column, sort=False)

    pairing = args.receptor_template is not None
    layout = args.layout or ("pairs" if pairing else "flat")
    if pairing:
        print(f"[i] pairing receptors via: {args.receptor_template} "
              f"({args.receptor_mode}) | layout: {layout}")

    poses_dir = args.out_dir / "poses"
    records = []
    mol2_jobs: dict[Path, list[dict]] = {}
    missing_sources = 0
    missing_receptors = 0
    sdf_ok = 0
    sdf_fail = 0
    written_status = {"BEST_MODEL", "COPIED", "COPIED_WHOLE"}

    for _, trow in top.sort_values("Consensus_Rank").iterrows():
        lig = trow[args.ligand_column]
        # Carried only when the consensus was labelled; a column of empty strings
        # would misreport an unlabeled run as one whose labels went missing.
        label = {"Class": trow["Class"]} if labels_present else {}
        if lig not in by_ligand.groups:
            records.append({"Consensus_Rank": trow.get("Consensus_Rank"), "Ligand": lig,
                            **label, "Status": "NO_RECEPTOR_ROW"})
            continue
        best = pick_best_pose(by_ligand.get_group(lig), args.best_by, lower_is_better)
        receptor = best["__receptor__"]
        rank_tag = int(trow["Consensus_Rank"]) if pd.notna(trow.get("Consensus_Rank")) else 0
        src = Path(str(best[file_column])) if file_column and pd.notna(best.get(file_column)) else None

        rec = {
            "Consensus_Rank": rank_tag,
            "Ligand": lig,
            "Raw_Ligand": best.get("Raw_Ligand", lig),
            **label,
            "Consensus_RRF_Score": trow.get("Consensus_RRF_Score"),
            "Best_Receptor": receptor,
            "Best_Raw_Score": best["__score__"],
            "Within_Receptor_Rank": best.get("__rank__"),
            "Top_Support_Votes": trow.get("Top_Support_Votes"),
            "Source_File": str(src) if src else "",
        }
        if start_col and end_col and pd.notna(best.get(start_col)) and pd.notna(best.get(end_col)):
            rec["Start_Line"] = int(best[start_col])
            rec["End_Line"] = int(best[end_col])

        # destination layout: one folder per hit (pairs) or a flat poses/ dir
        ext = src.suffix if (src and src.suffix) else (".mol2" if extract_mode else ".pose")
        # Labelled runs keep the active/decoy segment they have always had. An
        # unlabeled run drops it rather than stamping every pose "na", which
        # would read as a lookup that failed.
        class_segment = f"{trow['Class'] or 'na'}_" if labels_present else ""
        stem = f"{rank_tag:05d}_{class_segment}{safe_name(lig)}_{safe_name(receptor)}"
        if layout == "pairs":
            hit_dir = poses_dir / stem
            dest = hit_dir / f"{safe_name(lig)}_pose{ext}"
        else:
            hit_dir = poses_dir
            dest = poses_dir / f"{stem}{ext}"
        rec["Collected_File"] = str(dest)
        if not args.dry_run:
            dest.parent.mkdir(parents=True, exist_ok=True)

        # place the ligand pose (extract line range for mol2, else copy whole file)
        if src is None:
            rec["Status"] = "NO_SOURCE_PATH"
        elif extract_mode and "Start_Line" in rec:
            if not src.exists():
                rec["Status"] = "SOURCE_MISSING"
                missing_sources += 1
            else:
                mol2_jobs.setdefault(src, []).append(
                    {"start": rec["Start_Line"], "end": rec["End_Line"], "dest": dest})
                rec["Status"] = "EXTRACT"
        elif args.pose_selection == "best_model":  # vina multi-model pdbqt -> best pose only
            if not src.exists():
                rec["Status"] = "SOURCE_MISSING"
                missing_sources += 1
            elif args.dry_run:
                rec["Status"] = "BEST_MODEL"
            else:
                all_lines = src.read_text(errors="replace").splitlines(keepends=True)
                best = extract_best_vina_model(all_lines, lower_is_better)
                if best is None:  # no MODEL records; nothing to trim
                    shutil.copy2(src, dest)
                    rec["Status"] = "COPIED_WHOLE"
                else:
                    affinity, block = best
                    text = "".join(block)
                    dest.write_text(text if text.endswith("\n") else text + "\n")
                    rec["Model_Affinity"] = affinity
                    rec["Status"] = "BEST_MODEL"
        else:  # whole-file copy (every pose)
            if not src.exists():
                rec["Status"] = "SOURCE_MISSING"
                missing_sources += 1
            elif not args.dry_run:
                shutil.copy2(src, dest)
                rec["Status"] = "COPIED"
            else:
                rec["Status"] = "COPY"

        # place the paired receptor structure
        if pairing:
            recv_src = Path(args.receptor_template.format(receptor=receptor))
            recv_dest = hit_dir / recv_src.name
            rec["Receptor_File"] = str(recv_src)
            rec["Collected_Receptor"] = str(recv_dest) if args.receptor_mode != "none" else ""
            if not recv_src.exists():
                rec["Receptor_Status"] = "RECEPTOR_MISSING"
                missing_receptors += 1
            else:
                place_file(recv_src, recv_dest, args.receptor_mode, args.dry_run)
                rec["Receptor_Status"] = args.receptor_mode.upper() if not args.dry_run else "PENDING"

        # optional: export the ligand pose to SDF via Meeko (pdbqt poses only)
        if args.export_sdf and dest.suffix.lower() == ".pdbqt":
            sdf = dest.with_suffix(".sdf")
            rec["Exported_SDF"] = str(sdf)
            if args.dry_run:
                rec["SDF_Status"] = "PENDING"
            elif rec.get("Status") in written_status:
                ok, msg = run_mk_export(args.mk_export_cmd, dest, sdf)
                rec["SDF_Status"] = "OK" if ok else f"FAILED: {msg}"
                sdf_ok += int(ok)
                sdf_fail += int(not ok)
            else:
                rec["SDF_Status"] = "SKIPPED_NO_POSE"

        records.append(rec)

    # run the grouped mol2 extractions (one streaming pass per shared source file)
    for source, jobs in mol2_jobs.items():
        extract_mol2_ranges(source, jobs, args.mol2_line_base, args.mol2_header, args.dry_run)

    manifest = pd.DataFrame(records)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.out_dir / "collection_manifest.csv"
    manifest.to_csv(manifest_path, index=False)

    # report
    status_counts = manifest["Status"].value_counts().to_dict() if "Status" in manifest else {}
    print(f"[i] selected ligands: {len(manifest)}")
    print(f"[i] status breakdown: {status_counts}")
    if "Receptor_Status" in manifest:
        print(f"[i] receptor status: {manifest['Receptor_Status'].value_counts().to_dict()}")
    if args.export_sdf:
        print(f"[i] SDF export: {sdf_ok} written, {sdf_fail} failed"
              + (" (dry run)" if args.dry_run else ""))
    if missing_sources:
        print(f"[!] {missing_sources} best-pose source files were not found on this filesystem "
              f"(expected if you are running where the docking outputs do not live).",
              file=sys.stderr)
    if missing_receptors:
        print(f"[!] {missing_receptors} receptor files were not found via --receptor_template.",
              file=sys.stderr)
    if args.dry_run:
        print(f"[i] DRY RUN: no pose files written. Manifest: {manifest_path}")
        print("[i] Manifest preview (first 10 rows):")
        cols = [c for c in ("Consensus_Rank", "Ligand", "Class", "Best_Receptor",
                            "Best_Raw_Score", "Status", "Receptor_Status") if c in manifest]
        with pd.option_context("display.max_colwidth", 60, "display.width", 200):
            print(manifest[cols].head(10).to_string(index=False))
    else:
        print(f"[+] poses dir: {poses_dir}")
        print(f"[+] manifest:  {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
