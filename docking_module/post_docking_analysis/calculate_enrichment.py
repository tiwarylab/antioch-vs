#!/usr/bin/env python3
"""Calculate DUD-E-compatible enrichment metrics from docking summary parquets."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd

DEFAULT_LOGAUC_MIN_DECOY_FRACTION = 0.001
DEFAULT_EF_DECOY_FRACTION = 0.01
SCORE_COLUMN_CANDIDATES = ["Grid_Score", "Vina_Score", "Total_Score", "Score", "score"]
LIGAND_COLUMN_CANDIDATES = ["Ligand", "Name", "ligand", "name"]
AGAIN_SUFFIX_REGEX = r"-again\d+$"


def discover_parquets(path: Path) -> list[Path]:
    if path.is_file():
        if path.suffix != ".parquet":
            raise ValueError(f"Input file is not a parquet file: {path}")
        return [path]

    if not path.is_dir():
        raise FileNotFoundError(f"Input path does not exist or is not a directory: {path}")

    top_level = sorted(path.glob("*.parquet"))
    if top_level:
        return top_level

    return sorted(path.rglob("*.parquet"))


def read_labeled_parquets(path: Path, label: str) -> pd.DataFrame:
    parquet_paths = discover_parquets(path)
    if not parquet_paths:
        raise FileNotFoundError(f"No parquet files found in {path}")

    frames = []
    for parquet_path in parquet_paths:
        frame = pd.read_parquet(parquet_path)
        frame["Source_Parquet"] = str(parquet_path.resolve())
        frames.append(frame)

    df = pd.concat(frames, ignore_index=True)
    df["Class"] = label
    df["Is_Active"] = label == "active"
    return df


def choose_column(df: pd.DataFrame, requested: str | None, candidates: list[str], kind: str) -> str:
    if requested:
        if requested not in df.columns:
            raise ValueError(f"Requested {kind} column {requested!r} was not found.")
        return requested

    for column in candidates:
        if column in df.columns:
            return column

    raise ValueError(
        f"Could not auto-detect a {kind} column. Available columns: {', '.join(df.columns)}"
    )


def parse_csv_floats(value: str | None) -> list[float]:
    if not value:
        return []
    fractions = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        fraction = float(item)
        if fraction <= 0:
            raise ValueError("Enrichment fractions must be positive.")
        fractions.append(fraction)
    return fractions


def parse_csv_ints(value: str | None) -> list[int]:
    if not value:
        return []
    counts = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        count = int(item)
        if count <= 0:
            raise ValueError("Top-count cutoffs must be positive integers.")
        counts.append(count)
    return counts


def rank_table(
    df: pd.DataFrame,
    score_column: str,
    ligand_column: str,
    lower_is_better: bool,
    best_per_ligand: bool,
) -> pd.DataFrame:
    ranked = df.copy()
    ranked[score_column] = pd.to_numeric(ranked[score_column], errors="coerce")
    ranked = ranked.dropna(subset=[score_column])

    if ligand_column not in ranked.columns:
        ranked[ligand_column] = ranked.index.astype(str)

    ranked = ranked.sort_values(score_column, ascending=lower_is_better, na_position="last")

    if best_per_ligand:
        dedupe_columns = ["Class", ligand_column] if "Class" in ranked.columns else [ligand_column]
        ranked = ranked.drop_duplicates(subset=dedupe_columns, keep="first")

    ranked = ranked.reset_index(drop=True)
    ranked["Rank"] = ranked.index + 1
    ranked["Cumulative_Actives"] = ranked["Is_Active"].cumsum()
    ranked["Cumulative_Decoys"] = ranked["Rank"] - ranked["Cumulative_Actives"]
    return ranked


def strip_again_suffix(df: pd.DataFrame, ligand_column: str) -> pd.DataFrame:
    """Collapse ANTIOCH-VS/Meeko duplicate filenames such as ``CHEMBL1-again2``."""
    normalized = df.copy()
    normalized["Raw_Ligand"] = normalized[ligand_column]
    normalized[ligand_column] = normalized[ligand_column].astype(str).str.replace(
        AGAIN_SUFFIX_REGEX,
        "",
        regex=True,
    )
    return normalized


def build_roc_curve(ranked: pd.DataFrame, score_column: str) -> pd.DataFrame:
    """Build a score-tie-aware ROC curve with decoy recall on the x-axis."""
    active_count = int(ranked["Is_Active"].sum())
    decoy_count = int((~ranked["Is_Active"]).sum())
    if active_count == 0 or decoy_count == 0:
        return pd.DataFrame()

    rows = [{"Decoy_Recall": 0.0, "Active_Recall": 0.0}]
    active_seen = 0
    decoy_seen = 0

    # Rows are already sorted by score. Processing equal scores together makes
    # AUC, LogAUC, and EF independent of arbitrary ordering within a score tie.
    for _, score_group in ranked.groupby(score_column, sort=False, dropna=False):
        active_seen += int(score_group["Is_Active"].sum())
        decoy_seen += int((~score_group["Is_Active"]).sum())
        rows.append(
            {
                "Decoy_Recall": decoy_seen / decoy_count,
                "Active_Recall": active_seen / active_count,
            }
        )

    curve = pd.DataFrame(rows)
    # A pure-active score group creates a vertical ROC segment. For interpolation
    # at a positive decoy fraction, retain the highest recall at that x value.
    return curve.groupby("Decoy_Recall", as_index=False)["Active_Recall"].max()


def interpolate_active_recall(roc_curve: pd.DataFrame, decoy_fraction: float) -> float | None:
    """Return ROC active recall at a decoy-recall cutoff by linear interpolation."""
    if roc_curve.empty:
        return None
    if not 0.0 <= decoy_fraction <= 1.0:
        raise ValueError("Decoy fraction must be between 0 and 1.")

    x_values = roc_curve["Decoy_Recall"].tolist()
    y_values = roc_curve["Active_Recall"].tolist()
    if decoy_fraction <= x_values[0]:
        return float(y_values[0])

    for index in range(1, len(x_values)):
        x_left, x_right = x_values[index - 1], x_values[index]
        y_left, y_right = y_values[index - 1], y_values[index]
        if decoy_fraction <= x_right:
            if math.isclose(x_left, x_right):
                return float(max(y_left, y_right))
            fraction = (decoy_fraction - x_left) / (x_right - x_left)
            return float(y_left + fraction * (y_right - y_left))
    return float(y_values[-1])


def calculate_auc(roc_curve: pd.DataFrame) -> float | None:
    if roc_curve.empty:
        return None

    area = 0.0
    for left, right in zip(roc_curve.itertuples(index=False), roc_curve.iloc[1:].itertuples(index=False)):
        area += 0.5 * (left.Active_Recall + right.Active_Recall) * (
            right.Decoy_Recall - left.Decoy_Recall
        )
    return area


def calculate_log_auc(
    roc_curve: pd.DataFrame,
    min_decoy_fraction: float = DEFAULT_LOGAUC_MIN_DECOY_FRACTION,
) -> tuple[float | None, float | None, pd.DataFrame]:
    if roc_curve.empty:
        return None, None, pd.DataFrame()
    if not 0 < min_decoy_fraction < 1:
        raise ValueError("--logauc_min_decoy_fraction must be between 0 and 1.")

    active_at_min = interpolate_active_recall(roc_curve, min_decoy_fraction)
    points = [(min_decoy_fraction, active_at_min)]
    points.extend(
        (float(row.Decoy_Recall), float(row.Active_Recall))
        for row in roc_curve.itertuples(index=False)
        if row.Decoy_Recall > min_decoy_fraction
    )

    log_min = math.log10(min_decoy_fraction)
    log_max = 0.0
    denominator = log_max - log_min
    area = 0.0

    for (x1, y1), (x2, y2) in zip(points, points[1:]):
        log_x1 = math.log10(x1)
        log_x2 = math.log10(x2)
        area += 0.5 * (y1 + y2) * (log_x2 - log_x1)

    raw_log_auc_fraction = area / denominator
    random_log_auc_fraction = (1.0 - min_decoy_fraction) / (math.log(10.0) * denominator)
    adjusted_log_auc_fraction = raw_log_auc_fraction - random_log_auc_fraction

    curve = pd.DataFrame(
        {
            "Decoy_Recall": [point[0] for point in points],
            "Active_Recall": [point[1] for point in points],
            "Log10_Decoy_Recall": [math.log10(point[0]) for point in points],
        }
    )
    return (
        adjusted_log_auc_fraction * 100.0,
        random_log_auc_fraction * 100.0,
        curve,
    )


def diagnostic_row(ranked: pd.DataFrame, cutoff_name: str, cutoff_count: int) -> dict[str, float | int | str]:
    total_count = len(ranked)
    active_count = int(ranked["Is_Active"].sum())
    decoy_count = total_count - active_count
    cutoff_count = min(cutoff_count, total_count)
    top = ranked.head(cutoff_count)
    actives_found = int(top["Is_Active"].sum())
    decoys_found = cutoff_count - actives_found

    active_fraction = active_count / total_count if total_count else 0.0
    expected_actives = cutoff_count * active_fraction
    enrichment_factor = actives_found / expected_actives if expected_actives else None
    active_recall = actives_found / active_count if active_count else None
    decoy_recall = decoys_found / decoy_count if decoy_count else None

    return {
        "Diagnostic": cutoff_name,
        "Ranked_Ligands": cutoff_count,
        "Ranked_Ligand_Fraction": cutoff_count / total_count if total_count else None,
        "Actives_Recovered": actives_found,
        "Decoys_Recovered": decoys_found,
        "Active_Recall_Percent": active_recall * 100.0 if active_recall is not None else None,
        "Decoy_Recall_Percent": decoy_recall * 100.0 if decoy_recall is not None else None,
        "Conventional_Fold_Enrichment": enrichment_factor,
    }


def calculate_metrics(
    ranked: pd.DataFrame,
    fractions: list[float],
    top_counts: list[int],
    score_column: str,
    lower_is_better: bool,
    ef_decoy_fraction: float,
    logauc_min_decoy_fraction: float,
    input_summary: dict[str, object] | None = None,
) -> tuple[pd.DataFrame, dict[str, object], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    total_count = len(ranked)
    active_count = int(ranked["Is_Active"].sum())
    decoy_count = total_count - active_count

    diagnostic_rows = []
    for fraction in fractions:
        cutoff_count = max(1, math.ceil(total_count * fraction))
        diagnostic_rows.append(diagnostic_row(ranked, f"top_{fraction:g}_of_all_ligands", cutoff_count))

    for count in top_counts:
        diagnostic_rows.append(diagnostic_row(ranked, f"top_{count}_ligands", count))

    roc_curve = build_roc_curve(ranked, score_column)
    auc = calculate_auc(roc_curve)
    dude_log_auc, _random_log_auc, log_auc_curve = calculate_log_auc(
        roc_curve,
        min_decoy_fraction=logauc_min_decoy_fraction,
    )
    ef_active_recall = interpolate_active_recall(roc_curve, ef_decoy_fraction)

    metrics = pd.DataFrame(
        [
            {
                "Metric": "ROC_AUC_Percent",
                "Value": auc * 100.0 if auc is not None else None,
                "Definition": "Area under the ROC curve: active recall versus decoy recall across the full ranking.",
            },
            {
                "Metric": f"DUD_E_LogAUC_{logauc_min_decoy_fraction:g}_Percent",
                "Value": dude_log_auc,
                "Definition": (
                    "DUD-E adjusted LogAUC: semilog ROC area from the stated minimum decoy recall "
                    "to 100%, minus the random baseline."
                ),
            },
            {
                "Metric": f"EF{ef_decoy_fraction * 100:g}_Percent",
                "Value": ef_active_recall * 100.0 if ef_active_recall is not None else None,
                "Definition": "Percent of actives recovered when the stated percent of decoys has been recovered.",
            },
        ]
    )
    ef_metric_key = f"ef{ef_decoy_fraction * 100:g}_percent"
    summary = {
        "total_ligands": total_count,
        "active_ligands": active_count,
        "decoy_ligands": decoy_count,
        "active_fraction": active_count / total_count if total_count else None,
        "score_column": score_column,
        "score_direction": "lower_is_better" if lower_is_better else "higher_is_better",
        "roc_auc_percent": auc * 100.0 if auc is not None else None,
        "dude_logauc_percent": dude_log_auc,
        ef_metric_key: ef_active_recall * 100.0 if ef_active_recall is not None else None,
        "ef_decoy_recall_percent": ef_decoy_fraction * 100.0,
        "logauc_min_decoy_fraction": logauc_min_decoy_fraction,
    }
    if input_summary:
        summary.update(input_summary)
    return metrics, summary, log_auc_curve, pd.DataFrame(diagnostic_rows), roc_curve


def summarize_inputs(
    active_df: pd.DataFrame,
    decoy_df: pd.DataFrame,
    combined_df: pd.DataFrame,
    score_column: str,
    ligand_column: str,
) -> dict[str, object]:
    active_scored = active_df[pd.to_numeric(active_df[score_column], errors="coerce").notna()]
    decoy_scored = decoy_df[pd.to_numeric(decoy_df[score_column], errors="coerce").notna()]
    combined_scored = combined_df[pd.to_numeric(combined_df[score_column], errors="coerce").notna()]

    return {
        "input_active_rows": len(active_df),
        "input_decoy_rows": len(decoy_df),
        "input_active_unique_ligands": int(active_df[ligand_column].nunique(dropna=True)),
        "input_decoy_unique_ligands": int(decoy_df[ligand_column].nunique(dropna=True)),
        "scored_active_rows": len(active_scored),
        "scored_decoy_rows": len(decoy_scored),
        "scored_active_unique_ligands": int(active_scored[ligand_column].nunique(dropna=True)),
        "scored_decoy_unique_ligands": int(decoy_scored[ligand_column].nunique(dropna=True)),
        "combined_scored_rows": len(combined_scored),
    }


def write_outputs(
    ranked: pd.DataFrame,
    metrics: pd.DataFrame,
    diagnostics: pd.DataFrame,
    roc_curve: pd.DataFrame,
    log_auc_curve: pd.DataFrame,
    summary: dict[str, object],
    out_dir: Path,
    out_prefix: str,
    csv_n: int | None,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    ranked_parquet = out_dir / f"{out_prefix}_ranked.parquet"
    ranked_csv = out_dir / f"{out_prefix}_ranked.csv"
    metrics_csv = out_dir / f"{out_prefix}_metrics.csv"
    diagnostics_csv = out_dir / f"{out_prefix}_diagnostics.csv"
    roc_curve_csv = out_dir / f"{out_prefix}_roc_curve.csv"
    log_auc_curve_csv = out_dir / f"{out_prefix}_logauc_curve.csv"
    summary_json = out_dir / f"{out_prefix}_summary.json"

    ranked.to_parquet(ranked_parquet, index=False)
    ranked.head(csv_n).to_csv(ranked_csv, index=False)
    metrics.to_csv(metrics_csv, index=False)
    if not diagnostics.empty:
        diagnostics.to_csv(diagnostics_csv, index=False)
    roc_curve.to_csv(roc_curve_csv, index=False)
    if not log_auc_curve.empty:
        log_auc_curve.to_csv(log_auc_curve_csv, index=False)
    summary_json.write_text(json.dumps(summary, indent=2) + "\n")

    print(f"[+] Ranked parquet: {ranked_parquet}")
    print(f"[+] Ranked CSV: {ranked_csv}")
    print(f"[+] Enrichment metrics: {metrics_csv}")
    if not diagnostics.empty:
        print(f"[+] Optional cutoff diagnostics: {diagnostics_csv}")
    print(f"[+] ROC curve: {roc_curve_csv}")
    if not log_auc_curve.empty:
        print(f"[+] LogAUC curve: {log_auc_curve_csv}")
    print(f"[+] Summary JSON: {summary_json}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Calculate DUD-E-compatible enrichment metrics from active and decoy docking summaries."
    )
    parser.add_argument("--actives_dir", required=True, help="Folder or parquet file containing active summaries.")
    parser.add_argument("--decoys_dir", required=True, help="Folder or parquet file containing decoy summaries.")
    parser.add_argument("--out_dir", default="enrichment_summary", help="Output directory.")
    parser.add_argument("--out_prefix", default="enrichment", help="Prefix for output files.")
    parser.add_argument("--score_column", default=None, help="Score column. Auto-detected by default.")
    parser.add_argument("--ligand_column", default=None, help="Ligand/name column. Auto-detected by default.")
    parser.add_argument(
        "--higher_is_better",
        action="store_true",
        help="Use this if larger scores are better. Default assumes lower docking scores are better.",
    )
    parser.add_argument(
        "--keep_all_poses",
        action="store_true",
        help="Do not collapse to the best pose per ligand before enrichment calculation.",
    )
    parser.add_argument(
        "--strip_again_suffix",
        action="store_true",
        help=(
            "Treat names ending in -againN as duplicate preparations of the same ligand, "
            "then retain the best score per parent ligand. For example, CHEMBL1-again2 "
            "is collapsed into CHEMBL1."
        ),
    )
    parser.add_argument(
        "--diagnostic_fractions",
        "--fractions",
        dest="diagnostic_fractions",
        default=None,
        help=(
            "Optional comma-separated top fractions of all ranked ligands for diagnostic reporting. "
            "These are conventional fold-enrichment diagnostics, not DUD-E EF1."
        ),
    )
    parser.add_argument(
        "--diagnostic_top_counts",
        "--top_counts",
        dest="diagnostic_top_counts",
        default=None,
        help="Optional absolute top-N cutoffs for diagnostic reporting, for example 100,500,1000.",
    )
    parser.add_argument(
        "--ef_decoy_fraction",
        type=float,
        default=DEFAULT_EF_DECOY_FRACTION,
        help=(
            "Decoy recall at which to report active recovery. The DUD-E EF1 convention is 0.01 "
            "(1%% of decoys); this is the default."
        ),
    )
    parser.add_argument(
        "--roc_ef_decoy_fractions",
        dest="legacy_roc_ef_decoy_fractions",
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--csv_n",
        type=int,
        default=None,
        help="Number of ranked rows to write to CSV. Default writes all rows.",
    )
    parser.add_argument(
        "--logauc_min_decoy_fraction",
        type=float,
        default=DEFAULT_LOGAUC_MIN_DECOY_FRACTION,
        help=(
            "Lower decoy-fraction bound for DUD-E-style semilog ROC integration. "
            f"Default: {DEFAULT_LOGAUC_MIN_DECOY_FRACTION}."
        ),
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    active_df = read_labeled_parquets(Path(args.actives_dir), "active")
    decoy_df = read_labeled_parquets(Path(args.decoys_dir), "decoy")
    combined = pd.concat([active_df, decoy_df], ignore_index=True, sort=False)

    score_column = choose_column(combined, args.score_column, SCORE_COLUMN_CANDIDATES, "score")
    ligand_column = choose_column(combined, args.ligand_column, LIGAND_COLUMN_CANDIDATES, "ligand")
    if args.strip_again_suffix:
        combined = strip_again_suffix(combined, ligand_column)
    lower_is_better = not args.higher_is_better
    ranked = rank_table(
        combined,
        score_column=score_column,
        ligand_column=ligand_column,
        lower_is_better=lower_is_better,
        best_per_ligand=not args.keep_all_poses,
    )

    if ranked.empty:
        raise RuntimeError("No scored rows were available after filtering missing scores.")
    if args.legacy_roc_ef_decoy_fractions:
        legacy_fractions = parse_csv_floats(args.legacy_roc_ef_decoy_fractions)
        if len(legacy_fractions) != 1:
            raise ValueError(
                "--roc_ef_decoy_fractions is deprecated and now accepts exactly one value. "
                "Use --ef_decoy_fraction instead."
            )
        args.ef_decoy_fraction = legacy_fractions[0]
    if not 0 < args.ef_decoy_fraction <= 1:
        raise ValueError("--ef_decoy_fraction must be between 0 and 1.")

    metrics, summary, log_auc_curve, diagnostics, roc_curve = calculate_metrics(
        ranked,
        fractions=parse_csv_floats(args.diagnostic_fractions),
        top_counts=parse_csv_ints(args.diagnostic_top_counts),
        score_column=score_column,
        lower_is_better=lower_is_better,
        ef_decoy_fraction=args.ef_decoy_fraction,
        logauc_min_decoy_fraction=args.logauc_min_decoy_fraction,
        input_summary=summarize_inputs(
            active_df=active_df,
            decoy_df=decoy_df,
            combined_df=combined,
            score_column=score_column,
            ligand_column=ligand_column,
        ),
    )
    summary["ligand_name_normalization"] = (
        "strip_again_suffix" if args.strip_again_suffix else "none"
    )

    print("\nEnrichment summary")
    print("==================")
    for key, value in summary.items():
        print(f"{key}: {value}")
    print("")
    print(metrics.to_string(index=False))

    write_outputs(
        ranked=ranked,
        metrics=metrics,
        diagnostics=diagnostics,
        roc_curve=roc_curve,
        log_auc_curve=log_auc_curve,
        summary=summary,
        out_dir=Path(args.out_dir),
        out_prefix=args.out_prefix,
        csv_n=args.csv_n,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
