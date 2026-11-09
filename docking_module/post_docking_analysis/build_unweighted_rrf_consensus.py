#!/usr/bin/env python3
"""Build an unweighted multi-receptor consensus ranking with reciprocal rank fusion.

Inputs are one per-receptor ranking table per receptor, supplied as
``--input NAME=PATH``. ``PATH`` may be a Parquet file directly or a directory
containing ``enrichment_ranked.parquet``.

Two kinds of run are supported and they differ only in what is passed in:

* **Unlabeled screening** — hand it the extractor output directly
  (``docking_summary.parquet`` from Vina, ``dock6_summary.parquet`` from DOCK6).
  There are no active/decoy labels and none are needed: fusion depends on rank
  order, not on labels. The label-derived evaluation inputs are simply not
  written, and ``manifest.json`` records ``labels_present: false``.
* **Benchmark / DUD-E** — hand it ``enrichment_ranked.parquet`` from
  ``calculate_enrichment.py``, as before. Behaviour on this path is unchanged.

Labels are never used to weight or construct the consensus ranking.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

SCORE_CANDIDATES = ("Vina_Score", "Grid_Score", "Total_Score", "Score", "score")

VALID_CLASSES = frozenset({"active", "decoy"})

# Everything this stage writes. The guard below refuses a directory already
# holding any of them, and --overwrite clears all of them — including the ones
# this particular run will not produce. Removing only what is about to be
# rewritten is the defect: an unlabeled rerun would leave a labelled
# consensus_actives.parquet next to a manifest reading labels_present: false.
OWNED_ARTIFACTS = (
    "consensus_all.parquet",
    "consensus_top_hits.csv",
    "consensus_actives.parquet",
    "consensus_decoys.parquet",
    "consensus_manifest.csv",
    "manifest.json",
    "failures.csv",
)

# The engine is not recorded in the extractor tables, so the score column is the
# only evidence available. Where it says nothing, the manifest says null rather
# than guessing.
ENGINE_BY_SCORE_COLUMN = {"Vina_Score": "vina", "Grid_Score": "dock6"}


def parse_named_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("Each --input must use NAME=PATH.")
    name, path_text = value.split("=", 1)
    name = name.strip()
    path = Path(path_text.strip()).expanduser()
    if not name or not path_text.strip():
        raise argparse.ArgumentTypeError("Both NAME and PATH are required.")
    return name, path


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", value).strip("_") or "receptor"


def resolve_parquet(path: Path) -> Path:
    if path.is_file():
        if path.suffix != ".parquet":
            raise ValueError(f"Input is not a Parquet file: {path}")
        return path
    candidate = path / "enrichment_ranked.parquet"
    if candidate.is_file():
        return candidate
    raise FileNotFoundError(
        f"Expected an enrichment_ranked.parquet directory or Parquet file, got: {path}"
    )


def choose_score_column(frame: pd.DataFrame, requested: str | None) -> str:
    if requested:
        if requested not in frame.columns:
            raise ValueError(
                f"Score column {requested!r} was not found: {', '.join(frame.columns)}"
            )
        return requested
    for column in SCORE_CANDIDATES:
        if column in frame.columns:
            return column
    raise ValueError(f"Could not find a docking score column: {', '.join(frame.columns)}")


def read_ranking(
    name: str,
    path: Path,
    ligand_column: str,
    score_column: str | None,
    lower_is_better: bool,
) -> tuple[pd.DataFrame, str]:
    parquet = resolve_parquet(path)
    frame = pd.read_parquet(parquet)
    if ligand_column not in frame.columns:
        raise ValueError(f"{name}: missing required column(s): {ligand_column}")

    selected_score = choose_score_column(frame, score_column)

    # Class is optional. Benchmark tables carry it; extractor tables from an
    # unlabeled screen do not, and that is a normal run rather than an error.
    has_class = "Class" in frame.columns
    columns = [ligand_column, "Class", selected_score] if has_class else [ligand_column, selected_score]
    frame = frame[columns].copy()
    frame[selected_score] = pd.to_numeric(frame[selected_score], errors="coerce")
    frame = frame.dropna(subset=[ligand_column, selected_score])
    frame = frame.sort_values(selected_score, ascending=lower_is_better, kind="stable")

    if has_class:
        present = set(frame["Class"].dropna().unique())

        unrecognised = present - VALID_CLASSES
        if unrecognised:
            raise ValueError(
                f"{name}: Class must be one of {sorted(VALID_CLASSES)}; found "
                f"{sorted(unrecognised)}. A value outside the two recognised classes "
                "cannot be counted as either, and guessing would silently reassign "
                "ligands between them."
            )

        # Restored after the unlabeled slice dropped it. Labels being optional
        # for a *run* does not make them optional per receptor: enrichment over a
        # single class is undefined rather than merely poor — ROC-AUC needs both
        # — so a one-class table is rejected here rather than propagated into
        # metrics that look computable.
        missing_classes = VALID_CLASSES - present
        if missing_classes:
            raise ValueError(
                f"{name}: labelled tables must contain both classes; "
                f"{sorted(missing_classes)} is absent. Enrichment over one class "
                "has no meaning."
            )

        # De-duplicating on the ligand alone would quietly collapse a ligand
        # labelled both ways into whichever row sorted first. That is a data
        # error, and before this it was caught only as a side effect of Class
        # sitting in the de-duplication key.
        class_count = frame.groupby(ligand_column)["Class"].nunique()
        if (class_count > 1).any():
            raise ValueError(f"{name}: a ligand occurs in both active and decoy classes.")

    # One row per ligand, best score first. A DOCK6 table has one row per pose,
    # so this is where several poses collapse to the ligand's best.
    frame = frame.drop_duplicates(subset=[ligand_column], keep="first")
    if frame.empty:
        raise ValueError(f"{name}: no scored ligands remained after filtering.")

    return frame.rename(columns={ligand_column: "Ligand"}), selected_score


def prepare_out_dir(out_dir: Path, overwrite: bool) -> None:
    """Make sure the run writes into a directory holding nothing of its own.

    A missing directory is created; absence is not a dirty directory. Files this
    stage does not own are never touched — an --out_dir holding unrelated work
    is not ours to clean.
    """
    existing = [name for name in OWNED_ARTIFACTS if (out_dir / name).exists()]

    # A directory sitting where an artifact belongs cannot be overwritten in any
    # meaningful sense, and writing over it fails obscurely part-way through.
    directories = [name for name in existing if (out_dir / name).is_dir()]
    if directories:
        raise ValueError(
            "Output path(s) exist as directories where this stage writes files: "
            f"{', '.join(sorted(directories))} in {out_dir}. Remove or rename them; "
            "--overwrite does not delete directory trees."
        )

    if existing and not overwrite:
        raise ValueError(
            f"{out_dir} already holds consensus outputs: {', '.join(sorted(existing))}. "
            "Refusing to write a run on top of another one — a partial overlap leaves "
            "artifacts from two runs side by side with no way to tell them apart. "
            "Pass --overwrite to replace them, or choose an empty --out_dir."
        )

    if overwrite:
        for name in OWNED_ARTIFACTS:
            (out_dir / name).unlink(missing_ok=True)


def read_expected_ligands(path: Path) -> list[str]:
    """Read the expected-ligand inventory, preserving order and dropping blanks."""
    seen: dict[str, None] = {}
    for line in path.read_text().splitlines():
        identifier = line.strip()
        if identifier:
            seen.setdefault(identifier, None)
    if not seen:
        raise ValueError(f"Expected-ligand inventory is empty: {path}")
    return list(seen)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        type=parse_named_path,
        metavar="RECEPTOR=SUMMARY_DIR_OR_PARQUET",
        help="Repeat once per receptor. Inputs must all use the same docking engine.",
    )
    parser.add_argument("--out_dir", required=True, type=Path, help="Output directory.")
    parser.add_argument("--ligand_column", default="Ligand")
    parser.add_argument("--score_column", default=None, help="Auto-detected by default.")
    parser.add_argument(
        "--expected_ligands",
        default=None,
        type=Path,
        help=(
            "Optional inventory of ligands the run intended to dock, one identifier "
            "per line. Only with this can a ligand that failed in EVERY receptor be "
            "reported: consensus otherwise never sees it. Supplying it enables "
            "failures.csv."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help=(
            "Replace consensus outputs already present in --out_dir. Every artifact "
            "this stage owns is removed first, including ones this run will not "
            "write, so a rerun cannot leave the previous run's files behind. Files "
            "this stage does not own are never touched."
        ),
    )
    parser.add_argument(
        "--higher_is_better",
        action="store_true",
        help="Use only for source score columns where larger values are favorable.",
    )
    parser.add_argument(
        "--rrf_k",
        type=float,
        default=60.0,
        help="Reciprocal-rank fusion smoothing constant (default: 60).",
    )
    parser.add_argument(
        "--support_fraction",
        type=float,
        default=0.01,
        help="Fraction used for Top_Support_Votes (default: 0.01, i.e. top 1%%).",
    )
    parser.add_argument(
        "--csv_n", type=int, default=10000, help="Rows to write to the preview CSV."
    )
    args = parser.parse_args()

    if args.rrf_k < 0:
        parser.error("--rrf_k must be non-negative.")
    if not 0 < args.support_fraction <= 1:
        parser.error("--support_fraction must be in (0, 1].")

    names = [name for name, _ in args.input]
    if len(names) != len(set(names)):
        parser.error("Each receptor name in --input must be unique.")

    expected_ligands = (
        read_expected_ligands(args.expected_ligands) if args.expected_ligands else None
    )

    lower_is_better = not args.higher_is_better

    # Fusion sums a float per receptor, and float addition is not associative:
    # the order the receptors are added in decides the last bits of the result,
    # which is enough to swap two otherwise-tied ligands. Sorting the inputs by
    # receptor name fixes that order regardless of how --input was supplied
    # (FR-014, SC-011).
    ordered_inputs = sorted(args.input, key=lambda pair: pair[0])

    rankings: dict[str, pd.DataFrame] = {}
    score_columns: dict[str, str] = {}
    labelled_receptors: list[str] = []
    unlabelled_receptors: list[str] = []
    classes_by_ligand: pd.Series | None = None
    resolved_score_column: str | None = None
    score_column_source: str | None = None

    for name, path in ordered_inputs:
        frame, score_column = read_ranking(
            name,
            path,
            args.ligand_column,
            args.score_column,
            lower_is_better,
        )

        # Fusing ranks derived from different score columns would compare a Vina
        # affinity against a DOCK6 grid score. The ranks would still compute, so
        # this has to be an error rather than a coercion (FR-015).
        if resolved_score_column is None:
            resolved_score_column, score_column_source = score_column, name
        elif score_column != resolved_score_column:
            raise ValueError(
                "Receptors resolved different score columns: "
                f"{score_column_source!r} uses {resolved_score_column!r} but "
                f"{name!r} uses {score_column!r}. Pass --score_column to state "
                "which one applies, or fuse these receptors separately."
            )

        ranking = frame[["Ligand", score_column]].copy()
        ranking["Rank"] = ranking[score_column].rank(
            method="min", ascending=lower_is_better
        )
        rankings[name] = ranking[["Ligand", "Rank"]]
        score_columns[name] = score_column

        if "Class" in frame.columns:
            labelled_receptors.append(name)
            current_classes = frame.set_index("Ligand")["Class"]
            if classes_by_ligand is None:
                classes_by_ligand = current_classes
            else:
                overlap = classes_by_ligand.index.intersection(current_classes.index)
                changed = classes_by_ligand.loc[overlap] != current_classes.loc[overlap]
                if changed.any():
                    raise ValueError(f"{name}: active/decoy class differs for one or more ligand IDs.")
                classes_by_ligand = classes_by_ligand.combine_first(current_classes)
        else:
            unlabelled_receptors.append(name)

    # Labels are a property of the run, not of each table. Neither fallback is
    # honest: inferring labels for the tables that lack them invents data, and
    # dropping them from the tables that carry them discards the benchmark.
    if labelled_receptors and unlabelled_receptors:
        raise ValueError(
            "Inputs mix labelled and unlabeled tables. With Class: "
            f"{', '.join(labelled_receptors)}. Without: "
            f"{', '.join(unlabelled_receptors)}. A run is wholly one or the other — "
            "supply labels for every receptor, or for none."
        )

    labels_present = classes_by_ligand is not None

    # The consensus covers the union of every ligand any receptor scored. Built
    # in sorted order so the frame does not inherit whichever receptor happened
    # to be read first.
    union: set[str] = set()
    for ranking in rankings.values():
        union.update(ranking["Ligand"])
    consensus = pd.DataFrame(index=pd.Index(sorted(union), name="Ligand"))

    if labels_present:
        consensus["Class"] = classes_by_ligand.reindex(consensus.index)

    # One penalty for the whole run, not one per receptor. Receptors finish
    # unevenly, so a per-table penalty made the same failure cost less on the
    # receptor that failed more often — the less complete the evidence, the more
    # forgiving the aggregation. FR-006 requires every ligand to be aggregated
    # "on the same denominator", which per-table values never satisfied.
    #
    # With an inventory the penalty counts what the run intended to dock. Without
    # one it can only count what was observed, so the basis is recorded alongside
    # the value: a reader cannot otherwise tell which question was answered.
    if expected_ligands is not None:
        absence_penalty = len(expected_ligands) + 1
        absence_penalty_basis = "inventory_plus_one"
    else:
        absence_penalty = max(len(ranking) for ranking in rankings.values()) + 1
        absence_penalty_basis = "max_scored_plus_one"

    receptor_columns = []
    coverage = pd.Series(0, index=consensus.index, dtype=int)
    failed_by_receptor: dict[str, int] = {}
    for name, ranking in rankings.items():
        rank_column = f"Rank_{safe_name(name)}"
        receptor_columns.append(rank_column)
        rank_series = ranking.set_index("Ligand")["Rank"].reindex(consensus.index)
        coverage += rank_series.notna().astype(int)
        failed_by_receptor[name] = int(rank_series.isna().sum())
        # A docking failure is scored worse than every ligand successfully ranked
        # anywhere in the run, rather than being rewarded by omission.
        consensus[rank_column] = rank_series.fillna(absence_penalty)

    consensus["Best_Receptor_Rank"] = consensus[receptor_columns].min(axis=1)
    consensus["Median_Receptor_Rank"] = consensus[receptor_columns].median(axis=1)
    # How much of the run each ligand actually got through. A ligand at the
    # absence penalty in three of four receptors is not the same finding as one
    # that genuinely ranked mid-table everywhere.
    consensus["coverage_succeeded"] = coverage
    consensus["coverage_attempted"] = len(rankings)
    consensus["Consensus_RRF_Score"] = sum(
        1.0 / (args.rrf_k + consensus[column]) for column in receptor_columns
    )

    support_columns = []
    for name, ranking in rankings.items():
        rank_column = f"Rank_{safe_name(name)}"
        cutoff = max(1, int(len(ranking) * args.support_fraction))
        support_column = f"Top_{args.support_fraction:g}_Support_{safe_name(name)}"
        consensus[support_column] = consensus[rank_column] <= cutoff
        support_columns.append(support_column)
    consensus["Top_Support_Votes"] = consensus[support_columns].sum(axis=1).astype(int)
    consensus["Supporting_Receptors"] = consensus.apply(
        lambda row: ", ".join(
            name
            for name in rankings
            if row[f"Top_{args.support_fraction:g}_Support_{safe_name(name)}"]
        ),
        axis=1,
    )

    # Tail ties are routine — every ligand that failed in every receptor shares a
    # score — and pandas' default sort is not stable, so equal scores could come
    # back in a different order run to run. Ligand is the explicit tie-break.
    consensus = consensus.reset_index()
    consensus = consensus.sort_values(
        ["Consensus_RRF_Score", "Ligand"],
        ascending=[False, True],
        kind="stable",
    ).reset_index(drop=True)
    consensus["Consensus_Rank"] = consensus.index + 1

    args.out_dir.mkdir(parents=True, exist_ok=True)
    prepare_out_dir(args.out_dir, args.overwrite)
    all_path = args.out_dir / "consensus_all.parquet"
    csv_path = args.out_dir / "consensus_top_hits.csv"
    active_path = args.out_dir / "consensus_actives.parquet"
    decoy_path = args.out_dir / "consensus_decoys.parquet"
    manifest_path = args.out_dir / "consensus_manifest.csv"
    run_manifest_path = args.out_dir / "manifest.json"
    failures_path = args.out_dir / "failures.csv"

    shortlist = consensus.head(args.csv_n)
    consensus.to_parquet(all_path, index=False)
    shortlist.to_csv(csv_path, index=False)

    artifacts = [
        {
            "path": all_path.name,
            "ordered_by": "Consensus_RRF_Score",
            "direction": "higher_is_better",
            "row_count": len(consensus),
        },
        {
            "path": csv_path.name,
            "ordered_by": "Consensus_RRF_Score",
            "direction": "higher_is_better",
            "row_count": len(shortlist),
        },
    ]

    # Label-derived evaluation inputs exist to be fed back to enrichment, which
    # is a benchmark step. With no labels there is nothing to write, and writing
    # them empty would misreport an unlabeled run as one with zero actives.
    if labels_present:
        consensus[consensus["Class"] == "active"].to_parquet(active_path, index=False)
        consensus[consensus["Class"] == "decoy"].to_parquet(decoy_path, index=False)
        artifacts += [
            {
                "path": active_path.name,
                "ordered_by": "Consensus_RRF_Score",
                "direction": "higher_is_better",
                "row_count": int((consensus["Class"] == "active").sum()),
            },
            {
                "path": decoy_path.name,
                "ordered_by": "Consensus_RRF_Score",
                "direction": "higher_is_better",
                "row_count": int((consensus["Class"] == "decoy").sum()),
            },
        ]

    pd.DataFrame(
        {
            "Receptor": list(rankings),
            "Score_Column": [score_columns[name] for name in rankings],
            "Scored_Ligands": [len(rankings[name]) for name in rankings],
            "Failed_Ligands": [failed_by_receptor[name] for name in rankings],
        }
    ).to_csv(manifest_path, index=False)

    # A ligand that failed in every receptor is absent from every input table and
    # is therefore invisible here by construction. Only an inventory of what the
    # run intended to dock can reveal it, so without one the manifest says the
    # question could not be answered rather than implying nothing failed.
    if expected_ligands is None:
        total_failure_reporting = "unavailable_no_inventory"
        total_failures: list[str] = []
    else:
        total_failure_reporting = "complete"
        seen = set(consensus["Ligand"])
        total_failures = [
            identifier for identifier in expected_ligands if identifier not in seen
        ]
        pd.DataFrame(
            {
                "Ligand": total_failures,
                "Reason": ["absent from every receptor"] * len(total_failures),
            }
        ).to_csv(failures_path, index=False)

    run_manifest = {
        "engine": ENGINE_BY_SCORE_COLUMN.get(resolved_score_column),
        "method": "reciprocal_rank_fusion",
        "rrf_k": args.rrf_k,
        "score_column": resolved_score_column,
        "score_direction": "lower_is_better" if lower_is_better else "higher_is_better",
        "labels_present": labels_present,
        "tie_break": ["Ligand"],
        "receptor_order": list(rankings),
        # The penalty is a free parameter of the aggregation: a shortlist cannot
        # be audited without knowing what an absence cost, or which rule set it.
        "absence_penalty": absence_penalty,
        "absence_penalty_basis": absence_penalty_basis,
        "expected_ligand_inventory": (
            str(args.expected_ligands) if args.expected_ligands else None
        ),
        "total_failure_reporting": total_failure_reporting,
        "artifacts": artifacts,
        # Refinement is a later CTMD stage that writes its own files; it never
        # rewrites the consensus shortlist. Stated rather than left to be
        # inferred from an absent refined_top<M>.csv.
        "refinement_ran": False,
    }
    run_manifest_path.write_text(json.dumps(run_manifest, indent=2) + "\n")

    print(f"Receptors fused: {len(rankings)}")
    print(f"Receptor order: {', '.join(rankings)}")
    print(f"Unique ligands: {len(consensus)}")
    if labels_present:
        print(f"Actives: {(consensus['Class'] == 'active').sum()}")
        print(f"Decoys: {(consensus['Class'] == 'decoy').sum()}")
        print(f"Active enrichment input: {active_path}")
        print(f"Decoy enrichment input: {decoy_path}")
    else:
        print("Labels: none present; evaluation inputs skipped.")
    print(f"Consensus: {all_path}")
    print(f"Top hits CSV: {csv_path}")
    if len(shortlist) < args.csv_n:
        print(
            f"Note: {len(shortlist)} rows written, fewer than the {args.csv_n} requested "
            "by --csv_n. The consensus covers every ligand; there were no more to list."
        )
    print(f"Manifest: {manifest_path}")
    print(f"Run manifest: {run_manifest_path}")
    if expected_ligands is None:
        print(
            "Total-failure reporting: unavailable (no --expected_ligands inventory). "
            "Ligands that failed in every receptor cannot be seen from these inputs."
        )
    else:
        print(f"Total failures: {len(total_failures)} -> {failures_path}")


if __name__ == "__main__":
    main()
