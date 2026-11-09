#!/usr/bin/env python3
"""Regenerate the synthetic ranking fixtures.

Every table this writes is CONSTRUCTED. No docking run produced any of it. The
scores are arithmetic (``-12.0 + 0.1 * (rank - 1)``), chosen so the rank
structure each test depends on is stated here rather than discovered by reading
Parquet bytes.

Run from the repository root::

    python test_files/synthetic_ranking_fixtures/make_fixtures.py

The output is deterministic: no RNG, no clock, no filesystem ordering. Rerunning
reproduces the same bytes, which is what lets the T007 baseline digests mean
anything.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent

# Ligand identifiers are deliberately not ZINC/CHEMBL-shaped. Constructed data
# should never be mistakable for a real library.
FILLERS = [f"SYN{index:04d}" for index in range(1, 11)] + [
    f"SYN{index:04d}" for index in range(15, 21)
]
A, B, C, D = "SYN0011", "SYN0012", "SYN0013", "SYN0014"
LIGANDS = sorted(FILLERS + [A, B, C, D])

# Ligands in the expected-ligand inventory that appear in NO receptor table.
# Consensus cannot see these by construction (contract §4) — they exist to prove
# the inventory path reports them and the no-inventory path admits it cannot.
NEVER_DOCKED = ["SYN0021", "SYN0022"]

F = FILLERS

# Each receptor is an ordered list of tie groups, best first. A group with two
# members means both ligands carry the same score and therefore the same
# ``method="min"`` rank.
#
# The A/B and C/D structures are the point of this fixture:
#
#   A ranks (6, 13, 20, 9) and B ranks (20, 13, 6, 9) — the same multiset in a
#   different order. Their RRF sums are equal in exact arithmetic but are
#   accumulated in a different association depending on --input order, so they
#   can differ by an ulp and swap places. This is what makes sorting the
#   receptor inputs (T019) necessary rather than cosmetic, and it needs three or
#   more receptors: with two, a + b == b + a exactly.
#
#   C and D hold identical ranks in every receptor, so their RRF scores are
#   exactly equal. No accumulation order can separate them; only the explicit
#   Ligand tie-break (T020) makes their relative order stable.
RECEPTOR_ORDERS: dict[str, list[list[str]]] = {
    "rec_01": [
        [F[0]], [F[1]], [F[2]], [F[3]], [F[4]],
        [A],
        [F[5]], [F[6]], [F[7]],
        [C, D],
        [F[8]], [F[9]], [F[10]], [F[11]], [F[12]], [F[13]], [F[14]], [F[15]],
        [B],
    ],
    "rec_02": [
        [F[0]], [F[1]], [F[2]], [F[3]], [F[4]], [F[5]],
        [F[6]], [F[7]], [F[8]], [F[9]], [F[10]], [F[11]],
        [A, B],
        [C, D],
        [F[12]], [F[13]], [F[14]], [F[15]],
    ],
    "rec_03": [
        [F[0]], [F[1]], [F[2]], [F[3]], [F[4]],
        [B],
        [F[5]], [F[6]], [F[7]],
        [C, D],
        [F[8]], [F[9]], [F[10]], [F[11]], [F[12]], [F[13]], [F[14]], [F[15]],
        [A],
    ],
    "rec_04": [
        [F[0]], [F[1]], [F[2]], [F[3]], [F[4]], [F[5]], [F[6]], [F[7]],
        [A, B],
        [C, D],
        [F[8]], [F[9]], [F[10]], [F[11]], [F[12]], [F[13]], [F[14]], [F[15]],
    ],
}

# Actives for the labeled benchmark fixture. Both classes must be present or the
# pre-change code rejects the table outright.
ACTIVES = {F[0], F[1], F[2], F[3], F[4], A}


def scored_rows(groups: list[list[str]]) -> list[tuple[str, float]]:
    """Flatten tie groups into (ligand, score) pairs.

    Rank is the 1-based position of the group's first member, matching pandas'
    ``method="min"``. Score steps by 0.1 per rank so tied ligands share a score
    and distinct ranks never collide at one decimal place.
    """
    rows: list[tuple[str, float]] = []
    position = 0
    for group in groups:
        rank = position + 1
        score = round(-12.0 + 0.1 * (rank - 1), 2)
        for ligand in sorted(group):
            rows.append((ligand, score))
        position += len(group)
    return rows


def write_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False)
    print(f"wrote {path.relative_to(HERE.parent.parent)}  ({len(frame)} rows)")


def build_vina(groups: list[list[str]], receptor: str) -> pd.DataFrame:
    """The shape extract_vina_results.py writes: Ligand, Vina_Score, File_Path."""
    rows = scored_rows(groups)
    frame = pd.DataFrame(
        {
            "Ligand": [ligand for ligand, _ in rows],
            "Vina_Score": [score for _, score in rows],
            "File_Path": [
                f"synthetic/{receptor}/{ligand}_out.pdbqt" for ligand, _ in rows
            ],
        }
    )
    return frame.sort_values("Vina_Score", ascending=True).reset_index(drop=True)


def build_enrichment_ranked(groups: list[list[str]], receptor: str) -> pd.DataFrame:
    """The shape calculate_enrichment.py writes: adds Class and Is_Active."""
    frame = build_vina(groups, receptor)
    frame["Raw_Ligand"] = frame["Ligand"]
    frame["Class"] = [
        "active" if ligand in ACTIVES else "decoy" for ligand in frame["Ligand"]
    ]
    frame["Is_Active"] = frame["Class"] == "active"
    frame["Rank"] = frame.index + 1
    frame["Cumulative_Actives"] = frame["Is_Active"].cumsum()
    frame["Cumulative_Decoys"] = frame["Rank"] - frame["Cumulative_Actives"]
    return frame


def build_dock6(groups: list[list[str]], receptor: str) -> pd.DataFrame:
    """One row per POSE, as extract_dock6_results.py writes.

    Three poses per ligand. Pose 1 carries the ligand's real score; poses 2 and 3
    are progressively worse, so collapsing on the ligand column alone must keep
    pose 1 and must not be reachable by keeping the first row of the file.
    """
    rows = scored_rows(groups)
    records = []
    for ligand, score in rows:
        for pose_index in (1, 2, 3):
            records.append(
                {
                    "Ligand": ligand,
                    "Grid_Score": round(score + 0.7 * (pose_index - 1), 2),
                    "Pose_Index": pose_index,
                    "Source_File": f"synthetic/{receptor}/{ligand}_scored.mol2",
                    "Mol2_Start_Line": 1 + (pose_index - 1) * 120,
                    "Mol2_End_Line": pose_index * 120,
                }
            )
    frame = pd.DataFrame.from_records(records)
    # Worst pose first, so a de-duplication that ignores the sort would keep the
    # wrong row and the test would catch it.
    return frame.sort_values(
        ["Ligand", "Grid_Score"], ascending=[True, False]
    ).reset_index(drop=True)


# --- uneven completion (T033, FR-006, SC-013) --------------------------------
#
# Three receptors that scored DIFFERENT numbers of ligands. That difference is
# the whole point: under the superseded per-receptor rule the penalty was the
# length of the table a ligand was absent from, so the same failure cost a
# different amount in each receptor. Here that would be 10 in rec_b and 8 in
# rec_c — and the shared rule gives 13 in both, so the two rules are
# numerically distinguishable rather than coincidentally equal.
#
# SYN0007 is absent from two of the three receptors, which is what makes
# "the same penalty in every column" an assertion with something to compare.
UNEVEN_LIGANDS = [f"SYN{index:04d}" for index in range(1, 13)]
UNEVEN_MISSING = {
    "unv_a": [],                                    # scores all 12
    "unv_b": ["SYN0007", "SYN0011", "SYN0012"],     # scores 9
    "unv_c": ["SYN0007", "SYN0008", "SYN0009", "SYN0010", "SYN0011"],  # scores 7
}
# A larger inventory than any receptor saw, so len(inventory) + 1 and
# max(scored) + 1 cannot collide: 15 + 1 = 16 against 12 + 1 = 13.
UNEVEN_INVENTORY = UNEVEN_LIGANDS + ["SYN0013", "SYN0014", "SYN0015"]


def build_uneven(receptor: str) -> pd.DataFrame:
    kept = [lig for lig in UNEVEN_LIGANDS if lig not in UNEVEN_MISSING[receptor]]
    groups = [[ligand] for ligand in kept]
    return build_vina(groups, receptor)


def main() -> None:
    # --- unlabeled_vina: the P1 case. Four receptors; see RECEPTOR_ORDERS. ---
    for receptor, groups in RECEPTOR_ORDERS.items():
        write_parquet(
            build_vina(groups, receptor),
            HERE / "unlabeled_vina" / f"{receptor}_docking_summary.parquet",
        )

    # --- labeled_benchmark: same ligands, plus labels. Directory form, because
    # resolve_parquet() accepts a directory containing enrichment_ranked.parquet
    # and the benchmark path is the one that exercises it. ---
    for receptor in ("rec_01", "rec_02"):
        write_parquet(
            build_enrichment_ranked(RECEPTOR_ORDERS[receptor], receptor),
            HERE / "labeled_benchmark" / receptor / "enrichment_ranked.parquet",
        )

    # --- unlabeled_dock6: multi-pose collapse. ---
    for receptor in ("rec_01", "rec_02"):
        write_parquet(
            build_dock6(RECEPTOR_ORDERS[receptor], receptor),
            HERE / "unlabeled_dock6" / f"{receptor}_dock6_summary.parquet",
        )

    # --- partial_failure: rec_02 docked only half the library. ---
    write_parquet(
        build_vina(RECEPTOR_ORDERS["rec_01"], "rec_01"),
        HERE / "partial_failure" / "rec_01_docking_summary.parquet",
    )
    half = set(LIGANDS[:10])
    partial_groups = [
        [ligand for ligand in group if ligand in half]
        for group in RECEPTOR_ORDERS["rec_02"]
    ]
    partial_groups = [group for group in partial_groups if group]
    write_parquet(
        build_vina(partial_groups, "rec_02"),
        HERE / "partial_failure" / "rec_02_docking_summary.parquet",
    )

    # The inventory is a superset of what any receptor saw: it carries the two
    # ligands that appear nowhere, which is the only way consensus can learn
    # they exist at all.
    inventory = HERE / "partial_failure" / "expected_ligands.txt"
    inventory.write_text("\n".join(LIGANDS + NEVER_DOCKED) + "\n")
    print(f"wrote {inventory.relative_to(HERE.parent.parent)}  ({len(LIGANDS) + 2} ids)")

    # --- uneven_completion: receptors of differing lengths. ---
    for receptor in UNEVEN_MISSING:
        write_parquet(
            build_uneven(receptor),
            HERE / "uneven_completion" / f"{receptor}_docking_summary.parquet",
        )
    uneven_inventory = HERE / "uneven_completion" / "expected_ligands.txt"
    uneven_inventory.write_text("\n".join(UNEVEN_INVENTORY) + "\n")
    print(
        f"wrote {uneven_inventory.relative_to(HERE.parent.parent)}  "
        f"({len(UNEVEN_INVENTORY)} ids)"
    )


if __name__ == "__main__":
    main()
