#!/usr/bin/env python3
"""Write rank-1 pocket definitions for docking from a ranked-pocket CSV."""

from __future__ import annotations

import csv
import re
from pathlib import Path


def numeric_rank(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Pocket rank must be an integer, found {value!r}.") from exc


def safe_filename(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    return cleaned or "receptor"


def float_value(row: dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Pocket row is missing a numeric {key!r} value: {row}") from exc


def score_from_row(row: dict[str, str], method: str) -> tuple[str, float]:
    score_key = "p2rank_score" if method == "p2rank" else "score_sum"
    return score_key, float_value(row, score_key)


def dump_yaml(path: Path, content: dict) -> None:
    try:
        import yaml
    except ImportError as exc:
        raise ImportError("Writing docking pocket definitions requires pyyaml.") from exc
    with path.open("w") as handle:
        yaml.safe_dump(content, handle, sort_keys=False)


def export_best_pockets(ranked_pockets_csv: Path, method: str, output_dir: Path) -> tuple[Path, int]:
    """Export one interoperable docking-box YAML for each receptor's rank-1 pocket."""
    if method not in {"grasp", "p2rank"}:
        raise ValueError("method must be either 'grasp' or 'p2rank'.")

    with ranked_pockets_csv.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    rank_key = "pocket_rank" if method == "p2rank" else "rank"
    best_by_receptor = {}
    for row in rows:
        receptor = row.get("receptor", "")
        if not receptor:
            raise ValueError(f"Pocket row has no receptor value: {row}")
        if receptor not in best_by_receptor or numeric_rank(row[rank_key]) < numeric_rank(best_by_receptor[receptor][rank_key]):
            best_by_receptor[receptor] = row

    definitions_dir = output_dir / "docking_boxes"
    definitions_dir.mkdir(parents=True, exist_ok=True)
    index_rows = []
    for receptor in sorted(best_by_receptor):
        row = best_by_receptor[receptor]
        score_name, score = score_from_row(row, method)
        score = round(score, 2)
        center = [round(float_value(row, key), 2) for key in ("center_x", "center_y", "center_z")]
        size = [round(float_value(row, key), 2) for key in ("size_x", "size_y", "size_z")]
        definition_path = definitions_dir / f"{safe_filename(receptor)}_best_pocket.yaml"
        dump_yaml(
            definition_path,
            {
                "schema_version": 1,
                "receptor": receptor,
                "pocket": {"method": method, "rank": numeric_rank(row[rank_key]), score_name: score},
                "box": {
                    "center": list(center),
                    "size": list(size),
                    "center_string": " ".join(f"{value:.2f}" for value in center),
                    "size_string": " ".join(f"{value:.2f}" for value in size),
                },
                "vina": {"box": {"center": " ".join(f"{value:.2f}" for value in center), "size": " ".join(f"{value:.2f}" for value in size)}},
                "dock6": {"box": {"center": list(center), "size": list(size)}},
            },
        )
        index_rows.append(
            {
                "receptor": receptor,
                "method": method,
                "pocket_rank": numeric_rank(row[rank_key]),
                "score_name": score_name,
                "score": f"{score:.2f}",
                "center_x": f"{center[0]:.2f}",
                "center_y": f"{center[1]:.2f}",
                "center_z": f"{center[2]:.2f}",
                "size_x": f"{size[0]:.2f}",
                "size_y": f"{size[1]:.2f}",
                "size_z": f"{size[2]:.2f}",
                "definition_file": str(definition_path),
            }
        )

    index_path = definitions_dir / "best_pockets.csv"
    with index_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "receptor", "method", "pocket_rank", "score_name", "score",
                "center_x", "center_y", "center_z", "size_x", "size_y", "size_z", "definition_file",
            ],
        )
        writer.writeheader()
        writer.writerows(index_rows)
    return index_path, len(index_rows)
