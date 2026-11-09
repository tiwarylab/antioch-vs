#!/usr/bin/env python3
"""Measure kinase A-loop foldedness for AF2RAVE-selected structures.

This standalone post-analysis helper reads AF2RAVE/ANTIOCH-VS output files such
as filtered_structures.txt, cluster_centers.txt, and cluster_assignments.csv,
then measures an A-loop foldedness distance. For the JAK2 model mapping used in
the ANTIOCH-VS examples, the default metric is a D102-K162 salt-bridge style
minimum heavy-atom distance:

  min(dist(D102 OD1/OD2, K162 NZ))

Structures with distance <= cutoff are labeled "folded"; the rest are labeled
"extended". Cluster centers are overlaid on the output plots.
"""

from __future__ import annotations

import argparse
import shlex
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Atom:
    name: str
    resname: str
    chain: str
    resseq: int
    xyz: np.ndarray


def parse_pdb_atoms(path: Path) -> list[Atom]:
    atoms = []
    with path.open(errors="replace") as handle:
        for line in handle:
            if line[:6].strip() not in {"ATOM", "HETATM"}:
                continue
            try:
                atoms.append(
                    Atom(
                        name=line[12:16].strip(),
                        resname=line[17:20].strip(),
                        chain=line[21].strip(),
                        resseq=int(line[22:26]),
                        xyz=np.array(
                            [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                            dtype=float,
                        ),
                    )
                )
            except ValueError:
                continue
    return atoms


def parse_cif_atoms(path: Path) -> list[Atom]:
    atoms = []
    headers = []
    in_atom_site = False
    with path.open(errors="replace") as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith("_atom_site."):
                headers.append(line.split()[0])
                in_atom_site = True
                continue
            if not in_atom_site:
                continue
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                if line.startswith("#") and atoms:
                    break
                continue
            try:
                values = shlex.split(line)
                row = dict(zip(headers, values))
                atom_name = row.get("_atom_site.auth_atom_id") or row.get("_atom_site.label_atom_id")
                resname = row.get("_atom_site.auth_comp_id") or row.get("_atom_site.label_comp_id")
                chain = row.get("_atom_site.auth_asym_id") or row.get("_atom_site.label_asym_id") or ""
                resseq = row.get("_atom_site.auth_seq_id") or row.get("_atom_site.label_seq_id")
                atoms.append(
                    Atom(
                        name=str(atom_name).strip(),
                        resname=str(resname).strip(),
                        chain=str(chain).strip(),
                        resseq=int(resseq),
                        xyz=np.array(
                            [
                                float(row["_atom_site.Cartn_x"]),
                                float(row["_atom_site.Cartn_y"]),
                                float(row["_atom_site.Cartn_z"]),
                            ],
                            dtype=float,
                        ),
                    )
                )
            except (KeyError, TypeError, ValueError):
                continue
    return atoms


def parse_structure_atoms(path: Path) -> list[Atom]:
    if path.suffix.lower() in {".cif", ".mmcif"}:
        return parse_cif_atoms(path)
    return parse_pdb_atoms(path)


def read_path_list(path: Path) -> list[Path]:
    if not path.exists():
        return []
    return [Path(line.strip()) for line in path.read_text().splitlines() if line.strip()]


def basename_key(path: str | Path) -> str:
    return Path(str(path)).name


def load_cluster_assignments(cluster_assignments_csv: Path | None) -> dict[str, dict]:
    if not cluster_assignments_csv or not cluster_assignments_csv.exists():
        return {}
    df = pd.read_csv(cluster_assignments_csv)
    assignments = {}
    for _, row in df.iterrows():
        key = str(row.get("pdb", ""))
        if not key:
            continue
        value = {
            "cluster": row.get("cluster", ""),
            "distance_to_centroid": row.get("distance_to_centroid", ""),
            "is_center": bool(row.get("is_center", False)),
        }
        assignments[key] = value
        assignments[basename_key(key)] = value
    return assignments


def atoms_for_residue(atoms: list[Atom], resseq: int, atom_names: list[str], chain: str | None) -> list[Atom]:
    wanted = set(atom_names)
    return [
        atom
        for atom in atoms
        if atom.resseq == resseq
        and atom.name in wanted
        and (chain is None or atom.chain == chain)
    ]


def parse_atom_names(value: str) -> list[str]:
    atoms = [item.strip() for item in value.split(",") if item.strip()]
    if not atoms:
        raise argparse.ArgumentTypeError("Atom list cannot be empty.")
    return atoms


def measure_min_distance(
    atoms: list[Atom],
    res1: int,
    atoms1: list[str],
    res2: int,
    atoms2: list[str],
    chain: str | None,
) -> tuple[float, str, str]:
    group1 = atoms_for_residue(atoms, res1, atoms1, chain)
    group2 = atoms_for_residue(atoms, res2, atoms2, chain)
    if not group1:
        raise ValueError(f"Missing atoms for residue {res1}: {','.join(atoms1)}")
    if not group2:
        raise ValueError(f"Missing atoms for residue {res2}: {','.join(atoms2)}")

    best = None
    for atom1 in group1:
        for atom2 in group2:
            dist = float(np.linalg.norm(atom1.xyz - atom2.xyz))
            if best is None or dist < best[0]:
                best = (dist, atom1.name, atom2.name)
    if best is None:
        raise ValueError("No atom pairs found for A-loop distance.")
    return best


def state_from_distance(distance: float, cutoff: float) -> str:
    return "folded" if distance <= cutoff else "extended"


def compute_metrics(args) -> pd.DataFrame:
    structures = read_path_list(args.structures_list)
    if not structures:
        raise FileNotFoundError(f"No structures found in list: {args.structures_list}")

    centers = set(map(str, read_path_list(args.cluster_centers)))
    centers_by_name = {basename_key(path) for path in centers}
    assignments = load_cluster_assignments(args.cluster_assignments)

    rows = []
    for pdb in structures:
        try:
            atoms = parse_structure_atoms(pdb)
            distance, atom1, atom2 = measure_min_distance(
                atoms,
                args.res1 + args.residue_offset,
                args.atoms1,
                args.res2 + args.residue_offset,
                args.atoms2,
                args.chain,
            )
            assignment = assignments.get(str(pdb), assignments.get(basename_key(pdb), {}))
            rows.append(
                {
                    "pdb": str(pdb),
                    "pdb_name": pdb.name,
                    "aloop_distance_A": distance,
                    "aloop_state": state_from_distance(distance, args.cutoff),
                    "residue1": args.res1 + args.residue_offset,
                    "residue2": args.res2 + args.residue_offset,
                    "atom1": atom1,
                    "atom2": atom2,
                    "is_center": str(pdb) in centers or pdb.name in centers_by_name or bool(assignment.get("is_center", False)),
                    "cluster": assignment.get("cluster", ""),
                    "distance_to_centroid": assignment.get("distance_to_centroid", ""),
                    "error": "",
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "pdb": str(pdb),
                    "pdb_name": pdb.name,
                    "aloop_distance_A": np.nan,
                    "aloop_state": "error",
                    "residue1": args.res1 + args.residue_offset,
                    "residue2": args.res2 + args.residue_offset,
                    "atom1": "",
                    "atom2": "",
                    "is_center": str(pdb) in centers or pdb.name in centers_by_name,
                    "cluster": "",
                    "distance_to_centroid": "",
                    "error": str(exc),
                }
            )
    return pd.DataFrame(rows)


def write_center_summary(df: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    centers = df[(df["is_center"]) & (df["error"] == "")].copy()
    centers = centers.sort_values(["aloop_state", "aloop_distance_A"], na_position="last")
    centers.to_csv(output_dir / "aloop_foldedness_cluster_centers.csv", index=False)
    return centers


def plot_foldedness(df: pd.DataFrame, output_dir: Path, args) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    valid = df[df["error"] == ""].copy()
    centers = valid[valid["is_center"]].copy()

    fig, ax = plt.subplots(figsize=(7.0, 4.2), layout="compressed")
    ax.hist(valid["aloop_distance_A"], bins=args.bins, color="#8fbcd4", edgecolor="white", alpha=0.9)
    ax.axvline(args.cutoff, color="black", linestyle="--", linewidth=1.2, label=f"cutoff = {args.cutoff:g} A")
    if not centers.empty:
        ymax = ax.get_ylim()[1]
        ax.scatter(
            centers["aloop_distance_A"],
            np.full(len(centers), ymax * 0.06),
            marker="H",
            facecolors="none",
            edgecolors="#c23b22",
            s=85,
            linewidths=1.4,
            label="Cluster centers",
            zorder=3,
        )
    ax.set_xlabel(args.xlabel)
    ax.set_ylabel("Structure count")
    ax.set_title(args.title)
    ax.legend()
    fig.savefig(output_dir / "aloop_foldedness_histogram.png", dpi=args.dpi)
    plt.close(fig)

    ranked = valid.sort_values("aloop_distance_A").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.2), layout="compressed")
    colors = ranked["aloop_state"].map({"folded": "#2b8cbe", "extended": "#f03b20"}).fillna("#999999")
    ax.scatter(np.arange(len(ranked)), ranked["aloop_distance_A"], c=colors, marker="x", alpha=0.7, s=20)
    center_ranked = ranked[ranked["is_center"]]
    if not center_ranked.empty:
        ax.scatter(
            center_ranked.index,
            center_ranked["aloop_distance_A"],
            marker="H",
            facecolors="none",
            edgecolors="#c23b22",
            s=95,
            linewidths=1.4,
            label="Cluster centers",
            zorder=3,
        )
    ax.axhline(args.cutoff, color="black", linestyle="--", linewidth=1.2, label=f"cutoff = {args.cutoff:g} A")
    ax.set_xlabel("Structures sorted by A-loop distance")
    ax.set_ylabel(args.xlabel)
    ax.set_title(f"{args.title}: ranked distances")
    ax.legend()
    fig.savefig(output_dir / "aloop_foldedness_ranked.png", dpi=args.dpi)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--af2rave_dir", required=True, type=Path, help="AF2RAVE/ANTIOCH-VS output directory.")
    parser.add_argument("--output_dir", type=Path, help="Output directory. Defaults to <af2rave_dir>/kinase_aloop_foldedness.")
    parser.add_argument("--structures_list", type=Path, help="Defaults to <af2rave_dir>/filtered_structures.txt.")
    parser.add_argument("--cluster_centers", type=Path, help="Defaults to <af2rave_dir>/cluster_centers.txt.")
    parser.add_argument("--cluster_assignments", type=Path, help="Defaults to <af2rave_dir>/cluster_assignments.csv.")
    parser.add_argument("--chain", help="Optional chain ID for atom selection.")
    parser.add_argument("--res1", type=int, default=102, help="First A-loop foldedness residue. Default: D102.")
    parser.add_argument("--res2", type=int, default=162, help="Second A-loop foldedness residue. Default: K162.")
    parser.add_argument("--atoms1", type=parse_atom_names, default=["OD1", "OD2"], help="Comma-separated atoms on residue 1. Default: OD1,OD2.")
    parser.add_argument("--atoms2", type=parse_atom_names, default=["NZ"], help="Comma-separated atoms on residue 2. Default: NZ.")
    parser.add_argument("--residue_offset", type=int, default=0, help="Added to both residue numbers before PDB/mmCIF lookup.")
    parser.add_argument("--cutoff", type=float, default=10.0, help="Folded cutoff in Angstrom. Distance <= cutoff is folded.")
    parser.add_argument("--bins", type=int, default=50)
    parser.add_argument("--xlabel", default="A-loop distance: min(D102 OD1/OD2, K162 NZ) (A)")
    parser.add_argument("--title", default="Kinase A-loop foldedness")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    args.af2rave_dir = args.af2rave_dir.expanduser().resolve()
    args.output_dir = (args.output_dir or args.af2rave_dir / "kinase_aloop_foldedness").expanduser().resolve()
    args.structures_list = args.structures_list or args.af2rave_dir / "filtered_structures.txt"
    args.cluster_centers = args.cluster_centers or args.af2rave_dir / "cluster_centers.txt"
    args.cluster_assignments = args.cluster_assignments or args.af2rave_dir / "cluster_assignments.csv"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = compute_metrics(args)
    df.to_csv(args.output_dir / "aloop_foldedness_metrics.csv", index=False)
    centers = write_center_summary(df, args.output_dir)
    plot_foldedness(df, args.output_dir, args)

    print(f"Structures analyzed: {len(df)}")
    print(f"Structures with errors: {int((df['error'] != '').sum())}")
    print("A-loop states:")
    print(df[df["error"] == ""]["aloop_state"].value_counts().to_string())
    print(f"Cluster centers: {len(centers)}")
    if not centers.empty:
        print("Cluster-center A-loop states:")
        print(centers["aloop_state"].value_counts().to_string())
    print(f"Metrics: {args.output_dir / 'aloop_foldedness_metrics.csv'}")
    print(f"Center summary: {args.output_dir / 'aloop_foldedness_cluster_centers.csv'}")
    print(f"Histogram: {args.output_dir / 'aloop_foldedness_histogram.png'}")
    print(f"Ranked plot: {args.output_dir / 'aloop_foldedness_ranked.png'}")


if __name__ == "__main__":
    main()
