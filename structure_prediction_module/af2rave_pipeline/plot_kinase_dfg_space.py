#!/usr/bin/env python3
"""Plot AF2RAVE-selected kinase structures in Dunbrack D1/D2 space.

This is a standalone post-analysis helper. It reads AF2RAVE/ANTIOCH-VS output
files such as filtered_structures.txt, cluster_centers.txt, feature_stats.csv,
and cluster_assignments.csv, then computes:

  D1 = distance(ChelX CA, DFG Phe CZ)
  D2 = distance(beta3 Lys CA, DFG Phe CZ)

The definitions match the D1/D2 plotting used in the AF2RAVE-Glide kinase
notebook, but this script operates directly on PDB files and does not require
the notebook-specific kinaseCVs.py helper.
"""

from __future__ import annotations

import argparse
import shlex
import shutil
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
    """Parse atom records from an mmCIF file using author chain/residue IDs."""
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
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        return parse_cif_atoms(path)
    return parse_pdb_atoms(path)


def atom_xyz(atoms: list[Atom], resseq: int, atom_name: str, chain: str | None = None) -> np.ndarray:
    matches = [
        atom.xyz
        for atom in atoms
        if atom.resseq == resseq
        and atom.name == atom_name
        and (chain is None or atom.chain == chain)
    ]
    if not matches:
        chain_text = f" chain {chain}" if chain else ""
        raise ValueError(f"Missing atom: residue {resseq}{chain_text} atom {atom_name}")
    return matches[0]


def distance_angstrom(atoms: list[Atom], resseq_a: int, atom_a: str, resseq_b: int, atom_b: str, chain: str | None) -> float:
    return float(np.linalg.norm(atom_xyz(atoms, resseq_a, atom_a, chain) - atom_xyz(atoms, resseq_b, atom_b, chain)))


def kabsch_rmsd(moving: np.ndarray, reference: np.ndarray) -> float:
    moving_centroid = moving.mean(axis=0)
    reference_centroid = reference.mean(axis=0)
    x = moving - moving_centroid
    y = reference - reference_centroid
    cov = x.T @ y
    v, _, wt = np.linalg.svd(cov)
    fix = np.diag([1.0, 1.0, np.sign(np.linalg.det(v @ wt))])
    rotation = v @ fix @ wt
    aligned = moving @ rotation + (reference_centroid - moving_centroid @ rotation)
    return float(np.sqrt(np.mean(np.sum((aligned - reference) ** 2, axis=1))))


def ca_points_by_resseq(path: Path, chain: str | None = None) -> dict[int, np.ndarray]:
    points = {}
    for atom in parse_structure_atoms(path):
        if atom.name == "CA" and (chain is None or atom.chain == chain):
            points[atom.resseq] = atom.xyz
    return points


def calculate_ca_rmsd(path: Path, reference: Path, chain: str | None = None) -> float:
    moving = ca_points_by_resseq(path, chain)
    ref = ca_points_by_resseq(reference, chain)
    common = sorted(set(moving) & set(ref))
    if len(common) < 3:
        raise ValueError(f"Too few common CA residues for RMSD: {path}")
    return kabsch_rmsd(
        np.array([moving[resid] for resid in common]),
        np.array([ref[resid] for resid in common]),
    )


def read_path_list(path: Path) -> list[Path]:
    if not path.exists():
        return []
    return [Path(line.strip()) for line in path.read_text().splitlines() if line.strip()]


def basename_key(path: str | Path) -> str:
    return Path(str(path)).name


def load_rmsd_map(feature_stats_csv: Path | None) -> dict[str, float]:
    if not feature_stats_csv or not feature_stats_csv.exists():
        return {}
    df = pd.read_csv(feature_stats_csv)
    if "pdb" not in df.columns or "rmsd" not in df.columns:
        return {}
    rmsd_map = {}
    for _, row in df.iterrows():
        try:
            rmsd = float(row["rmsd"])
        except (TypeError, ValueError):
            continue
        rmsd_map[str(row["pdb"])] = rmsd
        rmsd_map[basename_key(row["pdb"])] = rmsd
    return rmsd_map


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


def classify_dfg_state(d1: float, d2: float, d1_cutoff: float, d2_in_cutoff: float, d2_out_cutoff: float) -> str:
    if d1 > d1_cutoff and d2 < d2_out_cutoff:
        return "DFG-out-like"
    if d1 < d1_cutoff and d2 < d2_in_cutoff:
        return "DFG-intermediate-like"
    if d1 < d1_cutoff and d2 > d2_in_cutoff:
        return "DFG-in-like"
    return "unassigned"


def compute_metrics(args) -> pd.DataFrame:
    structures = read_path_list(args.structures_list)
    if not structures:
        raise FileNotFoundError(f"No structures found in list: {args.structures_list}")

    centers = set(map(str, read_path_list(args.cluster_centers)))
    centers_by_name = {basename_key(path) for path in centers}
    rmsd_map = load_rmsd_map(args.feature_stats)
    assignments = load_cluster_assignments(args.cluster_assignments)

    chelx = args.chelx if args.chelx is not None else args.chele + 4
    dfgphe = args.dfgphe if args.dfgphe is not None else args.dfgasp + 1

    rows = []
    for pdb in structures:
        try:
            atoms = parse_structure_atoms(pdb)
            d1 = distance_angstrom(atoms, chelx + args.residue_offset, "CA", dfgphe + args.residue_offset, "CZ", args.chain)
            d2 = distance_angstrom(atoms, args.sbridgek + args.residue_offset, "CA", dfgphe + args.residue_offset, "CZ", args.chain)
            rmsd = rmsd_map.get(str(pdb), rmsd_map.get(basename_key(pdb), np.nan))
            if (not np.isfinite(rmsd)) and args.reference_pdb:
                rmsd = calculate_ca_rmsd(pdb, args.reference_pdb, args.chain)
            assignment = assignments.get(str(pdb), assignments.get(basename_key(pdb), {}))
            rows.append(
                {
                    "pdb": str(pdb),
                    "pdb_name": pdb.name,
                    "D1_A": d1,
                    "D2_A": d2,
                    "RMSD_A": rmsd,
                    "is_center": str(pdb) in centers or pdb.name in centers_by_name or bool(assignment.get("is_center", False)),
                    "cluster": assignment.get("cluster", ""),
                    "distance_to_centroid": assignment.get("distance_to_centroid", ""),
                    "dfg_state": classify_dfg_state(d1, d2, args.d1_cutoff, args.d2_in_cutoff, args.d2_out_cutoff),
                    "is_xtal_reference": False,
                    "error": "",
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "pdb": str(pdb),
                    "pdb_name": pdb.name,
                    "D1_A": np.nan,
                    "D2_A": np.nan,
                    "RMSD_A": np.nan,
                    "is_center": str(pdb) in centers or pdb.name in centers_by_name,
                    "cluster": "",
                    "distance_to_centroid": "",
                    "dfg_state": "error",
                    "is_xtal_reference": False,
                    "error": str(exc),
                }
            )
    return pd.DataFrame(rows)


def compute_xtal_metric(args) -> dict | None:
    if args.xtal_structure is None:
        return None

    model_chelx = args.chelx if args.chelx is not None else args.chele + 4
    model_dfgphe = args.dfgphe if args.dfgphe is not None else args.dfgasp + 1

    xtal_chain = args.xtal_chain if args.xtal_chain is not None else args.chain
    xtal_chelx = args.xtal_chelx
    if xtal_chelx is None and args.xtal_chele is not None:
        xtal_chelx = args.xtal_chele + 4
    if xtal_chelx is None:
        xtal_chelx = model_chelx

    xtal_dfgphe = args.xtal_dfgphe
    if xtal_dfgphe is None and args.xtal_dfgasp is not None:
        xtal_dfgphe = args.xtal_dfgasp + 1
    if xtal_dfgphe is None:
        xtal_dfgphe = model_dfgphe

    xtal_sbridgek = args.xtal_sbridgek if args.xtal_sbridgek is not None else args.sbridgek

    try:
        atoms = parse_structure_atoms(args.xtal_structure)
        d1 = distance_angstrom(
            atoms,
            xtal_chelx + args.xtal_residue_offset,
            "CA",
            xtal_dfgphe + args.xtal_residue_offset,
            "CZ",
            xtal_chain,
        )
        d2 = distance_angstrom(
            atoms,
            xtal_sbridgek + args.xtal_residue_offset,
            "CA",
            xtal_dfgphe + args.xtal_residue_offset,
            "CZ",
            xtal_chain,
        )
        state = classify_dfg_state(d1, d2, args.d1_cutoff, args.d2_in_cutoff, args.d2_out_cutoff)
        error = ""
    except Exception as exc:
        d1 = np.nan
        d2 = np.nan
        state = "error"
        error = str(exc)

    return {
        "pdb": str(args.xtal_structure),
        "pdb_name": args.xtal_label,
        "D1_A": d1,
        "D2_A": d2,
        "RMSD_A": np.nan,
        "is_center": False,
        "cluster": "xtal",
        "distance_to_centroid": "",
        "dfg_state": state,
        "is_xtal_reference": True,
        "error": error,
    }


def write_state_representatives(df: pd.DataFrame, output_dir: Path, top_n: int) -> pd.DataFrame:
    rows = []
    center_df = df[(df["is_center"]) & (~df.get("is_xtal_reference", False)) & (df["error"] == "")].copy()
    for state, group in center_df.groupby("dfg_state", sort=True):
        if state in {"error", "unassigned"}:
            continue
        sort_cols = ["RMSD_A"] if group["RMSD_A"].notna().any() else ["distance_to_centroid"]
        ranked = group.sort_values(sort_cols, na_position="last").head(top_n).copy()
        ranked.insert(0, "state_rank", range(1, len(ranked) + 1))
        rows.append(ranked)

    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    out.to_csv(output_dir / "kinase_d1_d2_state_representatives.csv", index=False)
    return out


def copy_representatives(representatives: pd.DataFrame, output_dir: Path) -> None:
    if representatives.empty:
        return
    rep_dir = output_dir / "state_representative_pdbs"
    rep_dir.mkdir(parents=True, exist_ok=True)
    for _, row in representatives.iterrows():
        state = str(row["dfg_state"]).replace("-", "_")
        dst = rep_dir / f"{state}_{int(row['state_rank']):02d}_{Path(row['pdb']).name}"
        shutil.copy2(row["pdb"], dst)


def plot_d1_d2(df: pd.DataFrame, output_dir: Path, args) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    valid = df[df["error"] == ""].copy()
    xtal = valid[valid.get("is_xtal_reference", False)].copy()
    centers = valid[(valid["is_center"]) & (~valid.get("is_xtal_reference", False))].copy()
    models = valid[~valid.get("is_xtal_reference", False)].copy()

    fig, ax = plt.subplots(figsize=(6, 5.2), layout="compressed")
    color = models["RMSD_A"].to_numpy(dtype=float)
    if np.isfinite(color).any():
        scatter = ax.scatter(
            models["D1_A"],
            models["D2_A"],
            marker="x",
            c=color,
            s=16,
            alpha=0.65,
            cmap=args.cmap,
            vmin=args.rmsd_vmin,
            vmax=args.rmsd_vmax if args.rmsd_vmax is not None else np.nanmax(color),
        )
        cbar = fig.colorbar(scatter, ax=ax)
        cbar.set_label("RMSD (A)")
    else:
        ax.scatter(models["D1_A"], models["D2_A"], marker="x", s=16, alpha=0.65)

    ax.plot(
        centers["D1_A"],
        centers["D2_A"],
        "H",
        markerfacecolor="none",
        markeredgecolor="blue",
        markersize=8,
        linestyle="None",
        label="Cluster centers",
    )
    if not xtal.empty:
        ax.plot(
            xtal["D1_A"],
            xtal["D2_A"],
            marker="*",
            markerfacecolor=args.xtal_color,
            markeredgecolor="black",
            markersize=14,
            linestyle="None",
            label=args.xtal_label,
        )

    ax.axvline(args.d1_cutoff, linestyle="--", color="k", linewidth=1)
    ax.hlines(args.d2_in_cutoff, args.xlim[0], args.d1_cutoff, linestyles="--", color="k", linewidth=1)
    ax.hlines(args.d2_out_cutoff, args.d1_cutoff, args.xlim[1], linestyles="--", color="k", linewidth=1)
    ax.set_xlim(args.xlim)
    ax.set_ylim(args.ylim)
    ax.set_xlabel("D1: ChelX CA - DFG Phe CZ (A)")
    ax.set_ylabel("D2: beta3 Lys CA - DFG Phe CZ (A)")
    ax.set_title(args.title)
    ax.legend()
    fig.savefig(output_dir / "kinase_d1_d2_clusters.png", dpi=args.dpi)
    plt.close(fig)


def parse_range(value: str) -> tuple[float, float]:
    parts = [float(x.strip()) for x in value.split(",")]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("Range must be two comma-separated numbers, e.g. 2,22")
    return parts[0], parts[1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--af2rave_dir", required=True, type=Path, help="AF2RAVE/ANTIOCH-VS output directory.")
    parser.add_argument("--output_dir", type=Path, help="Output directory. Defaults to <af2rave_dir>/kinase_d1_d2_analysis.")
    parser.add_argument("--structures_list", type=Path, help="Defaults to <af2rave_dir>/filtered_structures.txt.")
    parser.add_argument("--cluster_centers", type=Path, help="Defaults to <af2rave_dir>/cluster_centers.txt.")
    parser.add_argument("--feature_stats", type=Path, help="Defaults to <af2rave_dir>/feature_stats.csv.")
    parser.add_argument("--cluster_assignments", type=Path, help="Defaults to <af2rave_dir>/cluster_assignments.csv.")
    parser.add_argument("--reference_pdb", type=Path, help="Optional reference PDB to compute CA RMSD if feature_stats.csv has no RMSD.")
    parser.add_argument("--chain", help="Optional chain ID for D1/D2 and RMSD atom selection.")
    parser.add_argument("--chele", type=int, default=61, help="Conserved alphaC Glu residue number.")
    parser.add_argument("--sbridgek", type=int, default=45, help="beta3 Lys residue number.")
    parser.add_argument("--dfgasp", type=int, default=157, help="DFG Asp residue number.")
    parser.add_argument("--dfgphe", type=int, help="DFG Phe residue number. Defaults to dfgasp + 1.")
    parser.add_argument("--chelx", type=int, help="ChelX residue number. Defaults to chele + 4.")
    parser.add_argument("--residue_offset", type=int, default=0, help="Added to all residue numbers before PDB lookup.")
    parser.add_argument("--xtal_structure", type=Path, help="Optional PDB/mmCIF crystal structure to overlay in D1/D2 space.")
    parser.add_argument("--xtal_label", default="Xtal reference", help="Label for the overlaid crystal reference point.")
    parser.add_argument("--xtal_chain", help="Crystal structure chain ID. Defaults to --chain.")
    parser.add_argument("--xtal_chele", type=int, help="Crystal conserved alphaC Glu residue number.")
    parser.add_argument("--xtal_chelx", type=int, help="Crystal alphaC-Glu(+4) residue number. Overrides --xtal_chele + 4.")
    parser.add_argument("--xtal_sbridgek", type=int, help="Crystal beta3 Lys residue number.")
    parser.add_argument("--xtal_dfgasp", type=int, help="Crystal DFG Asp residue number.")
    parser.add_argument("--xtal_dfgphe", type=int, help="Crystal DFG Phe residue number. Overrides --xtal_dfgasp + 1.")
    parser.add_argument("--xtal_residue_offset", type=int, default=0, help="Added to crystal residue numbers before lookup.")
    parser.add_argument("--xtal_color", default="cyan")
    parser.add_argument("--d1_cutoff", type=float, default=11.0)
    parser.add_argument("--d2_in_cutoff", type=float, default=11.0)
    parser.add_argument("--d2_out_cutoff", type=float, default=14.0)
    parser.add_argument("--top_per_state", type=int, default=2)
    parser.add_argument("--copy_representatives", action="store_true")
    parser.add_argument("--xlim", type=parse_range, default=(2.0, 22.0))
    parser.add_argument("--ylim", type=parse_range, default=(2.0, 22.0))
    parser.add_argument("--rmsd_vmin", type=float, default=None)
    parser.add_argument("--rmsd_vmax", type=float, default=None)
    parser.add_argument("--cmap", default="hot")
    parser.add_argument("--title", default="AF2RAVE selected structures in kinase D1/D2 space")
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    args.af2rave_dir = args.af2rave_dir.expanduser().resolve()
    args.output_dir = (args.output_dir or args.af2rave_dir / "kinase_d1_d2_analysis").expanduser().resolve()
    args.structures_list = args.structures_list or args.af2rave_dir / "filtered_structures.txt"
    args.cluster_centers = args.cluster_centers or args.af2rave_dir / "cluster_centers.txt"
    args.feature_stats = args.feature_stats or args.af2rave_dir / "feature_stats.csv"
    args.cluster_assignments = args.cluster_assignments or args.af2rave_dir / "cluster_assignments.csv"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = compute_metrics(args)
    xtal_row = compute_xtal_metric(args)
    if xtal_row is not None:
        df = pd.concat([df, pd.DataFrame([xtal_row])], ignore_index=True)
    df.to_csv(args.output_dir / "kinase_d1_d2_metrics.csv", index=False)
    representatives = write_state_representatives(df, args.output_dir, args.top_per_state)
    if args.copy_representatives:
        copy_representatives(representatives, args.output_dir)
    plot_d1_d2(df, args.output_dir, args)

    print(f"Structures analyzed: {len(df)}")
    if "error" in df:
        print(f"Structures with D1/D2 errors: {int((df['error'] != '').sum())}")
    print(f"Cluster centers: {int(df['is_center'].sum()) if 'is_center' in df else 0}")
    if xtal_row is not None:
        print(f"Xtal reference: D1={xtal_row['D1_A']:.3f}, D2={xtal_row['D2_A']:.3f}, state={xtal_row['dfg_state']}")
    print("State counts among centers:")
    if not df.empty:
        print(df[df["is_center"]]["dfg_state"].value_counts().to_string())
    print(f"Metrics: {args.output_dir / 'kinase_d1_d2_metrics.csv'}")
    print(f"Representatives: {args.output_dir / 'kinase_d1_d2_state_representatives.csv'}")
    print(f"Plot: {args.output_dir / 'kinase_d1_d2_clusters.png'}")


if __name__ == "__main__":
    main()
