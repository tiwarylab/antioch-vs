#!/usr/bin/env python3
"""Convert GrASP atom scores (stored in PDB B-factors) into ranked pockets + docking boxes.

Pocket boundaries, ranking, and centers are produced by GrASP's OWN postprocessing functions, vendored
verbatim below (see the attribution block): atoms above a probability threshold are grouped with
average-linkage clustering, ranked by the sum of squared atom scores, and each site's center is the
volume-weighted center of its convex hull. Everything else here — reading the scored PDB, writing one
atom-selection PDB per site, and building a padded Vina/DOCK6 box per site — is ANTIOCH-VS's own
docking-support layer.

Requires: numpy, scikit-learn, scipy (for the vendored GrASP clustering).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull, Delaunay
from sklearn.cluster import AgglomerativeClustering


# ===========================================================================
# VENDORED FROM GrASP  —  copied verbatim, do not edit the function bodies.
#
#   Graph Attention Site Prediction (GrASP), file `site_metrics.py`.
#   Smith, Z.; Strobel, M.; Vani, B. P.; Tiwary, P. "Graph Attention Site
#   Prediction (GrASP): Identifying Druggable Binding Sites Using Graph Neural
#   Networks with Attention." J. Chem. Inf. Model. 2024.
#   Source:  https://github.com/tiwarylab/GrASP   (site_metrics.py)
#   License: MIT.
#
# The five functions below (sort_clusters, cluster_atoms_average, hull_center,
# get_centroid, hulls_from_clusters) are COPIED VERBATIM from GrASP's
# site_metrics.py so that our pocket boundaries, ranking, and centers match
# GrASP's published postprocessing exactly. Only the imports above are ours.
# To update, re-copy from upstream rather than editing these bodies.
# ===========================================================================

def sort_clusters(cluster_ids, probs, labels, score_type='mean'):
    """Sort clusters according to binding site scores.

    Parameters
    ----------
    cluster_ids : numpy.ndarray
        Cluster label for each atom.

    probs: numpy.ndarray
        Model predicted binding site probabilities for each atom.

    labels: numpy.ndarray
        Model predicting classes (binding/non-binding) for each atom.

    score_type: str
        Option for atom score aggregation.

    Returns
    -------
    numpy.ndarray
        Resorted clusters labels with highest scoring first.
    """

    c_probs = []
    unique_ids = np.unique(cluster_ids)

    for c_id in unique_ids[unique_ids >= 0]:
        if score_type == 'mean':
            c_prob = np.mean(probs[:, 1][labels][cluster_ids == c_id])
        elif score_type == 'sum':
            c_prob = np.sum(probs[:, 1][labels][cluster_ids == c_id])
        elif score_type == 'square':
            c_prob = np.sum((probs[:, 1][labels][cluster_ids == c_id]) ** 2)
        else:
            print('sort_clusters score_type must be mean, sum, or square.')
        c_probs.append(c_prob)

    c_order = np.argsort(c_probs)

    sorted_ids = -1 * np.ones(cluster_ids.shape)
    for new_c in range(len(c_order)):
        old_c = c_order[new_c]
        sorted_ids[cluster_ids == old_c] = new_c

    return sorted_ids


def cluster_atoms_average(all_coords, predicted_probs, threshold=.5, score_type='mean', **kwargs):
    """Cluster atoms into sites with average linkage clustering.

    Parameters
    ----------
    all_coords: numpy.ndarray
        Protein atomic coordinates.

    predicted_probs: numpy.ndarray
        Model predicted binding site probabilities for each atom.

    threshold: float
        Probability threshold to classify atoms as site/non-site.

    score_type: str
        Option for atom score aggregation.

    **kwargs
        Additional agglomerative clustering kwargs.

    Returns
    -------
    numpy.ndarray
        Coordinates of all binding site atoms.

    numpy.ndarray
        Sorted cluster ids for binding site atoms.


    numpy.ndarray
        Sorted cluster ids for all atoms with -1 representing non-site.
    """

    predicted_labels = predicted_probs[:, 1] > threshold
    if np.sum(predicted_labels) == 0:
        # No positive predictions were made with specified cutoff
        return None, None, None
    bind_coords = all_coords[predicted_labels]
    if bind_coords.shape[0] != 1:
        link_clustering = AgglomerativeClustering(linkage='average', **kwargs).fit(bind_coords)
        cluster_ids = link_clustering.labels_
        sorted_ids = sort_clusters(cluster_ids, predicted_probs, predicted_labels, score_type=score_type)
    else:
        # Under rare circumstances only one atom may be predicted as the binding pocket. In this case
        # the clustering fails so we'll just call this one atom our best 'cluster'.
        sorted_ids = np.ones(1)

    all_ids = -1 * np.ones(predicted_labels.shape)
    all_ids[predicted_labels] = sorted_ids

    return bind_coords, sorted_ids, all_ids


def hull_center(hull):
    """Compute the center of a convex hull.

    Parameters
    ----------
    hull: scipy.spatial.ConvexHull
        Convex hull to compute the center of.

    Returns
    -------
    numpy.ndarray
        Convex hull center of mass.
    """

    hull_com = np.zeros(3)
    tetras = Delaunay(hull.points[hull.vertices])

    for i in range(len(tetras.simplices)):
        tetra_verts = tetras.points[tetras.simplices][i]

        a, b, c, d = tetra_verts
        a, b, c = a - d, b - d, c - d
        tetra_vol = np.abs(np.linalg.det([a, b, c])) / 6

        tetra_com = np.mean(tetra_verts, axis=0)

        hull_com += tetra_com * tetra_vol

    hull_com = hull_com / hull.volume
    return hull_com


def get_centroid(coords):
    """Compute the centroid of a set of coordinates.

    Parameters
    ----------
    coords: numpy.ndarray
        Array of coordinates.

    Returns
    -------
    numpy.ndarray
        Centroid of input coordinates.
    """

    return np.mean(coords, axis=0)


def hulls_from_clusters(bind_coords, sorted_ids, n_sites):
    """Compute the convex hull for each binding site cluster.

    Parameters
    ----------
    bind_coords: numpy.ndarray
        Coordinates for all binding site atoms.

    sorted_ids: numpy.ndarray
        Cluster labels for all binding site atoms.

    n_sites: int
        Maximum number of binding sites.

    Returns
    -------
    list of np.ndarrays
        List of atom coordinates for each site.

    list of scipy.spatial.ConvexHulls
        Convex hull for each site.

    list of np.ndarrays
        Convex hull center for each site.
    """

    top_ids = np.unique(sorted_ids)[::-1][:n_sites]
    predicted_points_list = []
    predicted_hull_list = []
    predicted_center_list = []

    for c_id in top_ids:
        if c_id is not None:
            if c_id >= 0:
                predicted_points = bind_coords[sorted_ids == c_id]
                predicted_points_list.append(predicted_points)
                if len(predicted_points) < 4:  # You need four points to define a convex hull
                    predicted_center_list.append(get_centroid(predicted_points))
                    predicted_hull_list.append(None)
                else:
                    predicted_hull = ConvexHull(predicted_points)
                    predicted_hull_list.append(ConvexHull(predicted_points))
                    predicted_center_list.append(hull_center(predicted_hull))

    return predicted_points_list, predicted_hull_list, predicted_center_list

# ========================= end vendored GrASP code =========================


@dataclass
class Atom:
    line: str
    serial: int
    name: str
    resname: str
    chain: str
    resseq: str
    x: float
    y: float
    z: float
    score: float


@dataclass
class Pocket:
    receptor: str
    rank: int
    n_atoms: int
    n_residues: int
    score_sum: float
    score_mean: float
    score_max: float
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float
    pocket_pdb: str


def receptor_name_from_scored_pdb(pdb_path: Path) -> str:
    """Recover the input receptor name from GrASP's ``*_probs.pdb`` artifact."""
    stem = pdb_path.stem
    return stem[:-len("_probs")] if stem.endswith("_probs") else stem


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="YAML configuration file.")
    parser.add_argument("--scored-pdbs", type=Path, help="One scored PDB or a directory of scored PDBs.")
    parser.add_argument("--out-dir", type=Path, help="Output directory.")
    parser.add_argument("--score-threshold", type=float, help="Min GrASP probability for a pocket atom (GrASP default 0.3).")
    parser.add_argument("--distance-threshold", type=float, help="Average-linkage distance_threshold, in Angstrom (GrASP default 15).")
    parser.add_argument("--score-type", choices=["mean", "sum", "square"], help="Atom-score aggregation for ranking (GrASP default 'square').")
    parser.add_argument("--n-sites", type=int, help="Maximum number of ranked sites to keep per receptor.")
    parser.add_argument("--min-atoms", type=int, help="ANTIOCH-VS-only filter: discard sites with fewer atoms than this (default 1 = keep all, matching GrASP).")
    parser.add_argument("--box-padding", type=float, help="Padding around each pocket bounding box, in Angstrom.")
    parser.add_argument("--min-box-size", type=float, help="Minimum docking-box side length, in Angstrom.")
    parser.add_argument("--fixed-box-size", nargs=3, type=float, metavar=("X", "Y", "Z"), help="Optional fixed docking-box size, in Angstrom.")
    return parser.parse_args()


def load_config(path: Path | None) -> dict:
    if path is None:
        return {}
    try:
        import yaml
    except ImportError as exc:
        raise ImportError("Reading a GrASP YAML config requires pyyaml.") from exc
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def setting(args: argparse.Namespace, config: dict, arg_name: str, section: str, key: str, default):
    value = getattr(args, arg_name)
    if value is not None:
        return value
    return config.get(section, {}).get(key, default)


def discover_pdbs(path: Path) -> list[Path]:
    if path.is_file() and path.suffix.lower() == ".pdb":
        return [path]
    if path.is_dir():
        pdbs = sorted(candidate for candidate in path.rglob("*.pdb") if candidate.is_file())
        if pdbs:
            return pdbs
    raise FileNotFoundError(f"No PDB files found at: {path}")


def parse_atoms(pdb_path: Path) -> list[Atom]:
    atoms = []
    for line in pdb_path.read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            atoms.append(
                Atom(
                    line=line,
                    serial=int(line[6:11]),
                    name=line[12:16].strip(),
                    resname=line[17:20].strip(),
                    chain=line[21:22].strip(),
                    resseq=line[22:27].strip(),
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    score=float(line[60:66]),
                )
            )
        except ValueError as exc:
            raise ValueError(f"Could not parse fixed-width PDB atom record in {pdb_path}: {line}") from exc
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records found in {pdb_path}")
    return atoms


def build_pocket(receptor: str, atoms: list[Atom], rank: int, center, padding: float, min_size: float,
                 fixed_size, pocket_pdb: Path) -> Pocket:
    """Package one GrASP cluster into a Pocket record. `center` is GrASP's site center (convex-hull
    center, or centroid for <4 atoms); the box is the atoms' bounding box expanded by `padding`."""
    xs = [atom.x for atom in atoms]
    ys = [atom.y for atom in atoms]
    zs = [atom.z for atom in atoms]
    scores = [atom.score for atom in atoms]
    size = tuple(fixed_size) if fixed_size else tuple(
        max(max(values) - min(values) + 2 * padding, min_size) for values in (xs, ys, zs))
    residues = {(atom.chain, atom.resseq, atom.resname) for atom in atoms}
    return Pocket(
        receptor=receptor,
        rank=rank,
        n_atoms=len(atoms),
        n_residues=len(residues),
        score_sum=sum(scores),
        score_mean=sum(scores) / len(scores),
        score_max=max(scores),
        center_x=float(center[0]),
        center_y=float(center[1]),
        center_z=float(center[2]),
        size_x=size[0],
        size_y=size[1],
        size_z=size[2],
        pocket_pdb=str(pocket_pdb),
    )


def write_pocket_pdb(path: Path, atoms: list[Atom]) -> None:
    with path.open("w") as handle:
        for atom in atoms:
            handle.write(atom.line.rstrip() + "\n")
        handle.write("END\n")


def process_receptor(
    pdb_path: Path,
    output_dir: Path,
    score_threshold: float,
    distance_threshold: float,
    score_type: str,
    n_sites: int,
    min_atoms: int,
    box_padding: float,
    min_box_size: float,
    fixed_box_size,
) -> list[Pocket]:
    atoms = parse_atoms(pdb_path)
    coords = np.array([[atom.x, atom.y, atom.z] for atom in atoms], dtype=float)
    scores = np.array([atom.score for atom in atoms], dtype=float)
    predicted_probs = np.column_stack([1.0 - scores, scores])   # GrASP expects [p(non-site), p(site)]

    # GrASP's own clustering + ranking + hull centers (vendored functions above).
    bind_coords, sorted_ids, all_ids = cluster_atoms_average(
        coords, predicted_probs, threshold=score_threshold, score_type=score_type,
        n_clusters=None, distance_threshold=distance_threshold)

    receptor = receptor_name_from_scored_pdb(pdb_path)
    receptor_dir = output_dir / receptor
    receptor_dir.mkdir(parents=True, exist_ok=True)
    if bind_coords is None:                       # no atom above threshold
        return []

    _, _, center_list = hulls_from_clusters(bind_coords, sorted_ids, n_sites)
    top_ids = [c_id for c_id in np.unique(sorted_ids)[::-1][:n_sites] if c_id >= 0]  # same order as center_list

    pockets = []
    rank = 0
    for c_id, center in zip(top_ids, center_list):
        cluster_atoms = [atom for atom, aid in zip(atoms, all_ids) if aid == c_id]
        if len(cluster_atoms) < min_atoms:        # ANTIOCH-VS-only optional filter (default 1 = keep all)
            continue
        rank += 1
        pocket_pdb = receptor_dir / f"pocket_{rank:02d}.pdb"
        write_pocket_pdb(pocket_pdb, cluster_atoms)
        pockets.append(build_pocket(receptor, cluster_atoms, rank, center, box_padding, min_box_size,
                                    fixed_box_size, pocket_pdb))
    return pockets


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    paths = config.get("paths", {})
    scored_pdbs_value = args.scored_pdbs or paths.get("scored_pdbs")
    output_dir_value = args.out_dir or paths.get("output_dir")
    if not scored_pdbs_value:
        raise ValueError("Provide --scored-pdbs or paths.scored_pdbs in the config.")
    if not output_dir_value:
        raise ValueError("Provide --out-dir or paths.output_dir in the config.")
    scored_pdbs = Path(scored_pdbs_value).expanduser()
    output_dir = Path(output_dir_value).expanduser()

    # GrASP published postprocessing defaults: threshold 0.3, average-linkage distance_threshold 15 A,
    # 'square' site scoring, convex-hull center. min_atoms is ANTIOCH-VS's own optional filter.
    score_threshold = float(setting(args, config, "score_threshold", "pocket_clustering", "score_threshold", 0.3))
    distance_threshold = float(setting(args, config, "distance_threshold", "pocket_clustering", "distance_threshold", 15.0))
    score_type = str(setting(args, config, "score_type", "pocket_clustering", "score_type", "square"))
    n_sites = int(setting(args, config, "n_sites", "pocket_clustering", "n_sites", 100))
    min_atoms = int(setting(args, config, "min_atoms", "pocket_clustering", "min_atoms", 1))
    box_padding = float(setting(args, config, "box_padding", "docking_box", "padding", 5.0))
    min_box_size = float(setting(args, config, "min_box_size", "docking_box", "min_size", 15.0))
    fixed_box_size = tuple(args.fixed_box_size) if args.fixed_box_size else config.get("docking_box", {}).get("size")
    if fixed_box_size is not None:
        if not isinstance(fixed_box_size, (list, tuple)) or len(fixed_box_size) != 3:
            raise ValueError("docking_box.size must contain exactly three numbers or be null.")
        fixed_box_size = tuple(float(value) for value in fixed_box_size)
    if not 0.0 <= score_threshold <= 1.0:
        raise ValueError("score_threshold must be between 0 and 1.")
    if score_type not in ("mean", "sum", "square"):
        raise ValueError("score_type must be mean, sum, or square.")
    if distance_threshold <= 0 or n_sites < 1 or min_atoms < 1 or box_padding < 0 or min_box_size <= 0 \
            or (fixed_box_size and any(value <= 0 for value in fixed_box_size)):
        raise ValueError("distance_threshold, n_sites, min_atoms, and docking-box dimensions must be positive.")

    output_dir.mkdir(parents=True, exist_ok=True)
    pockets = []
    processing_rows = []
    for pdb_path in discover_pdbs(scored_pdbs):
        try:
            receptor_pockets = process_receptor(
                pdb_path, output_dir, score_threshold, distance_threshold, score_type, n_sites,
                min_atoms, box_padding, min_box_size, fixed_box_size
            )
            pockets.extend(receptor_pockets)
            processing_rows.append(
                {"scored_pdb": str(pdb_path), "status": "SUCCESS", "n_pockets": len(receptor_pockets), "error_message": ""}
            )
            print(f"{pdb_path.name}: exported {len(receptor_pockets)} ranked pocket(s).")
        except (OSError, ValueError) as exc:
            message = f"{type(exc).__name__}: {exc}"
            processing_rows.append(
                {"scored_pdb": str(pdb_path), "status": "FAILED", "n_pockets": 0, "error_message": message}
            )
            print(f"WARNING: GrASP pocket postprocessing failed for {pdb_path.name}: {message}", file=sys.stderr)

    fieldnames = list(Pocket.__dataclass_fields__)
    csv_path = output_dir / "grasp_ranked_pockets.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(asdict(pocket) for pocket in pockets)

    processing_manifest = output_dir / "grasp_postprocess_manifest.csv"
    with processing_manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["scored_pdb", "status", "n_pockets", "error_message"])
        writer.writeheader()
        writer.writerows(processing_rows)

    json_path = output_dir / "grasp_pocket_run.json"
    json_path.write_text(
        json.dumps(
            {
                "scored_pdbs": str(scored_pdbs.resolve()),
                "clustering": "GrASP site_metrics.py (vendored): average-linkage, convex-hull center",
                "score_threshold": score_threshold,
                "distance_threshold": distance_threshold,
                "score_type": score_type,
                "n_sites": n_sites,
                "min_atoms": min_atoms,
                "box_padding": box_padding,
                "min_box_size": min_box_size,
                "fixed_box_size": fixed_box_size,
                "n_pockets": len(pockets),
                "n_successful_receptors": sum(row["status"] == "SUCCESS" for row in processing_rows),
                "n_failed_receptors": sum(row["status"] == "FAILED" for row in processing_rows),
                "ranked_pockets_csv": str(csv_path.resolve()),
                "processing_manifest": str(processing_manifest.resolve()),
            },
            indent=2,
        )
        + "\n"
    )
    failures = [row for row in processing_rows if row["status"] == "FAILED"]
    if failures:
        print(f"WARNING: GrASP pocket postprocessing failed for {len(failures)} receptor(s); see {processing_manifest}.")
    print(f"Ranked pockets: {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
