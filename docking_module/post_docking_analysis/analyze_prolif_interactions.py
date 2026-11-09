#!/usr/bin/env python3
"""
Run ProLIF interaction analysis on top-ranked docking poses.

The script starts from a ranked CSV produced by the ANTIOCH-VS post-docking
extractors, converts the requested top Vina PDBQT poses to SDF with Meeko's
mk_export.py, and runs ProLIF against a receptor PDB. Outputs include a
per-ligand interaction count table, a long-form ligand/residue/interaction
table, and residue interaction frequencies.
"""

from __future__ import annotations

import argparse
import copy
import json
import multiprocessing as mp
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

DEFAULT_INTERACTIONS = [
    "Hydrophobic",
    "HBDonor",
    "HBAcceptor",
    "PiStacking",
    "Anionic",
    "Cationic",
    "CationPi",
    "PiCation",
    "VdWContact",
]

SCHRODINGER_INTERACTIONS = [
    "Hydrophobic",
    "HBDonor",
    "HBAcceptor",
    "XBDonor",
    "XBAcceptor",
    "PiStacking",
    "Anionic",
    "Cationic",
    "CationPi",
    "PiCation",
    "VdWContact",
]

SCHRODINGER_MAESTRO_PARAMETERS = {
    # Schrödinger Maestro H-bond criteria:
    # H...A <= 2.8 A and D-H...A >= 120 deg. ProLIF's explicit-H H-bond
    # classes use donor-acceptor heavy-atom distance plus D-H...A angle, not
    # H...A distance, so 3.5 A is used as the closest ProLIF approximation.
    # ProLIF does not expose Maestro's H...A-B acceptor angle cutoff.
    "HBAcceptor": {"distance": 3.5, "DHA_angle": (120, 180)},
    "HBDonor": {"distance": 3.5, "DHA_angle": (120, 180)},
    # Schrödinger salt bridge cutoff: oppositely charged atoms <= 5 A.
    "Anionic": {"distance": 5.0},
    "Cationic": {"distance": 5.0},
    # Schrödinger pi criteria:
    # face-to-face centroid distance <= 4.4 A and plane angle <= 30 deg;
    # edge-to-face centroid distance <= 5.5 A and plane angle >= 30 deg;
    # pi-cation centroid distance <= 6.6 A and angle <= 30 deg.
    "PiStacking": {
        "ftf_kwargs": {"distance": 4.4, "plane_angle": (0, 30)},
        "etf_kwargs": {"distance": 5.5, "plane_angle": (30, 90)},
    },
    "CationPi": {"distance": 6.6, "angle": (0, 30)},
    "PiCation": {"distance": 6.6, "angle": (0, 30)},
}

SCHRODINGER_GLIDE_PARAMETERS = {
    **SCHRODINGER_MAESTRO_PARAMETERS,
    # Schrödinger Glide H-bond criteria:
    # H...A <= 2.5 A and D-H...A >= 90 deg. ProLIF's explicit-H H-bond
    # classes use donor-acceptor heavy-atom distance, so keep the 3.5 A
    # ProLIF distance approximation and relax the D-H...A angle.
    "HBAcceptor": {"distance": 3.5, "DHA_angle": (90, 180)},
    "HBDonor": {"distance": 3.5, "DHA_angle": (90, 180)},
}

pd = None


def require_pandas():
    global pd
    if pd is None:
        try:
            import pandas as pandas_module
        except ImportError as exc:
            raise ImportError(
                "This script requires pandas. Create/activate the post-docking "
                "environment from environment.yaml."
            ) from exc
        pd = pandas_module
    return pd


try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **_kwargs):
        return iterable


def safe_name(value: object) -> str:
    name = str(value).strip() or "ligand"
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: json_safe(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(json_safe(data), handle, indent=2, sort_keys=False)
        handle.write("\n")


def require_columns(df: pd.DataFrame, columns: list[str], csv_path: Path) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"{csv_path} is missing required columns: {', '.join(missing)}")


def read_ranked_csv(args) -> pd.DataFrame:
    require_pandas()
    csv_path = Path(args.scores_csv).resolve()
    df = pd.read_csv(csv_path)
    require_columns(df, [args.file_path_column, args.ligand_column], csv_path)

    if args.score_column and args.score_column not in df.columns:
        raise ValueError(f"{csv_path} is missing score column: {args.score_column}")

    if args.top_n:
        df = df.head(args.top_n)

    df = df.copy().reset_index(drop=True)
    df.insert(0, "Input_Rank", df.index + 1)
    return df


def check_mk_export(command: str) -> None:
    if shutil.which(command) is None:
        raise FileNotFoundError(
            f"'{command}' was not found on PATH. Activate the post-docking "
            "environment or an environment with Meeko installed."
        )


def write_first_model_pdbqt(source: Path, destination: Path) -> None:
    saw_content = False
    with source.open() as input_handle, destination.open("w") as output_handle:
        for line in input_handle:
            saw_content = True
            output_handle.write(line)
            if line.startswith("ENDMDL"):
                break

    if not saw_content:
        raise ValueError(f"Empty PDBQT file: {source}")


def convert_pdbqt_pose_to_sdf(
    pdbqt_path: Path,
    sdf_dir: Path,
    ligand_name: str,
    rank: int,
    mk_export_command: str,
    keep_pose_pdbqt: bool,
    overwrite_sdf: bool,
) -> Path:
    label = f"{rank:05d}_{safe_name(ligand_name)}"
    final_sdf = sdf_dir / f"{label}.sdf"

    if final_sdf.exists() and not overwrite_sdf:
        return final_sdf.resolve()

    tmp_root = sdf_dir / "_conversion_tmp"
    tmp_dir = tmp_root / label
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    pose_pdbqt = tmp_dir / f"{label}_pose1.pdbqt"

    try:
        write_first_model_pdbqt(pdbqt_path, pose_pdbqt)

        result = subprocess.run(
            [mk_export_command, pose_pdbqt.name],
            cwd=tmp_dir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
        if result.returncode != 0:
            message = result.stderr.strip() or f"{mk_export_command} returned {result.returncode}"
            raise RuntimeError(message)

        expected_sdf = pose_pdbqt.with_suffix(".sdf")
        if expected_sdf.exists():
            generated_sdf = expected_sdf
        else:
            generated_sdfs = sorted(tmp_dir.glob("*.sdf"), key=lambda path: path.stat().st_mtime)
            if not generated_sdfs:
                raise FileNotFoundError(f"mk_export.py did not create an SDF for {pose_pdbqt}")
            generated_sdf = generated_sdfs[-1]

        if final_sdf.exists():
            final_sdf.unlink()
        generated_sdf.rename(final_sdf)

        if keep_pose_pdbqt:
            shutil.copy2(pose_pdbqt, sdf_dir / pose_pdbqt.name)
    finally:
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir)

    return final_sdf.resolve()


def convert_pose_worker(task: dict) -> dict:
    row = task["row"]
    ligand = row[task["ligand_column"]]
    pdbqt_path = Path(row[task["file_path_column"]]).expanduser().resolve()
    expected_sdf = Path(task["sdf_dir"]) / f"{int(row['Input_Rank']):05d}_{safe_name(ligand)}.sdf"
    was_cached = expected_sdf.exists() and not task["overwrite_sdf"]

    try:
        if not pdbqt_path.exists():
            raise FileNotFoundError(pdbqt_path)

        sdf_path = convert_pdbqt_pose_to_sdf(
            pdbqt_path=pdbqt_path,
            sdf_dir=Path(task["sdf_dir"]),
            ligand_name=ligand,
            rank=int(row["Input_Rank"]),
            mk_export_command=task["mk_export_command"],
            keep_pose_pdbqt=task["keep_pose_pdbqt"],
            overwrite_sdf=task["overwrite_sdf"],
        )
        row["Pose_SDF"] = str(sdf_path)
        row["Conversion_Status"] = "cached" if was_cached else "converted"
        return {"ok": True, "row": row}
    except Exception as exc:
        row["Conversion_Error"] = str(exc)
        return {"ok": False, "row": row}


def convert_top_poses(args, ranked_df: pd.DataFrame, out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    require_pandas()
    check_mk_export(args.mk_export_command)

    sdf_dir = out_dir / args.sdf_dir
    sdf_dir.mkdir(parents=True, exist_ok=True)

    conversion_cpus = args.conversion_cpus or args.cpus or 1
    conversion_cpus = max(1, int(conversion_cpus))
    converted_rows = []
    failed_rows = []

    tasks = [
        {
            "row": row,
            "ligand_column": args.ligand_column,
            "file_path_column": args.file_path_column,
            "sdf_dir": str(sdf_dir),
            "mk_export_command": args.mk_export_command,
            "keep_pose_pdbqt": args.keep_pose_pdbqt,
            "overwrite_sdf": args.overwrite_sdf,
        }
        for row in ranked_df.to_dict(orient="records")
    ]

    print(f"Converting PDBQT pose 1 to SDF with {conversion_cpus} worker(s).")
    if conversion_cpus == 1:
        iterator = map(convert_pose_worker, tasks)
        for result in tqdm(iterator, total=len(tasks), desc="Converting PDBQT pose 1 to SDF", unit="ligand"):
            if result["ok"]:
                converted_rows.append(result["row"])
            else:
                failed_rows.append(result["row"])
                print(f"WARNING: failed to convert {result['row'].get(args.ligand_column)}: "
                      f"{result['row'].get('Conversion_Error')}")
    else:
        with mp.Pool(processes=conversion_cpus) as pool:
            iterator = pool.imap_unordered(convert_pose_worker, tasks)
            for result in tqdm(iterator, total=len(tasks), desc="Converting PDBQT pose 1 to SDF", unit="ligand"):
                if result["ok"]:
                    converted_rows.append(result["row"])
                else:
                    failed_rows.append(result["row"])
                    print(f"WARNING: failed to convert {result['row'].get(args.ligand_column)}: "
                          f"{result['row'].get('Conversion_Error')}")

    converted_df = pd.DataFrame(converted_rows)
    failed_df = pd.DataFrame(failed_rows)
    if not converted_df.empty:
        converted_df = converted_df.sort_values("Input_Rank").reset_index(drop=True)
    if not failed_df.empty:
        failed_df = failed_df.sort_values("Input_Rank").reset_index(drop=True)
    return converted_df, failed_df


def parse_residue_option(residues: list[str] | None):
    if not residues:
        return None
    if len(residues) == 1 and residues[0].lower() == "all":
        return "all"
    return residues


def deep_update(base: dict, updates: dict) -> dict:
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            deep_update(base[key], value)
        else:
            base[key] = value
    return base


def load_custom_parameters_json(value: str | None) -> dict:
    if not value:
        return {}

    path = Path(value).expanduser()
    if path.exists():
        with path.open() as handle:
            data = json.load(handle)
    else:
        data = json.loads(value)

    if not isinstance(data, dict):
        raise ValueError("--custom_parameters_json must define a JSON object.")
    return data


def custom_parameter_overrides(args) -> dict:
    overrides = load_custom_parameters_json(args.custom_parameters_json)

    hbond_overrides = {}
    if args.hbond_distance is not None:
        hbond_overrides["distance"] = args.hbond_distance
    if args.hbond_dha_angle is not None:
        hbond_overrides["DHA_angle"] = tuple(args.hbond_dha_angle)
    if hbond_overrides:
        overrides.setdefault("HBAcceptor", {}).update(hbond_overrides)
        overrides.setdefault("HBDonor", {}).update(hbond_overrides)

    if args.ionic_distance is not None:
        overrides.setdefault("Anionic", {})["distance"] = args.ionic_distance
        overrides.setdefault("Cationic", {})["distance"] = args.ionic_distance

    pi_cation_overrides = {}
    if args.pication_distance is not None:
        pi_cation_overrides["distance"] = args.pication_distance
    if args.pication_angle is not None:
        pi_cation_overrides["angle"] = tuple(args.pication_angle)
    if pi_cation_overrides:
        overrides.setdefault("CationPi", {}).update(pi_cation_overrides)
        overrides.setdefault("PiCation", {}).update(pi_cation_overrides)

    pi_stacking_overrides = overrides.setdefault("PiStacking", {})
    if args.pistacking_ftf_distance is not None:
        pi_stacking_overrides.setdefault("ftf_kwargs", {})["distance"] = args.pistacking_ftf_distance
    if args.pistacking_ftf_plane_angle is not None:
        pi_stacking_overrides.setdefault("ftf_kwargs", {})["plane_angle"] = tuple(args.pistacking_ftf_plane_angle)
    if args.pistacking_etf_distance is not None:
        pi_stacking_overrides.setdefault("etf_kwargs", {})["distance"] = args.pistacking_etf_distance
    if args.pistacking_etf_plane_angle is not None:
        pi_stacking_overrides.setdefault("etf_kwargs", {})["plane_angle"] = tuple(args.pistacking_etf_plane_angle)
    if not pi_stacking_overrides:
        overrides.pop("PiStacking", None)

    if args.hydrophobic_distance is not None:
        overrides.setdefault("Hydrophobic", {})["distance"] = args.hydrophobic_distance

    if args.vdw_tolerance is not None:
        overrides.setdefault("VdWContact", {})["tolerance"] = args.vdw_tolerance

    return overrides


def interaction_configuration(args) -> tuple[list[str] | None, dict[str, dict]]:
    if args.interactions:
        interactions = args.interactions
    elif args.interaction_preset == "prolif":
        interactions = None
    else:
        interactions = SCHRODINGER_INTERACTIONS

    if args.interaction_preset == "schrodinger_maestro":
        parameters = copy.deepcopy(SCHRODINGER_MAESTRO_PARAMETERS)
    elif args.interaction_preset == "schrodinger_glide":
        parameters = copy.deepcopy(SCHRODINGER_GLIDE_PARAMETERS)
    else:
        parameters = {}

    deep_update(parameters, custom_parameter_overrides(args))

    if interactions is not None:
        parameters = {
            name: value
            for name, value in parameters.items()
            if name in interactions
        }

    return interactions, parameters


def load_receptor(receptor_pdb: Path, guess_bonds: bool):
    try:
        import MDAnalysis as mda
        import prolif as plf
    except ImportError as exc:
        raise ImportError(
            "ProLIF analysis requires prolif and MDAnalysis. Create/activate "
            "the post-docking environment from environment.yaml."
        ) from exc

    universe = mda.Universe(str(receptor_pdb))
    if guess_bonds:
        universe.atoms.guess_bonds()
    return plf.Molecule.from_mda(universe)


def warm_rdkit_ring_info(mol) -> None:
    """
    Force RDKit ring caches to be populated before ProLIF SMARTS matching.

    Some fused aromatic systems can be missed by SMARTS using ring-size
    decorators until ring membership has been queried explicitly. Warming the
    cache makes ProLIF pi-ring matching deterministic and avoids notebook-order
    dependent behavior.
    """
    try:
        from rdkit import Chem

        Chem.GetSymmSSSR(mol)
        for atom in mol.GetAtoms():
            atom.IsInRing()
            atom.IsInRingSize(5)
            atom.IsInRingSize(6)
    except Exception:
        # Ring-cache warming is a best-effort stabilization step. If RDKit
        # rejects a molecule here, normal downstream parsing/reporting should
        # still handle the failure path.
        pass


def warm_prolif_ligand_ring_info(ligand) -> None:
    warm_rdkit_ring_info(ligand)
    try:
        for residue_id in ligand.residues:
            warm_rdkit_ring_info(ligand[residue_id])
    except Exception:
        pass


def load_ligands(converted_df: pd.DataFrame):
    try:
        import prolif as plf
        from rdkit import Chem
    except ImportError as exc:
        raise ImportError(
            "ProLIF analysis requires rdkit and prolif. Create/activate the "
            "post-docking environment from environment.yaml."
        ) from exc

    ligands = []
    metadata = []

    for row in converted_df.to_dict(orient="records"):
        sdf_path = Path(row["Pose_SDF"])
        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mol = next((item for item in supplier if item is not None), None)
        if mol is None:
            print(f"WARNING: could not read SDF molecule: {sdf_path}")
            continue

        warm_rdkit_ring_info(mol)
        ligand = plf.Molecule(mol)
        warm_prolif_ligand_ring_info(ligand)
        ligands.append(ligand)
        metadata.append(row)

    return ligands, metadata


def interaction_count(value) -> int:
    if isinstance(value, (list, tuple)):
        return len(value)
    return 1 if value else 0


def fingerprint_to_long_dataframe(fp, metadata: list[dict], args) -> pd.DataFrame:
    require_pandas()
    rows = []
    columns = [
        "Input_Rank",
        "Ligand",
        args.score_column,
        "File_Path",
        "Pose_SDF",
        "Ligand_Residue",
        "Protein_Residue",
        "Interaction",
        "Count",
    ]

    for pose_index, pair_dict in fp.ifp.items():
        if pose_index >= len(metadata):
            continue
        meta = metadata[pose_index]
        base = {
            "Input_Rank": meta.get("Input_Rank"),
            "Ligand": meta.get(args.ligand_column),
            "File_Path": meta.get(args.file_path_column),
            "Pose_SDF": meta.get("Pose_SDF"),
        }
        if args.score_column and args.score_column in meta:
            base[args.score_column] = meta.get(args.score_column)

        for (ligand_residue, protein_residue), interactions in pair_dict.items():
            for interaction_name, values in interactions.items():
                count = interaction_count(values)
                if count <= 0:
                    continue
                rows.append({
                    **base,
                    "Ligand_Residue": str(ligand_residue),
                    "Protein_Residue": str(protein_residue),
                    "Interaction": interaction_name,
                    "Count": count,
                })

    if rows:
        return pd.DataFrame(rows)
    return pd.DataFrame(columns=[column for column in columns if column])


def run_prolif(args, converted_df: pd.DataFrame):
    try:
        import prolif as plf
    except ImportError as exc:
        raise ImportError(
            "ProLIF analysis requires prolif. Create/activate the post-docking "
            "environment from environment.yaml."
        ) from exc

    receptor = load_receptor(Path(args.receptor_pdb).resolve(), args.guess_bonds)
    ligands, metadata = load_ligands(converted_df)

    if not ligands:
        raise ValueError("No readable ligand SDF files were available for ProLIF.")

    interactions, parameters = interaction_configuration(args)

    fp_kwargs = {"count": True}
    if interactions:
        fp_kwargs["interactions"] = interactions
    if parameters:
        fp_kwargs["parameters"] = parameters

    residues = parse_residue_option(args.residues)
    combined_ifp = {}
    iterator = enumerate(ligands)
    iterator = tqdm(
        iterator,
        total=len(ligands),
        desc="Running ProLIF per ligand",
        unit="ligand",
    ) if not args.no_progress else iterator
    for ligand_index, ligand in iterator:
        ligand_fp = plf.Fingerprint(**fp_kwargs)
        combined_ifp[ligand_index] = ligand_fp.generate(
            ligand,
            receptor,
            residues=residues,
            metadata=True,
        )
    fp = SimpleNamespace(ifp=combined_ifp)

    return fp, metadata


def build_ligand_counts(long_df: pd.DataFrame, metadata: list[dict], args) -> pd.DataFrame:
    require_pandas()
    metadata_columns = ["Input_Rank", args.ligand_column]
    if args.score_column:
        metadata_columns.append(args.score_column)
    metadata_columns.extend([args.file_path_column, "Pose_SDF"])

    meta_df = pd.DataFrame(metadata)
    present_columns = [column for column in metadata_columns if column in meta_df.columns]
    meta_df = meta_df[present_columns].copy()
    meta_df = meta_df.rename(columns={
        args.ligand_column: "Ligand",
        args.file_path_column: "File_Path",
    })

    if long_df.empty:
        counts_df = meta_df.copy()
        raw_interactions = DEFAULT_INTERACTIONS
        for interaction in raw_interactions:
            counts_df[interaction] = 0
    else:
        counts_df = (
            long_df
            .pivot_table(
                index=["Ligand"],
                columns="Interaction",
                values="Count",
                aggfunc="sum",
                fill_value=0,
            )
            .reset_index()
        )
        counts_df.columns.name = None
        counts_df = meta_df.merge(counts_df, on="Ligand", how="left")

    interaction_columns = [
        column for column in counts_df.columns
        if column not in {"Input_Rank", "Ligand", args.score_column, "File_Path", "Pose_SDF"}
    ]
    for column in interaction_columns:
        counts_df[column] = counts_df[column].fillna(0).astype(int)

    for interaction in DEFAULT_INTERACTIONS:
        if interaction not in counts_df.columns:
            counts_df[interaction] = 0

    hbond_cols = [
        column for column in counts_df.columns
        if column in {"HBDonor", "HBAcceptor", "XBAcceptor", "XBDonor"}
    ]
    salt_bridge_cols = [
        column for column in counts_df.columns
        if column in {"Anionic", "Cationic"}
    ]
    pi_stacking_cols = [
        column for column in counts_df.columns
        if column in {"PiStacking", "FaceToFace", "EdgeToFace"}
    ]
    pi_cation_cols = [
        column for column in counts_df.columns
        if column in {"CationPi", "PiCation"}
    ]

    counts_df["HBond_Count"] = counts_df[hbond_cols].sum(axis=1) if hbond_cols else 0
    counts_df["SaltBridge_Count"] = counts_df[salt_bridge_cols].sum(axis=1) if salt_bridge_cols else 0
    counts_df["PiStacking_Count"] = counts_df[pi_stacking_cols].sum(axis=1) if pi_stacking_cols else 0
    counts_df["PiCation_Count"] = counts_df[pi_cation_cols].sum(axis=1) if pi_cation_cols else 0

    raw_count_columns = [
        column for column in interaction_columns
        if column not in {"HBond_Count", "SaltBridge_Count", "PiStacking_Count", "PiCation_Count"}
    ]
    counts_df["Total_Interaction_Count"] = counts_df[raw_count_columns].sum(axis=1) if raw_count_columns else 0

    sort_columns = ["Input_Rank"] if "Input_Rank" in counts_df.columns else ["Ligand"]
    return counts_df.sort_values(sort_columns).reset_index(drop=True)


def build_residue_frequencies(long_df: pd.DataFrame) -> pd.DataFrame:
    require_pandas()
    if long_df.empty:
        return pd.DataFrame(columns=["Protein_Residue", "Interaction", "Ligand_Count", "Total_Count"])

    return (
        long_df
        .groupby(["Protein_Residue", "Interaction"], as_index=False)
        .agg(
            Ligand_Count=("Ligand", "nunique"),
            Total_Count=("Count", "sum"),
        )
        .sort_values(["Ligand_Count", "Total_Count"], ascending=[False, False])
        .reset_index(drop=True)
    )


def write_outputs(
    args,
    out_dir: Path,
    converted_df: pd.DataFrame,
    failed_df: pd.DataFrame,
    long_df: pd.DataFrame,
    ligand_counts: pd.DataFrame,
    residue_frequencies: pd.DataFrame,
) -> None:
    interactions, interaction_parameters = interaction_configuration(args)
    ligand_counts_parquet = out_dir / "prolif_ligand_interaction_counts.parquet"
    ligand_counts_csv = out_dir / "prolif_ligand_interaction_counts.csv"
    long_parquet = out_dir / "prolif_interactions_long.parquet"
    long_csv = out_dir / "prolif_interactions_long.csv"
    residue_csv = out_dir / "prolif_residue_interaction_frequencies.csv"
    converted_csv = out_dir / "prolif_converted_ligands.csv"
    failed_csv = out_dir / "prolif_failed_conversions.csv"
    summary_json = out_dir / args.json_out

    ligand_counts.to_parquet(ligand_counts_parquet, index=False)
    ligand_counts.to_csv(ligand_counts_csv, index=False)
    long_df.to_parquet(long_parquet, index=False)
    long_df.to_csv(long_csv, index=False)
    residue_frequencies.to_csv(residue_csv, index=False)
    converted_df.to_csv(converted_csv, index=False)
    if not failed_df.empty:
        failed_df.to_csv(failed_csv, index=False)

    summary = {
        "stage": "post_docking_analysis",
        "analysis": "prolif_interactions",
        "status": "success",
        "timestamp": datetime.now().isoformat(),
        "inputs": {
            "scores_csv": str(Path(args.scores_csv).resolve()),
            "receptor_pdb": str(Path(args.receptor_pdb).resolve()),
        },
        "configuration": {
            "top_n": args.top_n,
            "cpus": args.cpus,
            "conversion_cpus": args.conversion_cpus or args.cpus,
            "overwrite_sdf": args.overwrite_sdf,
            "interaction_preset": args.interaction_preset,
            "fingerprint_mode": "per_ligand_generate",
            "interactions": interactions or "prolif_defaults",
            "interaction_parameters": interaction_parameters,
            "residues": args.residues,
            "guess_bonds": args.guess_bonds,
        },
        "outputs": {
            "out_dir": str(out_dir),
            "sdf_dir": str(out_dir / args.sdf_dir),
            "ligand_counts_parquet": str(ligand_counts_parquet),
            "ligand_counts_csv": str(ligand_counts_csv),
            "long_interactions_parquet": str(long_parquet),
            "long_interactions_csv": str(long_csv),
            "residue_frequencies_csv": str(residue_csv),
            "converted_ligands_csv": str(converted_csv),
            "failed_conversions_csv": str(failed_csv) if not failed_df.empty else None,
            "summary_json": str(summary_json),
        },
        "metrics": {
            "ranked_rows_requested": int(args.top_n) if args.top_n else None,
            "converted_ligands": int(len(converted_df)),
            "failed_conversions": int(len(failed_df)),
            "interaction_rows": int(len(long_df)),
            "ligands_with_interactions": int(long_df["Ligand"].nunique()) if not long_df.empty else 0,
        },
    }
    write_json(summary_json, summary)

    print("\nProLIF interaction analysis complete.")
    print(f"Ligands analyzed: {len(ligand_counts)}")
    print(f"Interaction rows: {len(long_df)}")
    print(f"Ligand counts: {ligand_counts_csv}")
    print(f"Long interactions: {long_csv}")
    print(f"Residue frequencies: {residue_csv}")
    print(f"Summary JSON: {summary_json}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert top Vina PDBQT poses from a ranked CSV and calculate "
            "protein-ligand interaction fingerprints with ProLIF."
        )
    )
    parser.add_argument("--scores_csv", required=True, help="Ranked CSV from the Vina extractor.")
    parser.add_argument("--receptor_pdb", required=True, help="Prepared receptor PDB with hydrogens.")
    parser.add_argument("--out_dir", required=True, help="Output directory for ProLIF results.")
    parser.add_argument("--top_n", type=int, default=10000, help="Number of top-ranked ligands to analyze.")
    parser.add_argument(
        "-c",
        "--cpus",
        type=int,
        default=1,
        help="Workers for PDBQT-to-SDF conversion unless --conversion_cpus is set.",
    )

    parser.add_argument("--ligand_column", default="Ligand", help="Ligand name column in the CSV.")
    parser.add_argument("--score_column", default="Vina_Score", help="Docking score column in the CSV.")
    parser.add_argument("--file_path_column", default="File_Path", help="Docked PDBQT path column in the CSV.")

    parser.add_argument("--sdf_dir", default="top_ligands_sdf", help="Subfolder for converted ligand SDFs.")
    parser.add_argument("--mk_export_command", default="mk_export.py", help="Meeko PDBQT export command.")
    parser.add_argument("--keep_pose_pdbqt", action="store_true", help="Keep temporary pose-1 PDBQT files.")
    parser.add_argument(
        "--conversion_cpus",
        type=int,
        default=None,
        help="Workers for PDBQT-to-SDF conversion. Defaults to --cpus.",
    )
    parser.add_argument(
        "--overwrite_sdf",
        action="store_true",
        help="Recreate converted SDF files even when they already exist.",
    )

    parser.add_argument(
        "--interaction_preset",
        choices=["schrodinger_maestro", "schrodinger_glide", "prolif"],
        default="schrodinger_maestro",
        help=(
            "Interaction cutoff preset. Schrödinger presets map the supported "
            "Maestro/Glide cutoffs onto ProLIF parameters; 'prolif' uses raw "
            "ProLIF defaults."
        ),
    )
    parser.add_argument(
        "--custom_parameters_json",
        default=None,
        help=(
            "Optional JSON string or JSON file path with ProLIF interaction "
            "parameter overrides, e.g. '{\"HBAcceptor\": {\"distance\": 3.2}}'."
        ),
    )
    parser.add_argument("--hbond_distance", type=float, default=None, help="Override HBAcceptor/HBDonor distance.")
    parser.add_argument(
        "--hbond_dha_angle",
        type=float,
        nargs=2,
        metavar=("MIN", "MAX"),
        default=None,
        help="Override HBAcceptor/HBDonor D-H...A angle range.",
    )
    parser.add_argument("--ionic_distance", type=float, default=None, help="Override Anionic/Cationic distance.")
    parser.add_argument("--pication_distance", type=float, default=None, help="Override CationPi/PiCation distance.")
    parser.add_argument(
        "--pication_angle",
        type=float,
        nargs=2,
        metavar=("MIN", "MAX"),
        default=None,
        help="Override CationPi/PiCation angle range.",
    )
    parser.add_argument(
        "--pistacking_ftf_distance",
        type=float,
        default=None,
        help="Override face-to-face PiStacking centroid distance.",
    )
    parser.add_argument(
        "--pistacking_ftf_plane_angle",
        type=float,
        nargs=2,
        metavar=("MIN", "MAX"),
        default=None,
        help="Override face-to-face PiStacking plane angle range.",
    )
    parser.add_argument(
        "--pistacking_etf_distance",
        type=float,
        default=None,
        help="Override edge-to-face PiStacking centroid distance.",
    )
    parser.add_argument(
        "--pistacking_etf_plane_angle",
        type=float,
        nargs=2,
        metavar=("MIN", "MAX"),
        default=None,
        help="Override edge-to-face PiStacking plane angle range.",
    )
    parser.add_argument("--hydrophobic_distance", type=float, default=None, help="Override Hydrophobic distance.")
    parser.add_argument("--vdw_tolerance", type=float, default=None, help="Override VdWContact tolerance.")
    parser.add_argument(
        "--interactions",
        nargs="+",
        default=None,
        help=(
            "Optional ProLIF interaction names. Omit to use the selected "
            "interaction preset."
        ),
    )
    parser.add_argument(
        "--residues",
        nargs="+",
        default=None,
        help="Optional protein residues for ProLIF, or 'all'. Omit for nearby residues.",
    )
    parser.add_argument(
        "--guess_bonds",
        action="store_true",
        help="Ask MDAnalysis to guess receptor bonds before ProLIF conversion.",
    )
    parser.add_argument("--no_progress", action="store_true", help="Disable ProLIF progress bar.")
    parser.add_argument("--json_out", default="prolif_interaction_summary.json", help="Summary JSON filename.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    ranked_df = read_ranked_csv(args)
    converted_df, failed_df = convert_top_poses(args, ranked_df, out_dir)

    if converted_df.empty:
        raise RuntimeError("No ligands were converted to SDF; ProLIF analysis cannot run.")

    fp, metadata = run_prolif(args, converted_df)
    long_df = fingerprint_to_long_dataframe(fp, metadata, args)
    ligand_counts = build_ligand_counts(long_df, metadata, args)
    residue_frequencies = build_residue_frequencies(long_df)

    write_outputs(args, out_dir, converted_df, failed_df, long_df, ligand_counts, residue_frequencies)


if __name__ == "__main__":
    main()
