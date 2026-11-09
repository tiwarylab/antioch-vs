#!/usr/bin/env python3
"""Stage collected consensus poses as one CTMD project per receptor-ligand pair.

``collect_top_consensus_poses.py`` preserves the selected docking pose and its
receptor provenance. CTMD, however, prepares one receptor with its ligand(s) in
one working directory. Consensus hits can originate from different receptors,
so they must be staged independently:

    output_dir/
      00001_active_ligand_receptor/
        rec_prepped.pdb
        molecules/
          ligand_pose.mol2
        run_ctmd_setup.sh

DOCK6 collections already contain MOL2 poses. For Vina collections, this script
converts the collected single-model PDBQT pose through Meeko to SDF and Open
Babel to MOL2. CTMD reassigns AM1-BCC charges during its own setup, so the MOL2
is used for the selected pose coordinates and molecular graph, not as the final
AMBER charge model.
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


def run(command: list[str], description: str) -> None:
    try:
        completed = subprocess.run(command, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError(f"{description}: executable not found: {command[0]}") from exc
    if completed.returncode:
        detail = (completed.stderr or completed.stdout).strip()
        raise RuntimeError(f"{description} failed: {detail or f'exit {completed.returncode}'}")


def copy_or_link(source: Path, destination: Path, mode: str) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        destination.unlink()
    if mode == "symlink":
        destination.symlink_to(source.resolve())
    else:
        shutil.copy2(source, destination)


def validate_single_mol2(path: Path) -> None:
    text = path.read_text(errors="replace")
    molecule_count = text.count("@<TRIPOS>MOLECULE")
    if molecule_count != 1:
        raise ValueError(f"Expected exactly one MOL2 molecule in {path}; found {molecule_count}.")
    if "@<TRIPOS>ATOM" not in text or "@<TRIPOS>BOND" not in text:
        raise ValueError(f"{path} is not a complete MOL2 pose.")


def materialize_ligand(source: Path, destination: Path, mk_export_cmd: str,
                        obabel_cmd: str) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    suffix = source.suffix.lower()
    if suffix == ".mol2":
        shutil.copy2(source, destination)
        validate_single_mol2(destination)
        return "copied_mol2"

    if suffix == ".sdf":
        run([obabel_cmd, "-isdf", str(source), "-omol2", "-O", str(destination)],
            "Open Babel SDF-to-MOL2 conversion")
        validate_single_mol2(destination)
        return "sdf_to_mol2"

    if suffix == ".pdbqt":
        temporary_sdf = destination.with_suffix(".sdf")
        run([mk_export_cmd, str(source), "-s", str(temporary_sdf)], "Meeko PDBQT-to-SDF export")
        try:
            run([obabel_cmd, "-isdf", str(temporary_sdf), "-omol2", "-O", str(destination)],
                "Open Babel SDF-to-MOL2 conversion")
        finally:
            temporary_sdf.unlink(missing_ok=True)
        validate_single_mol2(destination)
        return "pdbqt_to_sdf_to_mol2"

    raise ValueError(f"Unsupported collected ligand format: {source}")


def write_runner(path: Path, ctmd_module: Path) -> None:
    path.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "cd \"$(dirname \"$0\")\"\n"
        f"bash \"{ctmd_module / 'prepare_ctmd.sh'}\" rec_prepped.pdb \"${{1:-False}}\"\n"
    )
    path.chmod(0o755)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--collection_manifest", required=True, type=Path,
                        help="collection_manifest.csv produced by collect_top_consensus_poses.py")
    parser.add_argument("--out_dir", required=True, type=Path,
                        help="Directory that will contain one CTMD project per selected pair.")
    parser.add_argument("--ctmd_module", required=True, type=Path,
                        help="ANTIOCH-VS ctmd_module directory containing prepare_ctmd.sh.")
    parser.add_argument("--receptor_mode", choices=("copy", "symlink"), default="copy",
                        help="Copy (default, portable) or symlink receptor PDBs into CTMD projects.")
    parser.add_argument("--mk_export_cmd", default="mk_export.py",
                        help="Meeko exporter used only for collected Vina PDBQT poses.")
    parser.add_argument("--obabel_cmd", default="obabel",
                        help="Open Babel executable used for SDF-to-MOL2 conversion.")
    parser.add_argument("--dry_run", action="store_true")
    args = parser.parse_args()

    manifest_path = args.collection_manifest.expanduser().resolve()
    ctmd_module = args.ctmd_module.expanduser().resolve()
    if not manifest_path.is_file():
        parser.error(f"Collection manifest not found: {manifest_path}")
    if not (ctmd_module / "prepare_ctmd.sh").is_file():
        parser.error(f"CTMD module does not contain prepare_ctmd.sh: {ctmd_module}")

    with manifest_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    needed = {"Consensus_Rank", "Collected_File", "Collected_Receptor", "Status", "Receptor_Status"}
    missing = needed - set(rows[0] if rows else [])
    if missing:
        parser.error(f"Manifest lacks required columns: {sorted(missing)}")

    staged: list[dict[str, str]] = []
    failures = 0
    for row in rows:
        status = row.get("Status", "")
        receptor_status = row.get("Receptor_Status", "")
        if status not in {"EXTRACT", "BEST_MODEL", "COPIED", "COPIED_WHOLE"}:
            continue
        if receptor_status not in {"COPY", "SYMLINK"}:
            print(f"[!] rank {row['Consensus_Rank']}: receptor was not collected ({receptor_status})", file=sys.stderr)
            failures += 1
            continue

        pose = Path(row["Collected_File"])
        receptor = Path(row["Collected_Receptor"])
        pair_name = pose.parent.name
        project = args.out_dir / pair_name
        ligand_destination = project / "molecules" / "ligand_pose.mol2"
        receptor_destination = project / "rec_prepped.pdb"
        record = {
            "Consensus_Rank": row["Consensus_Rank"],
            "Ligand": row.get("Ligand", ""),
            "Class": row.get("Class", ""),
            "Best_Receptor": row.get("Best_Receptor", ""),
            "Source_Pose": str(pose),
            "Source_Receptor": str(receptor),
            "CTMD_Project": str(project),
            "CTMD_Ligand": str(ligand_destination),
            "CTMD_Receptor": str(receptor_destination),
            "Status": "",
        }
        try:
            if receptor.suffix.lower() not in {".pdb", ".ent"}:
                raise ValueError(f"CTMD receptor must be PDB, received {receptor}")
            if not pose.is_file():
                raise FileNotFoundError(f"Collected pose missing: {pose}")
            if not receptor.is_file():
                raise FileNotFoundError(f"Collected receptor missing: {receptor}")
            if not args.dry_run:
                project.mkdir(parents=True, exist_ok=True)
                materialized = materialize_ligand(pose, ligand_destination,
                                                  args.mk_export_cmd, args.obabel_cmd)
                copy_or_link(receptor, receptor_destination, args.receptor_mode)
                write_runner(project / "run_ctmd_setup.sh", ctmd_module)
            else:
                materialized = "pending"
            record["Ligand_Conversion"] = materialized
            record["Status"] = "READY" if not args.dry_run else "PENDING"
        except (OSError, RuntimeError, ValueError) as exc:
            record["Status"] = f"FAILED: {exc}"
            failures += 1
            print(f"[!] rank {row['Consensus_Rank']}: {exc}", file=sys.stderr)
        staged.append(record)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    staged_manifest = args.out_dir / "ctmd_projects_manifest.csv"
    if staged:
        with staged_manifest.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=staged[0].keys())
            writer.writeheader()
            writer.writerows(staged)
    print(f"[+] CTMD projects: {len(staged) - failures} ready, {failures} failed")
    print(f"[+] manifest: {staged_manifest}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
