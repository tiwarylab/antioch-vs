#!/usr/bin/env python3
"""Relax receptor PDBs with AlphaFold's restrained AmberRelaxation protocol.

This is a local energy minimization intended to resolve clashes before docking
receptor preparation. It is not molecular dynamics or full protein refinement.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

# Normalized residue names from the AlphaFold Amber-relax notebook workflow.
MODRES = {
    "MSE": "MET", "MLY": "LYS", "FME": "MET", "HYP": "PRO", "TPO": "THR",
    "CSO": "CYS", "SEP": "SER", "M3L": "LYS", "HSK": "HIS", "SAC": "SER",
    "PCA": "GLU", "DAL": "ALA", "CME": "CYS", "CSD": "CYS", "OCS": "CYS",
    "DPR": "PRO", "B3K": "LYS", "ALY": "LYS", "YCM": "CYS", "MLZ": "LYS",
    "4BF": "TYR", "KCX": "LYS", "B3E": "GLU", "B3D": "ASP", "HZP": "PRO",
    "CSX": "CYS", "BAL": "ALA", "HIC": "HIS", "DBZ": "ALA", "DCY": "CYS",
    "DVA": "VAL", "NLE": "LEU", "SMC": "CYS", "AGM": "ARG", "B3A": "ALA",
    "DAS": "ASP", "DLY": "LYS", "DSN": "SER", "DTH": "THR", "GL3": "GLY",
    "HY3": "PRO", "LLP": "LYS", "MGN": "GLN", "MHS": "HIS", "TRQ": "TRP",
    "B3Y": "TYR", "PHI": "PHE", "PTR": "TYR", "TYS": "TYR", "IAS": "ASP",
    "GPL": "LYS", "KYN": "TRP", "SEC": "CYS",
}


def import_alphafold():
    try:
        from alphafold.common import protein, residue_constants
        from alphafold.relax import relax
    except ImportError as exc:
        raise SystemExit(
            "AlphaFold relaxation modules are unavailable. Run this script with "
            "a local-ColabFold/AlphaFold Python environment that includes OpenMM."
        ) from exc
    return protein, relax, residue_constants


def pdb_to_string(pdb_file: Path, residue_constants) -> str:
    """Normalize PDB records before AlphaFold's protein parser reads them."""
    lines = []
    seen = set()
    modres = dict(MODRES)
    model = 1
    for raw_line in pdb_file.read_bytes().splitlines():
        line = raw_line.decode("utf-8", "ignore").rstrip()
        if line.startswith("MODEL"):
            model = int(line[5:])
        if line.startswith("MODRES"):
            source_residue = line[12:15]
            target_residue = line[24:27]
            if source_residue not in modres and target_residue in residue_constants.restype_3to1:
                modres[source_residue] = target_residue
        if line.startswith("HETATM") and line[17:20] in modres:
            line = "ATOM  " + line[6:17] + modres[line[17:20]] + line[20:]
        if line.startswith("ATOM"):
            atom_name = line[12:16].strip()
            residue_name = line[17:20]
            residue_id = line[22:27].strip()
            chain = line[21:22]
            if residue_id and residue_id[-1].isalpha():
                residue_id = residue_id[:-1]
                line = line[:26] + " " + line[27:]
            key = f"{model}_{chain}_{residue_id}_{residue_name}_{atom_name}"
            if key not in seen:
                lines.append(line)
                seen.add(key)
        elif line.startswith(("MODEL", "TER", "ENDMDL")):
            lines.append(line)
    return "\n".join(lines)


def relax_one(input_pdb: Path, output_pdb: Path, max_iterations: int,
              tolerance: float, stiffness: float, max_outer_iterations: int,
              use_gpu: bool) -> dict:
    protein, relax, residue_constants = import_alphafold()
    protein_obj = protein.from_pdb_string(pdb_to_string(input_pdb, residue_constants))
    amber_relaxer = relax.AmberRelaxation(
        max_iterations=max_iterations,
        tolerance=tolerance,
        stiffness=stiffness,
        exclude_residues=[],
        max_outer_iterations=max_outer_iterations,
        use_gpu=use_gpu,
    )
    relaxed_pdb, debug_data, violations = amber_relaxer.process(prot=protein_obj)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    output_pdb.write_text(relaxed_pdb)
    return {"debug": debug_data, "violations": violations}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--input-pdb", type=Path, help="One receptor PDB to relax.")
    source.add_argument("--input-dir", type=Path, help="Directory of receptor PDBs to relax.")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--glob", default="*.pdb", help="Pattern used with --input-dir (default: *.pdb).")
    parser.add_argument("--max-iterations", type=int, default=2000)
    parser.add_argument("--tolerance", type=float, default=2.39)
    parser.add_argument("--stiffness", type=float, default=10.0,
                        help="Heavy-atom restraint strength in kcal/mol/A^2 (default: 10).")
    parser.add_argument("--max-outer-iterations", type=int, default=3)
    parser.add_argument("--use-gpu", action="store_true", help="Use an OpenMM GPU platform when available.")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if args.max_iterations < 0:
        parser.error("--max-iterations cannot be negative")
    if args.tolerance <= 0:
        parser.error("--tolerance must be positive")
    if args.stiffness < 0:
        parser.error("--stiffness cannot be negative")
    if args.max_outer_iterations < 1:
        parser.error("--max-outer-iterations must be at least 1")

    # AlphaFold's Amber code uses JAX even for CPU minimization. Do not let a
    # CUDA-enabled JAX installation probe a nonexistent GPU on a CPU node.
    if not args.use_gpu:
        os.environ["JAX_PLATFORMS"] = "cpu"

    inputs = [args.input_pdb] if args.input_pdb else sorted(args.input_dir.glob(args.glob))
    if not inputs:
        raise SystemExit("No input PDBs found.")

    for input_pdb in inputs:
        if not input_pdb.is_file():
            raise SystemExit(f"Input PDB does not exist: {input_pdb}")
        output_pdb = args.output_dir / input_pdb.name
        if output_pdb.exists() and not args.overwrite:
            print(f"Skipping existing: {output_pdb}")
            continue
        print(f"Amber-relaxing: {input_pdb.name}")
        metadata = relax_one(
            input_pdb,
            output_pdb,
            args.max_iterations,
            args.tolerance,
            args.stiffness,
            args.max_outer_iterations,
            args.use_gpu,
        )
        report_path = output_pdb.with_suffix(".amber_relax.json")
        report_path.write_text(json.dumps(metadata, indent=2, default=str) + "\n")
        print(f"  wrote: {output_pdb}")
        print(f"  report: {report_path}")


if __name__ == "__main__":
    main()
