#!/usr/bin/env python3
"""Extract charges from a UMF database and write hit_data_mol2.spc for CTMD."""
import json
import os
import sys
from pathlib import Path

import argparse
import subprocess


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        prog=Path(__file__).stem,
        description="Write hit_data_mol2.spc by extracting charges from a UMF database.",
    )

    parser.add_argument(
        "pose_folder",
        type=str,
        help="Folder with all ligand poses (must be in mol2 format; one ligand per mol2 file)",
    )
    parser.add_argument(
        "umf_template",
        type=str,
        help="UMF file containing all ligands used",
    )
    parser.add_argument(
        "-p", "--protein_pdb",
        type=str,
        help="Protein PDB to start MD from",
        required=True)

    # --- Optional arguments ---
    parser.add_argument(
        "-o", "--output_dir",
        type=str,
        default="./",
        help="Path to the output file (default: %(default)s).",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="Enable verbose output.",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Locate the umfpropext binary (same tree as this script)
# ---------------------------------------------------------------------------
def _find_umfpropext() -> Path:
    """Return the absolute path to the umfpropext executable."""
    # This script lives under antioch-vs/ctmd_module/
    dockflow_root = Path(__file__).resolve().parent.parent
    candidate = dockflow_root.parent / "UMDFP" / "bin" / "umfpropext"
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return candidate
    # Fallback: check PATH
    import shutil
    found = shutil.which("umfpropext")
    if found:
        return Path(found)
    raise FileNotFoundError(
        "Cannot find 'umfpropext'.  Expected at"
        f" {candidate} or on PATH."
    )


UMFPROPEXT = _find_umfpropext()


# ---------------------------------------------------------------------------
# Read the base molecule name from a mol2 file (line 2 of the header)
# ---------------------------------------------------------------------------
def _mol2_base_name(mol2_path: Path) -> str:
    """
    A mol2 header looks like:
        @<TRIPOS>MOLECULE
        Z1128837175_pose1
        ...
    We return everything before the first '_pose' token (if present).
    """
    with open(mol2_path) as fh:
        # skip the @<TRIPOS>MOLECULE tag line
        next(fh)
        name_line = next(fh).strip()

    
    try:
        name_line = name_line[: name_line.index("_pose")]
    except ValueError:
        passqq
    return name_line.strip()


# ---------------------------------------------------------------------------
# Query umfpropext for charges and return {mol_name: net_charge, ...}
# ---------------------------------------------------------------------------
def _query_charges(mol_names: list[str], umf_path: Path) -> dict[str, str]:
    """
    Pipe molecule names into umfpropext on stdin, capture JSON output,
    and parse out the 'net_charge' field for each molecule.

    Returns a dict mapping molecule name -> net_charge string.
    """
    stdin_text = "\n".join(mol_names) + "\n"

    result = subprocess.run(
        [str(UMFPROPEXT), str(umf_path), "-q", "-silent", "-json"],
        input=stdin_text,
        capture_output=True,
        text=True,
        timeout=120,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"umfpropext exited with code {result.returncode}\n"
            f"stderr: {result.stderr}"
        )

    # The JSON output looks like:
    # {
    #   "exit_code": 0,
    #   "read_molecules": [
    #     {"Z5129818600": [{"net_charge": "0"}]},
    #     ...
    #   ],
    #   "total_molecules_read": "N"
    # }
    data = json.loads(result.stdout)
    charges: dict[str, str] = {}
    for entry in data.get("read_molecules", []):
        # each entry is a single-key dict like {"Z5129818600": [{"net_charge": "0"}]}
        for mol_name, prop_list in entry.items():
            for props in prop_list:
                charges[mol_name] = props.get("net_charge", "unknown")

    return charges


def main(args: argparse.Namespace) -> int:
    """Main entry point."""

    if args.verbose:
        print(f"Input      : {args.pose_folder}")
        print(f"UMF template: {args.umf_template}")
        print(f"Output     : {args.output_dir}")

    input_path = Path(args.pose_folder)
    if not input_path.exists():
        print(f"Error: '{args.pose_folder}' does not exist.", file=sys.stderr)
        return 1

    umf_path = Path(args.umf_template)
    if not umf_path.exists():
        print(f"Error: '{args.umf_template}' does not exist.", file=sys.stderr)
        return 1

    os.makedirs(args.output_dir+"/MD_components",exist_ok=True)
    output_path = Path(args.output_dir)

    # Collect mol2 files and extract base molecule names
    mol2_files = sorted(f for f in os.listdir(input_path) if f.endswith(".mol2"))
    if args.verbose:
        print(f"{len(mol2_files)} mol2 file(s) found in input path")

    # Build a mapping: mol2_filename -> (base_name, mol2_path)
    name_map: dict[str, tuple[str, Path]] = {}
    unique_names: list[str] = []
    seen_names: set[str] = set()

    for mf in mol2_files:
        mol2_full = input_path / mf
        base_name = _mol2_base_name(mol2_full)
        name_map[mf] = (base_name, mol2_full)
        if base_name not in seen_names:
            unique_names.append(base_name)
            seen_names.add(base_name)

    if args.verbose:
        print(f"{len(unique_names)} unique molecule name(s) to query")

    # Query umfpropext for charges
    try:
        charge_map = _query_charges(unique_names, umf_path)
    except (subprocess.TimeoutExpired, RuntimeError, json.JSONDecodeError) as exc:
        print(f"Error querying umfpropext: {exc}", file=sys.stderr)
        return 1

    if args.verbose:
        for name, charge in charge_map.items():
            print(f"  {name} -> net_charge = {charge}")

    # Write the output .spc file 
    with open(output_path/"hit_data_mol2.spc", "w") as outf:
        for mf in mol2_files:
            base_name, _ = name_map[mf]
            charge = charge_map.get(base_name, "unknown")
            line = f"{mf} True {charge}\n"
            print(mf, "True", charge)
            outf.write(line)

    # Copy over key files
    exec_path = str(Path(__file__).resolve().parent)
    os.system(f"cp {exec_path}/*.sh {exec_path}/setup.leap {str(output_path)}")
    os.system(f"rm -rf {str(output_path)}/extracted_poses && cp -r {args.pose_folder} {str(output_path)}/extracted_poses")
    os.system(f"cp {str(output_path)}/hit_data_mol2.spc {str(output_path)}/extracted_poses")
    os.system(f"cp {args.protein_pdb} {str(output_path)}/MD_components/rec_prepped.pdb")
    
    # Enter folder and setup structure
    os.system(f"cd {str(output_path)} && echo y|./1_folder_structure.sh")
    os.system(f"cd {str(output_path)} && ./2_setup.sh")
    return 0


if __name__ == "__main__":
    sys.exit(main(parse_args()))
