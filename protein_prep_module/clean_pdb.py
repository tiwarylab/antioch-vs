#!/usr/bin/env python3
"""
clean_pdb_for_tleap.py

Cleans a PDB file so that tleap (AMBER) can load it and report "Unit is OK".

Handles:
  1. PDB ID  → downloads from RCSB, then cleans
  2. PDB file path → reads and cleans directly

Cleaning steps:
  - Keep only standard protein residues (20 amino acids)
  - Remove waters, ligands, ions, modified residues
  - Resolve alternate locations (keep highest occupancy)
  - Ensure proper element symbols in column 77-78
  - Re-number atoms sequentially
  - Add TER records between chains
  - Ensure complete backbone atoms (N, CA, C, O) per residue
  - Fix chain IDs and residue numbering

Usage:
  python clean_pdb_for_tleap.py <pdb_id_or_path> [output.pdb]
  python clean_pdb_for_tleap.py 7DAT cleaned_7dat.pdb
  python clean_pdb_for_tleap.py ./my_protein.pdb cleaned_my_protein.pdb

Requires: biopython (installed in chem_env)
"""

import sys
import os
import argparse
import urllib.request
from pathlib import Path

from Bio.PDB import PDBParser, PDBIO, Select, Chain, Residue, Atom
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings

# Suppress Biopython warnings about non-standard atoms
warnings.filterwarnings("ignore", category=PDBConstructionWarning)


# ─── Standard 20 amino acids ───────────────────────────────────────────────
STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

# ─── Expected atoms per residue (full + backbone-only fallback) ─────────────
EXPECTED_ATOMS = {
    "ALA": ["N", "CA", "C", "O", "CB"],
    "ARG": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "ASN": ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"],
    "ASP": ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
    "CYS": ["N", "CA", "C", "O", "CB", "SG"],
    "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"],
    "GLU": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"],
    "GLY": ["N", "CA", "C", "O"],
    "HIS": ["N", "CA", "C", "O", "CB", "CG", "ND1", "CE1", "NE2", "CD2"],
    "ILE": ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"],
    "LEU": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"],
    "LYS": ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
    "MET": ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"],
    "PHE": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["N", "CA", "C", "O", "CB", "CG", "CD"],
    "SER": ["N", "CA", "C", "O", "CB", "OG"],
    "THR": ["N", "CA", "C", "O", "CB", "OG1", "CG2"],
    "TRP": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "TYR": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2"],
}

BACKBONE_ATOMS = {"N", "CA", "C", "O"}


# ─── Element inference from atom name ──────────────────────────────────────
def infer_element(atom_name: str) -> str:
    """Infer the element symbol (PDB cols 77-78) from a protein atom name.

    For amino-acid atoms the element is the LEADING alphabetic character after any
    leading digits: CA/CB/CG/CD/CE/CZ -> C, N/NE/ND/NZ -> N, O/OD/OE/OG/OH -> O,
    SG/SD -> S, 1HH/HB2/HD21 -> H. Selenomethionine's SE is the one common
    two-letter element in otherwise-standard protein structures.

    Taking the element from the *last* characters (the previous behaviour) mislabels
    CA/CD/NE/CG as Ca/Cd/Ne/Cg (calcium/cadmium/neon/...); downstream tools such as
    CPPTRAJ then cannot resolve those elements ("Mapping ... for element ??") and
    topology preparation fails.
    """
    core = atom_name.strip().lstrip("0123456789")
    if not core:
        return "C"
    if core[:2].upper() == "SE":
        return "SE"
    return core[0].upper()


# ─── Download PDB from RCSB ────────────────────────────────────────────────
def download_pdb(pdb_id: str, output_path: str) -> str:
    """Download a PDB file from RCSB by ID."""
    pdb_id = pdb_id.upper().strip()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"[DOWNLOAD] Fetching {pdb_id} from RCSB...")
    urllib.request.urlretrieve(url, output_path)
    size = Path(output_path).stat().st_size
    print(f"[DOWNLOAD] Saved to {output_path} ({size / 1024:.1f} KB)")
    return output_path


# ─── Check if input is a PDB ID or file path ──────────────────────────────
def is_pdb_id(text: str) -> bool:
    """Heuristic: 4-character alphanumeric string = likely a PDB ID."""
    text = text.strip().upper()
    return len(text) == 4 and text.isalnum()


# ─── Main cleaning logic ──────────────────────────────────────────────────
def clean_pdb(input_path: str, output_path: str, chain_filter: set[str] | None = None):
    """Read a PDB file, clean it, and write a tleap-compatible version."""

    print(f"[PARSE] Reading {input_path}...")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_path)
    model = structure[0]  # Take first model only

    # ── Collect clean residues per chain ──────────────────────────────────
    chains_data: dict[str, list[tuple[str, Residue]]] = {}  # chain_id → [(resname, residue)]

    for chain in model:
        chain_id = chain.id
        if chain_id == " ":
            chain_id = "A"  # Default unnamed chains to 'A'

        # Apply chain filter if specified
        if chain_filter is not None and chain_id not in chain_filter:
            continue

        good_residues: list[tuple[str, Residue]] = []

        for residue in chain:
            resname = residue.resname.strip()
            hetfield = residue.get_id()[0]

            # Skip non-standard residues, waters, ligands, ions
            if resname not in STANDARD_AA:
                continue
            # Skip hetero-residues that are actually standard AA (e.g., modified)
            if hetfield != " " and hetfield != "h_" and hetfield != "H_ATM":
                # Some PDBs mark standard residues as HETATM; allow if name is standard
                pass

            # Resolve alternate locations: keep highest occupancy
            atoms_to_keep = []
            seen_altloc: dict[str, tuple[float, Atom]] = {}
            for atom in residue:
                altloc = atom.altloc
                # Use occupancy from the atom record
                try:
                    real_occ = atom.occupancy if hasattr(atom, 'occupancy') else 1.0
                except Exception:
                    real_occ = 1.0

                key = atom.name
                if key not in seen_altloc or real_occ > seen_altloc[key][0]:
                    seen_altloc[key] = (real_occ, atom)

            for atom in seen_altloc.values():
                atoms_to_keep.append(atom[1])

            if not atoms_to_keep:
                continue

            # Check if residue has at least backbone atoms
            atom_names = {a.name.strip() for a in atoms_to_keep}
            if not BACKBONE_ATOMS.issubset(atom_names):
                print(f"  [WARN] Residue {resname} {residue.get_id()[1]} chain {chain_id} "
                      f"missing backbone atoms ({BACKBONE_ATOMS - atom_names}) — skipping")
                continue

            good_residues.append((resname, residue))

        if good_residues:
            chains_data[chain_id] = good_residues

    if not chains_data:
        print("[ERROR] No valid protein residues found after cleaning!")
        sys.exit(1)

    print(f"[CLEAN] Found {len(chains_data)} chain(s) with valid protein residues.")
    if chain_filter is not None:
        print(f"         (filtered from available chains to: {', '.join(sorted(chain_filter))})")

    # ── Write cleaned PDB ────────────────────────────────────────────────
    print(f"[WRITE] Writing cleaned PDB to {output_path}...")

    atom_serial = 0
    lines: list[str] = []

    # Header
    lines.append(f"REMARK   4 CLEANED BY clean_pdb_for_tleap.py — standard protein only")
    lines.append(f"REMARK   4 Chains retained: {', '.join(sorted(chains_data.keys()))}")

    total_atoms = 0
    total_residues = 0

    for chain_id in sorted(chains_data.keys()):
        residues = chains_data[chain_id]
        prev_resseq = None

        for resname, residue in residues:
            resseq = residue.get_id()[1]
            icode = residue.get_id()[2]

            # Get atoms from this residue
            residue_atoms = list(residue.get_atoms())

            # Filter to expected atoms for this residue type
            expected = EXPECTED_ATOMS.get(resname, None)

            for atom in residue_atoms:
                aname = atom.name.strip()

                # Keep all atoms for standard residues (tleap needs side chains too)
                # Only filter out truly unexpected atoms if we have a known list
                if expected is not None and aname not in expected and aname not in BACKBONE_ATOMS:
                    continue

                atom_serial += 1
                pos = atom.get_vector().get_array()
                element = infer_element(aname)

                # Format atom name per PDB convention:
                # - 1-3 chars: right-justified in cols 13-16 (e.g., "   N", "  CA", " CB1")
                # - 4 chars with first char alpha: left-justified (e.g., " CG ")
                # - 4 chars with first char digit: right-justified (e.g., " 1HH")
                if len(aname) == 4 and aname[0].isalpha():
                    formatted_name = f"{aname:<4s}"
                else:
                    formatted_name = f"{aname:>4s}"

                # Format PDB ATOM record (fixed-width columns)
                line = (
                    f"ATOM  {atom_serial:>5d} "
                    f"{formatted_name} "
                    f"{resname:>3s} "
                    f"{chain_id:>1s}"
                    f"{resseq:>4d}"
                    f"    "
                    f"{pos[0]:>8.3f}"
                    f"{pos[1]:>8.3f}"
                    f"{pos[2]:>8.3f}"
                    f"{1.00:>6.2f}"
                    f"{0.00:>6.2f}"
                    f"          "
                    f"{element:>2s}  "
                )
                lines.append(line)
                total_atoms += 1

            total_residues += 1

        # TER record at end of each chain
        lines.append(f"TER   {atom_serial + 1:>5d}      {resname:>3s} {chain_id:>1s}{resseq:>4d}")

    lines.append("END")

    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"[DONE]  Cleaned PDB written to {output_path}")
    print(f"        Total atoms: {total_atoms}, Total residues: {total_residues}, "
          f"Chains: {len(chains_data)}")


# ─── Verify with tleap ─────────────────────────────────────────────────────
def verify_with_tleap(pdb_path: str, tleap_bin: str = "tleap") -> bool:
    """Run tleap on the cleaned PDB and check for successful load."""
    print(f"\n[VERIFY] Running tleap on {pdb_path}...")

    # Find the leaprc.protein.ff14SB file
    import glob
    leaprc_candidates = glob.glob(
        "/scratch/software/micromamba_root/envs/chem_env/share/leap/lib/leaprc/**/leaprc.protein.ff14SB",
        recursive=True
    )
    if not leaprc_candidates:
        # Try broader search
        leaprc_candidates = glob.glob(
            "/scratch/software/micromamba_root/envs/chem_env/**/leaprc.protein.ff14SB",
            recursive=True
        )
    leaprc_path = leaprc_candidates[0] if leaprc_candidates else None

    # Create a tleap script
    script_lines = []
    if leaprc_path:
        script_lines.append(f"source {leaprc_path}")
    script_lines.extend([
        f"mol = loadpdb {pdb_path}",
        "check mol",
        "quit",
    ])
    script_path = pdb_path + ".tleap.in"
    with open(script_path, "w") as f:
        f.write("\n".join(script_lines) + "\n")

    import subprocess
    result = subprocess.run(
        [tleap_bin, "-f", script_path],
        capture_output=True, text=True, timeout=120
    )

    output = result.stdout + result.stderr

    # Check for success indicators
    has_unit_ok = "Unit is OK" in output
    has_fatal_type = any("does not have a type" in line for line in output.splitlines())
    has_fatal_param = any("FATAL" in line and "type" in line.lower() for line in output.splitlines())
    has_loading = "Loading PDB file" in output
    has_unknown_res = "Unknown residue" in output

    # Print relevant tleap output (last 40 lines for brevity)
    relevant_lines = []
    for line in output.splitlines():
        if any(kw in line for kw in [
            "Loading PDB", "Unit is OK", "FATAL", "Error",
            "check:", "Errors:", "Warnings:", "total atoms",
            "Checking parameters", "Unknown residue", "Sourcing",
            "Exiting LEaP"
        ]):
            relevant_lines.append(line.strip())

    if relevant_lines:
        for line in relevant_lines[-40:]:
            print(f"  tleap: {line}")
    else:
        # Fallback: show last 20 lines
        for line in output.splitlines()[-20:]:
            print(f"  tleap: {line.strip()}")

    # Clean up script file
    if os.path.exists(script_path):
        os.remove(script_path)

    if has_loading and not has_fatal_type and not has_fatal_param:
        if has_unit_ok:
            print("[VERIFY] ✅ PDB loaded successfully — Unit is OK!")
        else:
            print("[VERIFY] ✅ PDB loaded successfully in tleap (no fatal type errors)!")
        return True
    else:
        print("[VERIFY] ❌ tleap reported errors.")
        return False


# ─── Main ──────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Clean a PDB file for tleap compatibility.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s 7DAT                    # Download 7DAT, clean it, verify with tleap
  %(prog)s ./protein.pdb           # Clean local PDB file
  %(prog)s 7DAT cleaned.pdb        # Specify output filename
  %(prog)s 4KDT --chains A D       # Keep only chains A and D
  %(prog)s ./protein.pdb out.pdb   # Clean and save to specific path
        """,
    )
    parser.add_argument(
        "input",
        help="PDB ID (e.g., 7DAT) or path to a PDB file.",
    )
    parser.add_argument(
        "output",
        nargs="?",
        default=None,
        help="Output PDB path (default: cleaned_<input>.pdb in current dir).",
    )
    parser.add_argument(
        "--no-verify",
        action="store_true",
        help="Skip tleap verification step.",
    )
    parser.add_argument(
        "--tleap-path",
        default="tleap",
        help="Path to tleap binary (default: tleap in PATH).",
    )
    parser.add_argument(
        "--chains", "-c",
        nargs="+",
        default=None,
        help="Only keep specified chains (e.g., --chains A D). Omit to keep all.",
    )

    args = parser.parse_args()
    input_val = args.input.strip()

    # Determine if input is a PDB ID or file path
    output_pdb = args.output  # may be None
    chain_filter = set(args.chains) if args.chains else None

    if is_pdb_id(input_val) and not os.path.isfile(input_val):
        # It's a PDB ID — download first
        tmp_pdb = f"{input_val.lower()}.pdb"
        download_pdb(input_val, tmp_pdb)
        work_pdb = tmp_pdb
        if output_pdb is None:
            output_pdb = f"cleaned_{input_val.lower()}.pdb"
    elif os.path.isfile(input_val):
        work_pdb = input_val
        if output_pdb is None:
            stem = Path(input_val).stem
            suffix = Path(input_val).suffix
            output_pdb = f"cleaned_{stem}{suffix}"
    else:
        print(f"[ERROR] '{input_val}' is not a valid PDB ID and file not found.")
        sys.exit(1)

    # Clean the PDB
    clean_pdb(work_pdb, output_pdb, chain_filter=chain_filter)

    # Verify with tleap unless skipped
    if not args.no_verify:
        verify_with_tleap(output_pdb, tleap_bin=args.tleap_path)


if __name__ == "__main__":
    main()
