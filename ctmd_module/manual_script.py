#!/usr/bin/env python3
"""
manual_script.py — Generate hit_data_mol2.spc from a manual charge list.

Reads a two-column text file (mol2_filename <whitespace> net_charge),
verifies each mol2 file's actual net charge by summing atomic partial
charges, and writes hit_data_mol2.spc only when charges match.

Usage:
    python3 manual_script.py <input_txt> [-o output_dir]

Input format (tab or space separated, no header):
    Z12345678_pose1_docked.mol2  0
    Z87654321_pose1_docked.mol2 -1
    ...
"""
import argparse
import sys
import os
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        prog="manual_script",
        description="Generate hit_data_mol2.spc from a manual mol2+charge list, "
                    "verifying net charges against the mol2 files.",
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Two-column text file: <mol2_filename> <net_charge>",
    )
    parser.add_argument(
        "-o", "--output_dir",
        type=str,
        default=".",
        help="Directory to write hit_data_mol2.spc (default: current dir).",
    )
    return parser.parse_args()


def read_input(filepath):
    """Parse the two-column input file. Returns list of (mol2_name, expected_charge)."""
    entries = []
    with open(filepath) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                print(f"Warning: skipping malformed line {lineno}: {line!r}", file=sys.stderr)
                continue
            mol2_name = parts[0]
            try:
                charge = int(parts[1])
            except ValueError:
                # Allow float charges but round to nearest int for AMBER compatibility
                charge = round(float(parts[1]))
            entries.append((mol2_name, charge))
    return entries


def get_mol2_net_charge(mol2_path):
    """
    Sum partial charges from a mol2 file using the @<TRIPOS>ATOM block.
    In mol2 format, column 10 (0-indexed split) of each ATOM line is the charge.
    Returns the rounded integer net charge.
    """
    total = 0.0
    in_atom_block = False
    with open(mol2_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped == "@<TRIPOS>ATOM":
                in_atom_block = True
                continue
            if in_atom_block:
                # Next @<TRIPOS>... tag ends the atom block
                if stripped.startswith("@<TRIPOS>"):
                    break
                parts = line.split()
                if len(parts) >= 10:
                    try:
                        total += float(parts[9])
                    except ValueError:
                        pass
    return round(total)


def main():
    args = parse_args()
    input_path = Path(args.input_file)

    if not input_path.exists():
        print(f"Error: input file '{input_path}' not found.", file=sys.stderr)
        sys.exit(1)

    entries = read_input(input_path)
    if not entries:
        print("Error: no valid entries found in input file.", file=sys.stderr)
        sys.exit(1)

    # Determine mol2 search directory: same dir as input file by default
    mol2_dir = input_path.parent

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    errors = 0

    for mol2_name, expected_charge in entries:
        mol2_path = mol2_dir / mol2_name

        if not mol2_path.exists():
            print(f"MISSING: {mol2_name} not found in {mol2_dir}", file=sys.stderr)
            errors += 1
            continue

        actual_charge = get_mol2_net_charge(mol2_path)

        if actual_charge != expected_charge:
            print(
                f"MISMATCH: {mol2_name} — expected {expected_charge}, "
                f"got {actual_charge}",
                file=sys.stderr,
            )
            errors += 1
            continue

        results.append((mol2_name, expected_charge))
        print(f"OK: {mol2_name} charge={expected_charge}")

    if not results:
        print("\nNo entries passed verification. Nothing written.", file=sys.stderr)
        sys.exit(1)

    # Write hit_data_mol2.spc
    spc_path = output_dir / "hit_data_mol2.spc"
    with open(spc_path, "w") as f:
        for mol2_name, charge in results:
            f.write(f"{mol2_name} True {charge}\n")

    print(f"\nWrote {len(results)} entries to {spc_path}")
    if errors:
        print(f"Skipped {errors} entry/entries due to errors.", file=sys.stderr)


if __name__ == "__main__":
    main()
