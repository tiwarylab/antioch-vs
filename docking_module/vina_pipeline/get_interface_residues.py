import argparse
import csv
import os
import sys

import MDAnalysis as mda
import numpy as np


def extract_interface_residues(pdb_path, output_csv, cutoff, receptor_chains, ligand_chains):
    print(f"Loading structure: {pdb_path}")
    
    try:
        # Load Universe with just the PDB
        u = mda.Universe(pdb_path)
    except Exception as e:
        print(f"❌ Error loading PDB: {e}")
        raise RuntimeError(f"Failed to load input structure: {pdb_path}")

    # Safeguard: Re-index residues if all resids default to 1
    if np.all(u.residues.resids == 1):
        print("⚠️ Warning: Missing resid info in PDB. Automatically re-indexing residues.")
        u.add_TopologyAttr('resids', np.arange(1, len(u.residues) + 1))

    # Construct MDAnalysis selection strings dynamically
    rec_sel_str = " or ".join([f"chainID {c}" for c in receptor_chains])
    lig_sel_str = " or ".join([f"chainID {c}" for c in ligand_chains])

    print(f"Receptor Selection: {rec_sel_str}")
    print(f"Ligand Selection:   {lig_sel_str}")
    print(f"Distance Cutoff:    {cutoff} Å")

    # Find interacting residues
    # 1. Receptor residues near Ligand
    rec_interface = u.select_atoms(f"({rec_sel_str}) and around {cutoff} ({lig_sel_str})").residues
    
    # 2. Ligand residues near Receptor
    lig_interface = u.select_atoms(f"({lig_sel_str}) and around {cutoff} ({rec_sel_str})").residues

    # Write to CSV
    with open(output_csv, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Chain", "ResidueID", "ResidueName"])
        
        # Helper function to write residues and extract their specific chain ID
        def write_residues(residues):
            for res in residues:
                # Safely extract chain ID from the first atom of the residue
                chain = res.atoms[0].chainID if hasattr(res.atoms[0], 'chainID') else res.segid
                writer.writerow([chain, res.resid, res.resname])

        write_residues(rec_interface)
        write_residues(lig_interface)

    total_residues = len(rec_interface) + len(lig_interface)
    print(f"✅ Success! {total_residues} interacting residues saved to: {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract interface residues from a PDB complex based on a distance cutoff.")
    
    parser.add_argument("-i", "--input", required=True, help="Path to the input PDB file.")
    parser.add_argument("-o", "--output", default="interface_residues.csv", help="Path to save the output CSV (default: interface_residues.csv).")
    parser.add_argument("-c", "--cutoff", type=float, default=5.0, help="Distance cutoff in Angstroms (default: 5.0).")
    
    parser.add_argument("-r", "--receptor", nargs='+', required=True, 
                        help="Chain ID(s) of the receptor. Provide space-separated letters (e.g., -r A B).")
    parser.add_argument("-l", "--ligand", nargs='+', required=True, 
                        help="Chain ID(s) of the ligand. Provide space-separated letters (e.g., -l C).")
    
    args = parser.parse_args()
    
    # Ensure file exists before running
    if not os.path.exists(args.input):
        print(f"❌ Error: Input file '{args.input}' not found.")
        sys.exit(1)
    else:
        try:
            extract_interface_residues(
                pdb_path=args.input,
                output_csv=args.output,
                cutoff=args.cutoff,
                receptor_chains=args.receptor,
                ligand_chains=args.ligand
            )
        except Exception as e:
            print(f"❌ Error: {e}")
            sys.exit(1)
