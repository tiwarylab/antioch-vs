#!/usr/bin/env python3
"""Create collection_manifest_fixed.csv for dock6 top1pct pairs.

The original manifest has paths pointing to a previous machine's project mount.
This script rewrites:
  - Collected_File → local path to <pose_name>_pose.mol2 in poses/
  - Collected_Receptor → cleaned PDB in protein_structures_cleaned/
  - Adds Exported_SDF and SDF_Status columns (N/A for dock6, already MOL2)
"""
import csv
from pathlib import Path
import sys,os

pose_folder="dock6_top1pct_pairs/"
if len(sys.argv)>1:
    pose_folder=sys.argv[1]
BASE = Path(os.path.abspath("."))
POSES_DIR = BASE / pose_folder / "poses"
CLEANED_DIR = BASE / "protein_structures_cleaned"

src = BASE / pose_folder / "collection_manifest.csv"
dst = BASE / pose_folder / "collection_manifest_fixed.csv"

with open(src, newline="") as f:
    reader = csv.DictReader(f)
    rows = list(reader)

# New columns to match vina fixed manifest format
new_fields = [
    "Consensus_Rank", "Ligand", "Raw_Ligand", "Class",
    "Consensus_RRF_Score", "Best_Receptor", "Best_Raw_Score",
    "Within_Receptor_Rank", "Top_Support_Votes", "Source_File",
    "Collected_File", "Model_Affinity", "Status",
    "Receptor_File", "Collected_Receptor", "Receptor_Status",
    "Exported_SDF", "SDF_Status",
]

fixed_rows = []
for row in rows:
    rank = row["Consensus_Rank"]
    ligand = row["Ligand"]
    best_rec = row["Best_Receptor"]

    # Extract pair folder name from the original Collected_File path
    # e.g. .../poses/00001_decoy_ZINC00820863_rmsa08_r251_m1_s080/ZINC00820863_pose.mol2
    collected_file = row["Collected_File"]
    pair_name = Path(collected_file).parent.name

    # Local pose file
    local_pose = POSES_DIR / pair_name / (ligand + "_pose.mol2")

    # Cleaned receptor PDB
    rec_pdb = CLEANED_DIR / (best_rec + ".pdb")

    fixed = {
        "Consensus_Rank": rank,
        "Ligand": ligand,
        "Raw_Ligand": row["Raw_Ligand"],
        "Class": row["Class"],
        "Consensus_RRF_Score": row["Consensus_RRF_Score"],
        "Best_Receptor": best_rec,
        "Best_Raw_Score": row["Best_Raw_Score"],
        "Within_Receptor_Rank": row["Within_Receptor_Rank"],
        "Top_Support_Votes": row["Top_Support_Votes"],
        "Source_File": row["Source_File"],
        "Collected_File": str(local_pose),
        "Model_Affinity": row["Best_Raw_Score"],  # Use raw score as affinity
        "Status": row["Status"],  # EXTRACT for dock6
        "Receptor_File": str(rec_pdb),
        "Collected_Receptor": str(rec_pdb),
        "Receptor_Status": "COPY",
        "Exported_SDF": "",
        "SDF_Status": "N/A",
    }
    fixed_rows.append(fixed)

with open(dst, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=new_fields)
    writer.writeheader()
    writer.writerows(fixed_rows)

# Verify files exist
missing_poses = 0
missing_recs = 0
for r in fixed_rows:
    if not Path(r["Collected_File"]).is_file():
        missing_poses += 1
        print(f"[!] Missing pose: {r['Collected_File']}")
    if not Path(r["Collected_Receptor"]).is_file():
        missing_recs += 1
        print(f"[!] Missing receptor: {r['Collected_Receptor']}")

print(f"\nFixed manifest written to: {dst}")
print(f"Total rows: {len(fixed_rows)}")
print(f"Missing poses: {missing_poses}")
print(f"Missing receptors: {missing_recs}")
