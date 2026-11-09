#!/bin/bash
# ============================================================================
# clean_ctmd.sh — Remove all generated artifacts from raw_start/
#
# Usage:
#   bash clean_ctmd.sh
#
# Keeps only:
#   - molecules/          (original ligand mol2 files)
#   - prepare_ctmd.sh     (automation script)
#   - rec_prepped.pdb     (receptor)
#   - clean_ctmd.sh       (this script)
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "Cleaning generated artifacts from $SCRIPT_DIR ..."

# Per-ligand output folders
for d in Z*_pose1_docked/; do
    [ -d "$d" ] && rm -rf "$d" && echo "  Removed $d"
done

# Generated directories
[ -d "extracted_poses" ] && rm -rf extracted_poses && echo "  Removed extracted_poses/"
[ -d "MD_components" ]   && rm -rf MD_components   && echo "  Removed MD_components/"

# Generated files
for f in hit_data_mol2.spc ctmd_setup.log leap.log setup.leap _setup.sh \
         1_folder_structure.sh 2_setup.sh Z*_pose1_docked.mol2; do
    for actual in $f; do
        [ -f "$actual" ] && rm -f "$actual" && echo "  Removed $actual"
    done
done

echo ""
echo "Remaining files:"
ls -1
