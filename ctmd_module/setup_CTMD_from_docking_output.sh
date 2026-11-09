#!/usr/bin/env bash
# ============================================================================
# run_ctmd_pipeline.sh — One-command CTMD setup for JAK2 test data
#
# Wraps the three-step pipeline:
#   1. Clean receptor PDBs (fix non-standard residues for tleap)
#   2. Stage CTMD projects from manifest (SDF→MOL2, copy receptors)
#   3. Run prepare_ctmd.sh on all staged projects (antechamber + tleap)
#
# Usage:
#   cd /path/to/test_files/ctmd_module_test/data_JAK2
#   bash run_ctmd_pipeline.sh [step] [options]
#
# Steps (run individually or "all"):
#   clean    — Clean receptor PDBs using clean_pdb.py
#   stage    — Stage CTMD projects from collection_manifest_fixed.csv
#   prepare  — Run prepare_ctmd.sh on all staged projects
#   check    — Print success/failure summary for all projects
#   all      — Run clean → stage → prepare → check (default)
#
# Options:
#   -j N     Max parallel jobs for prepare step (default 12)
#   -r       Re-stage from scratch (delete existing ctmd_projects/)
#   -p ID    Operate on a single project ID (e.g. 00002_decoy_...)
#   -h       Show this help
# ============================================================================

set -euo pipefail

# ─── Paths ────────────────────────────────────────────────────────────────────
DOCKFLOW_ROOT="/usr/software/dockflow/"
CTMD_SRC="$DOCKFLOW_ROOT/ctmd_module"
CLEAN_PDB="$DOCKFLOW_ROOT/protein_prep_module/clean_pdb.py"

DATA_DIR="." #"$(cd "$(dirname "$0")" && pwd)"
PROTEIN_DIR="$DATA_DIR/protein_structures"
CLEANED_DIR="$DATA_DIR/protein_structures_cleaned"
if [ "$#" -lt 1 ]; then
	echo "Select the directory with poses"
	exit 1
else
	POSES_DIR="$1"
	shift
fi
MANIFEST="$POSES_DIR/collection_manifest_fixed.csv"
PROJECTS_DIR="$DATA_DIR/ctmd_projects"

# ─── Defaults ────────────────────────────────────────────────────────────────
MAX_JOBS=12
RE_STAGE=false
SINGLE_PROJECT=""
STEP="all"

usage() {
    sed -n '3,/^$/p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

while getopts "j:rp:h" opt; do
    case $opt in
        j) MAX_JOBS=$optarg ;;
        r) RE_STAGE=true ;;
        p) SINGLE_PROJECT=$optarg ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done
shift $((OPTIND - 1))

# Allow overriding step from positional arg
if [ $# -ge 1 ]; then
    STEP="$1"
fi

# ─── Helpers ──────────────────────────────────────────────────────────────────
ensure_env() {
    if ! command -v antechamber &>/dev/null; then
        echo "[!] AMBER tools not found in PATH. Activating chem_env..."
        source /scratch/software/micromamba_root/sourceme.sh 2>/dev/null || {
            export MAMBA_ROOT_PREFIX=/scratch/software/micromamba_root
            export CONDA_PKGS_DIRS=/scratch/software/micromamba_cache/pkgs
            export XDG_CACHE_HOME=/scratch/software/micromamba_cache/xdg
            eval "$(/scratch/software/micromamba_install/bin/micromamba shell hook -s bash)"
        }
        micromamba activate chem_env
    fi
}

# ─── Step 1: Clean receptors ─────────────────────────────────────────────────
do_clean() {
    echo "=== Step 1: Cleaning receptor PDBs ==="
    mkdir -p "$CLEANED_DIR"
    local count=0
    for pdb in "$PROTEIN_DIR"/*.pdb; do
        [ -f "$pdb" ] || continue
        name=$(basename "$pdb")
        echo "  Cleaning $name ..."
        python3 "$CLEAN_PDB" "$pdb" "$CLEANED_DIR/$name"
        count=$((count + 1))
    done
    echo "  → Cleaned $count receptor(s) into $CLEANED_DIR"
}

# ─── Step 2: Stage CTMD projects ─────────────────────────────────────────────
do_stage() {
    echo "=== Step 2: Staging CTMD projects ==="
    if [ "$RE_STAGE" = true ] && [ -d "$PROJECTS_DIR" ]; then
        echo "  Removing existing ctmd_projects/ (re-stage requested)"
        rm -rf "$PROJECTS_DIR"
    fi
    mkdir -p "$PROJECTS_DIR"

    python3 "$DOCKFLOW_ROOT/docking_module/post_docking_analysis/stage_consensus_ctmd_pairs.py" \
        --collection_manifest "$MANIFEST" \
        --out_dir "$PROJECTS_DIR" \
        --ctmd_module "$CTMD_SRC" \
        --receptor_mode copy \
        --mk_export_cmd mk_export.py \
        --obabel_cmd obabel

    local total=$(find "$PROJECTS_DIR" -maxdepth 1 -mindepth 1 -type d | wc -l)
    echo "  → $total project(s) staged in $PROJECTS_DIR"
}

# ─── Step 3: Run prepare_ctmd.sh on projects ─────────────────────────────────
do_prepare() {
    echo "=== Step 3: Running prepare_ctmd.sh ==="

    local projects=()
    if [ -n "$SINGLE_PROJECT" ]; then
        dir="$PROJECTS_DIR/$SINGLE_PROJECT"
        if [ ! -d "$dir" ]; then
            echo "[!] Project $SINGLE_PROJECT not found in $PROJECTS_DIR"
            return 1
        fi
        projects=("$dir")
    else
        while IFS= read -r -d '' d; do
            projects+=("$d")
        done < <(find "$PROJECTS_DIR" -maxdepth 1 -mindepth 1 -type d -print0)
    fi

    echo "  Running on ${#projects[@]} project(s), max $MAX_JOBS parallel ..."

    local running=0
    local pids=()

    for dir in "${projects[@]}"; do
        name=$(basename "$dir")
        echo "  [START] $name"

        (
            cd "$dir"
            bash "$CTMD_SRC/prepare_ctmd.sh" rec_prepped.pdb False 2>&1 | tail -3
            echo "EXIT_$name: $?"
        ) &
        pids+=($!)
        running=$((running + 1))

        if [ $running -ge $MAX_JOBS ]; then
            wait "${pids[0]}" 2>/dev/null
            pids=("${pids[@]:1}")
            running=$((running - 1))
        fi
    done

    # Wait for remaining jobs
    for pid in "${pids[@]}"; do
        wait "$pid" 2>/dev/null
    done

    echo "  → Prepare step complete."
}

# ─── Step 4: Check results ──────────────────────────────────────────────────
do_check() {
    echo "=== Checking CTMD project status ==="
    local success=0 fail=0 failed_list=""

    for dir in "$PROJECTS_DIR"/0*; do
        [ -d "$dir" ] || continue
        name=$(basename "$dir")
        if [ -f "$dir/ligand_pose/system.prmtop" ]; then
            success=$((success + 1))
        else
            fail=$((fail + 1))
            failed_list="$failed_list\n  $name"
        fi
    done

    echo "  Total:   $((success + fail))"
    echo "  Success: $success"
    echo "  Failed:  $fail"

    if [ $fail -gt 0 ]; then
        echo ""
        echo "  Failed projects:"
        echo -e "$failed_list"
    else
        echo ""
        echo "  ✅ All projects successful!"
    fi
}

# ─── Main ─────────────────────────────────────────────────────────────────────
case $STEP in
    clean)   ensure_env; do_clean ;;
    stage)   ensure_env; do_stage ;;
    prepare) ensure_env; do_prepare ;;
    check)   do_check ;;
    all)
        ensure_env
        do_clean
        do_stage
        do_prepare
        do_check
        echo ""
        echo "=== Pipeline Complete ==="
        ;;
    *)
        echo "[!] Unknown step: $STEP"
        usage 1
        ;;
esac
