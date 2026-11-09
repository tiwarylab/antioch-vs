#!/usr/bin/env python3
"""
generate_colvar_plumed.py
=========================
Generates PLUMED input files (.dat) for extracting collective-variable distances
from MD trajectories, then runs PLUMED driver to produce COLVAR files.

PBC Handling:
    PLUMED's DISTANCE action automatically applies minimum-image convention (MIC)
    when box dimensions are available from the trajectory (xtc files carry box info).
    The script verifies this at runtime by comparing raw vs PBC-corrected distances.
    For intra-protein distances in a ~80Å box, the protein typically stays within
    one box image, so raw and PBC-corrected distances are equivalent.

Usage:
    # Generate plumed.dat files only (dry run):
    python3 generate_colvar_plumed.py --dry-run

    # Test on starting PDB files only (single frame, no trajectory):
    python3 generate_colvar_plumed.py --test-pdb

    # Full run: generate + run PLUMED on XTC trajectories:
    python3 generate_colvar_plumed.py

    # Process a single PDB file directly:
    python3 generate_colvar_plumed.py --pdb path/to/structure.pdb

    # Process all PDBs in a directory:
    python3 generate_colvar_plumed.py --pdb-dir path/to/pdbs/

    # Help:
    python3 generate_colvar_plumed.py --help

Requires (set in environment before running):
    source /scratch/software/plumed-2.9.2/plumed2-2.9.2/sourceme.sh
    export PATH="/scratch/software/micromamba_root/envs/chem_env/bin:$PATH"
    export PYTHONPATH="/scratch/software/micromamba_root/envs/chem_env/lib/python3.12/site-packages:$PYTHONPATH"

Author: OpenClaw
Generated for AF2RAVE COLVAR extraction pipeline
"""

import argparse
import csv
import glob
import os
import re
import subprocess
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np


# ---------------------------------------------------------------------------
# Paths (edit these if your layout differs)
# ---------------------------------------------------------------------------
BASE_DIR=Path(".")
FEATURES_CSV = BASE_DIR / "cluster_center_structures" / "unique_features.csv"
PLUMED_BIN = "plumed"

# PLUMED driver timestep (ps) - matches Amber dt=0.002 with STRIDE handled by PLUMED
TIMESTEP_PS = 0.002

# How many trajectory frames to skip between printed COLVAR frames
# Amber ntwx=5000 writes every 5000 steps × 0.002 ps = 10 ps
# STRIDE=1 means print every frame read from trajectory
COLVAR_STRIDE = 1


# ---------------------------------------------------------------------------
# Feature name parsing
# ---------------------------------------------------------------------------
def parse_feature_name(name: str) -> tuple[dict, dict]:
    """Parse a feature name like 'LYS162A-CE-VAL173A-CA' into two atom dicts.

    Returns (atom1_info, atom2_info) where each dict has keys:
        resname, resid, chainid, atomname
    """
    # Pattern: RESNAME(3)-ATOMNAME-RESNAME(3)-ATOMNAME
    # e.g. LYS162A-CE-VAL173A-CA
    pattern = r"^([A-Z]{3})(\d+)([A-Za-z]+)-(.+)-([A-Z]{3})(\d+)([A-Za-z]+)-(.+)$"
    m = re.match(pattern, name)
    if not m:
        raise ValueError(f"Cannot parse feature name: {name}")

    atom1 = {
        "resname": m.group(1),
        "resid": int(m.group(2)),
        "chainid": m.group(3),
        "atomname": m.group(4),
    }
    atom2 = {
        "resname": m.group(5),
        "resid": int(m.group(6)),
        "chainid": m.group(7),
        "atomname": m.group(8),
    }
    return atom1, atom2


# Amber renames HIS based on protonation state; treat them all as HIS.
HIS_VARIANTS = {"HIS", "HIE", "HID", "HIP", "HSE", "HSP"}


def _resname_query(resname: str) -> str:
    """Return an MDAnalysis resname selector that matches Amber HIS variants.

    If *resname* is a HIS variant, expand to 'resname HIS HIE HID HIP HSE HSP'.
    Otherwise return the plain 'resname <name>' clause.
    """
    if resname.upper() in HIS_VARIANTS:
        # MDAnalysis accepts space-separated resnames as OR
        return "resname " + " ".join(HIS_VARIANTS)
    return f"resname {resname}"


def find_atom_index(universe: mda.Universe, atom_info: dict) -> int | None:
    """Find the 0-based atom index in an MDAnalysis Universe.

    Returns the index or None if not found.
    Handles Amber HIS protonation variants (HIE/HID/HIP/HSE/HSP → HIS).
    """
    res_query = _resname_query(atom_info["resname"])
    sel = f"{res_query} and resid {atom_info['resid']} and name {atom_info['atomname']}"
    atoms = universe.select_atoms(sel)
    if len(atoms) == 1:
        return atoms[0].ix
    return None


# ---------------------------------------------------------------------------
# Topology handling
# ---------------------------------------------------------------------------
def discover_folders(base_dir: Path) -> list[Path]:
    """Find all simulation folders (skip templates, backups, and temp dirs)."""
    skip = {"amber_MD_template", "spib_results_backup", "spib_results", "dunbrack.dat", "old_colvar_files.zip"}
    folders = []
    for d in sorted(base_dir.iterdir()):
        if not d.is_dir():
            continue
        if d.name in skip:
            continue
        if d.name.startswith("tmp_"):
            continue
        folders.append(d)
    return folders


def find_only_protein_pdb(folder: Path) -> Path | None:
    """Find the protein-only PDB that matches the stripped XTC trajectory."""
    p = folder / "only_protein.pdb"
    return p if p.exists() else None


def find_pdb(folder: Path) -> Path | None:
    """Find the starting PDB file in a folder (falls back to only_protein.pdb)."""
    pdbs = list(folder.glob("only*.pdb"))
    # Filter out template PDBs; keep the structure PDB (named after folder)
    for p in pdbs:
        if p.name.startswith(folder.name):
            return p
    if pdbs:
        return pdbs[0]
    # Last resort: any .pdb that looks like a structure file
    pdbs = sorted(folder.glob("*.pdb"), key=lambda x: x.name)
    for p in pdbs:
        if not p.name.startswith("_"):
            return p
    return None


def find_xtc(folder: Path) -> Path | None:
    """Find the XTC trajectory file."""
    xtc = folder / "06p5_Prod.xtc"
    return xtc if xtc.exists() else None


def find_nc(folder: Path) -> Path | None:
    """Find the Amber NetCDF trajectory file."""
    nc = folder / "06p5_Prod.nc"
    return nc if nc.exists() else None


def verify_pbc_handling(folder: Path, feature_names: list[str]) -> tuple[bool, str]:
    """Verify that PBC minimum-image convention is handled correctly.

    Checks whether any atom pair in the features can be split across periodic
    boundaries. If the protein stays compact within one box image, raw Euclidean
    distances are equivalent to PBC-corrected ones.

    Returns (ok, message).
    """
    xtc = find_xtc(folder)
    if xtc is None:
        return True, "No XTC trajectory — skipping PBC check"

    topo = find_only_protein_pdb(folder) or find_pdb(folder)
    if topo is None:
        return True, "No topology found — skipping PBC check"

    try:
        traj_u = mda.Universe(str(topo), str(xtc))
    except Exception:
        print("\tTrajectory:",str(xtc))
        print("\tTopology:",str(topo))
        return True, "Could not load trajectory for PBC check — skipping"

    dims = traj_u.dimensions
    if dims is None or dims[0] is None:
        return True, "No box info in trajectory — skipping PBC check"

    box = np.array([dims[0], dims[1], dims[2]])
    max_diff = 0.0
    sample_stride = max(1, traj_u.trajectory.n_frames // 50)

    for ts in traj_u.trajectory[::sample_stride]:
        coords = traj_u.atoms.positions
        for fname in feature_names[:20]:
            try:
                a1, a2 = parse_feature_name(fname)
                idx1 = find_atom_index(traj_u, a1)
                idx2 = find_atom_index(traj_u, a2)
                if idx1 is None or idx2 is None:
                    continue
                raw_d = np.linalg.norm(coords[idx1] - coords[idx2])
                diff = coords[idx1] - coords[idx2]
                diff = diff - box * np.round(diff / box)
                pbc_d = np.linalg.norm(diff)
                max_diff = max(max_diff, abs(raw_d - pbc_d))
            except Exception:
                continue

    if max_diff < 0.5:
        return True, f"PBC OK — max correction {max_diff:.2f} Å (negligible)"
    else:
        return False, f"PBC WARNING — max correction {max_diff:.2f} Å (distances may be wrong without MIC!)"


def find_rst(folder: Path) -> Path | None:
    """Find the last restart file (has full-system coordinates)."""
    rst = folder / "06p5_Prod.rst"
    return rst if rst.exists() else None


def find_prmtop(folder: Path) -> Path | None:
    """Find the Amber prmtop file."""
    p = folder / "system.prmtop"
    return p if p.exists() else None


def find_inpcrd(folder: Path) -> Path | None:
    """Find the Amber inpcrd file."""
    p = folder / "system.inpcrd"
    return p if p.exists() else None



def build_atom_index_map(folder: Path, temp_pdb_path: Path) -> mda.Universe | None:
    """Load a protein-only Universe for atom-index mapping.

    If *temp_pdb_path* is given (standalone PDB mode), load it directly.
    Otherwise prefer only_protein.pdb (matches stripped XTC), falls back to
    prmtop+inpcrd with solvent stripped, then to the starting PDB.
    Returns an MDAnalysis Universe or None on failure.
    """
    # Standalone PDB mode: load the provided file directly
    if temp_pdb_path and temp_pdb_path.exists():
        try:
            universe = mda.Universe(str(temp_pdb_path))
            print(f"  Topology source: standalone PDB {temp_pdb_path.name} ({len(universe.atoms)} atoms)")
            return universe
        except Exception:
            pass

    # 1. Prefer only_protein.pdb — it matches the stripped XTC exactly
    only_pdb = find_only_protein_pdb(folder)
    if only_pdb:
        try:
            universe = mda.Universe(str(only_pdb))
            print(f"  Topology source: only_protein.pdb ({len(universe.atoms)} atoms)")
            return universe
        except Exception:
            pass

    # 2. Fall back to prmtop+inpcrd (full system) — select protein only
    prmtop = find_prmtop(folder)
    inpcrd = find_inpcrd(folder)
    if prmtop and inpcrd:
        try:
            full = mda.Universe(str(prmtop), str(inpcrd))
            # Create a sub-Universe with protein atoms only
            protein_group = full.select_atoms("protein")
            universe = mda.Universe.agg(protein_group)
            print(f"  Topology source: prmtop+inpcrd (protein subset, {len(universe.atoms)} atoms)")
            return universe
        except Exception:
            pass

    # 3. Fall back to restart file
    rst = find_rst(folder)
    if rst:
        try:
            universe = mda.Universe(str(rst))
            print(f"  Topology source: restart file {rst.name} ({len(universe.atoms)} atoms)")
            return universe
        except Exception:
            pass

    # 4. Last resort: starting PDB
    pdb = find_pdb(folder)
    if pdb:
        try:
            universe = mda.Universe(str(pdb))
            print(f"  Topology source: protein PDB {pdb.name} ({len(universe.atoms)} atoms)")
            return universe
        except Exception:
            pass

    print(f"  [ERROR] No valid topology found for {folder.name}")
    return None


# ---------------------------------------------------------------------------
# PLUMED .dat generation
# ---------------------------------------------------------------------------
def generate_plumed_dat(
    feature_names: list[str],
    index_map: dict[str, tuple[int, int]],
    output_file: Path,
    colvar_output: str = "COLVAR",
) -> None:
    """Write a PLUMED input file with DISTANCE actions for all features."""

    lines = []
    lines.append("# Auto-generated PLUMED input for COLVAR extraction")
    lines.append(f"# {len(feature_names)} pairwise distance features")
    lines.append("")

    # DISTANCE actions - one per feature
    for i, fname in enumerate(feature_names, 1):
        if fname not in index_map:
            print(f"  [WARN] Skipping unmapped feature: {fname}")
            continue
        idx1, idx2 = index_map[fname]
        # PLUMED uses 1-based atom indexing
        lines.append(
            f"d{i}: DISTANCE ATOMS={idx1},{idx2}"
        )

    lines.append("")

    # PRINT action - output all distances as COLVAR
    # Build the ARG list
    arg_list = []
    field_list = ["time"]
    for i, fname in enumerate(feature_names, 1):
        if fname in index_map:
            arg_list.append(f"d{i}")
            field_list.append(fname)

    lines.append(
        f"PRINT ARG={','.join(arg_list)} "
        f"FILE={colvar_output} "
        f"STRIDE={COLVAR_STRIDE}"
    )

    output_file.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# PLUMED driver execution
# ---------------------------------------------------------------------------
def run_plumed_driver(
    plumed_dat: Path,
    topology_pdb: Path,
    trajectory_file: Path,
    traj_format: str,  # "xtc" or "pdb"
    work_dir: Path,
) -> bool:
    """Run PLUMED driver to extract COLVAR from trajectory.

    Returns True on success.
    """
    # Resolve to absolute paths so they work regardless of cwd.
    abs_plumed = plumed_dat.resolve()
    abs_pdb = topology_pdb.resolve()
    abs_traj = trajectory_file.resolve()

    cmd = [
        PLUMED_BIN, "driver",
        "--plumed", str(abs_plumed),
        "--pdb", str(abs_pdb),
    ]

    if traj_format == "xtc":
        cmd.extend(["--ixtc", str(abs_traj)])
    elif traj_format == "pdb":
        # For PDB trajectories, use --mf_pdb.
        # If topology and trajectory are the same file (standalone PDB mode),
        # PLUMED can use the same file for both.
        cmd.extend(["--mf_pdb", str(abs_traj)])

    cmd.extend([
        "--timestep", str(TIMESTEP_PS),
        "--trajectory-stride", "1",
    ])

    print(f"  Running: {' '.join(cmd[:8])} ...")

    result = subprocess.run(
        cmd,
        cwd=str(work_dir),
        capture_output=True,
        text=True,
        timeout=3600,  # 1 hour max per trajectory
    )

    if result.returncode != 0:
        print(f"  [ERROR] PLUMED driver failed (rc={result.returncode})")
        print(f"  stderr: {result.stderr[-1000:]}")
        return False

    # Check that COLVAR was written (use the name from the plumed.dat FILE= arg)
    # Derive it: last line of plumed.dat contains FILE=<name>
    colvar_name = "COLVAR"  # default
    try:
        dat_text = plumed_dat.read_text()
        for line in dat_text.splitlines():
            if line.startswith("PRINT"):
                import re as _re
                m = _re.search(r"FILE=(\S+)", line)
                if m:
                    colvar_name = m.group(1)
                    break
    except Exception:
        pass

    colvar = work_dir / colvar_name
    if colvar.exists():
        lines = colvar.read_text().strip().split("\n")
        print(f"  COLVAR written: {len(lines)} lines (header + {len(lines)-1} frames)")
        return True
    else:
        print(f"  [ERROR] COLVAR file not found after PLUMED run")
        return False


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------
def load_features(csv_path: Path) -> list[str]:
    """Load feature names from the unique_features.csv file."""
    features = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            features.append(row["feature_name"])
    return features


def process_folder(
    folder: Path,
    feature_names: list[str],
    universe: mda.Universe,
    dry_run: bool = False,
    test_pdb: bool = False,
    standalone_pdb: Path | None = None,
) -> dict:
    """Process a single simulation folder.

    Returns a status dict.
    """
    result = {
        "folder": folder.name,
        "status": "skipped",
        "plumed_dat": None,
        "colvar": None,
        "error": None,
    }

    print(f"\n{'=' * 60}")
    print(f"Processing: {folder.name}")
    print(f"{'=' * 60}")

    # --- Step 1: Map feature names to atom indices ---
    index_map = {}
    failures = []

    for fname in feature_names:
        try:
            atom1, atom2 = parse_feature_name(fname)
            idx1 = find_atom_index(universe, atom1)
            idx2 = find_atom_index(universe, atom2)

            if idx1 is None or idx2 is None:
                failures.append(fname)
                continue

            # PLUMED uses 1-based indexing; MDAnalysis uses 0-based
            index_map[fname] = (idx1 + 1, idx2 + 1)
        except ValueError as e:
            failures.append(f"{fname} ({e})")

    if failures:
        print(f"  [WARN] {len(failures)} features could not be mapped:")
        for f in failures[:5]:
            print(f"    - {f}")
        if len(failures) > 5:
            print(f"    ... and {len(failures)-5} more")

    if not index_map:
        result["status"] = "failed"
        result["error"] = "No features could be mapped to atom indices"
        return result

    print(f"  Mapped {len(index_map)}/{len(feature_names)} features")

    # --- Step 1.5: Verify PBC handling (skip for standalone PDBs — no box info) ---
    if not standalone_pdb:
        pbc_ok, pbc_msg = verify_pbc_handling(folder, feature_names)
        print(f"  PBC check: {pbc_msg}")
        if not pbc_ok:
            print("  [WARN] PBC correction may be needed. PLUMED DISTANCE uses minimum-image")
            print("         convention automatically when box info is available from trajectory.")
    else:
        print("  PBC check: skipped (standalone PDB — no periodic box)")

    # --- Step 2: Determine topology and trajectory files ---
    pdb = find_pdb(folder)
    xtc = find_xtc(folder)
    nc = find_nc(folder)

    # Determine which topology to use for PLUMED
    topo_pdb = None
    traj_file = None
    traj_format = None

    if standalone_pdb:
        # Standalone PDB mode: use the provided PDB as both topology and trajectory
        topo_pdb = standalone_pdb
        traj_file = standalone_pdb
        traj_format = "pdb"
    elif test_pdb:
        # Test mode: use the protein-only PDB as both topology and trajectory (single frame)
        only_pdb = find_only_protein_pdb(folder)
        if only_pdb:
            topo_pdb = only_pdb
            traj_file = only_pdb
            traj_format = "pdb"
        elif pdb:
            topo_pdb = pdb
            traj_file = pdb
            traj_format = "pdb"
        else:
            result["status"] = "failed"
            result["error"] = "No PDB file found for testing"
            return result
    else:
        # Full run mode: use protein-only topology matching the stripped XTC
        only_pdb = find_only_protein_pdb(folder)
        if xtc and only_pdb:
            topo_pdb = only_pdb
            traj_file = xtc
            traj_format = "xtc"
        elif nc and pdb:
            # Amber NetCDF trajectory — PLUMED does not natively read .nc
            result["status"] = "skipped"
            result["error"] = "Amber .nc trajectory not directly supported by PLUMED driver (needs conversion)"
            return result
        else:
            result["status"] = "skipped"
            result["error"] = "No suitable trajectory/topology pair found"
            return result

    if not standalone_pdb:
        print(f"  Topology: {topo_pdb.name}")
        print(f"  Trajectory: {traj_file.name} ({traj_format})")

    # --- Step 3: Generate plumed.dat ---
    if standalone_pdb:
        # Per-PDB plumed.dat and COLVAR names to avoid collisions
        stem = standalone_pdb.stem
        plumed_dat = standalone_pdb.parent / f"{stem}_extract_colvar.plumed.dat"
        colvar_name = f"{stem}_COLVAR"
    else:
        plumed_dat = folder / "extract_colvar.plumed.dat"
        colvar_name = "COLVAR"
    generate_plumed_dat(feature_names, index_map, plumed_dat, colvar_output=colvar_name)
    result["plumed_dat"] = str(plumed_dat)
    print(f"  Generated: {plumed_dat.name}")

    if dry_run:
        result["status"] = "dry_run"
        return result

    # --- Step 4: Run PLUMED driver ---
    if traj_format == "nc":
        # PLUMED doesn't natively read Amber .nc files.
        # We'd need to convert .nc to a format PLUMED understands.
        # For now, skip and note the limitation.
        result["status"] = "skipped"
        result["error"] = "Amber .nc trajectory not directly supported by PLUMED driver (needs conversion)"
        return result

    success = run_plumed_driver(plumed_dat, topo_pdb, traj_file, traj_format, folder)

    if success:
        result["status"] = "success"
        result["colvar"] = str(folder / colvar_name)
    else:
        result["status"] = "failed"
        result["error"] = "PLUMED driver execution failed"

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Generate PLUMED input files for COLVAR distance extraction",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Dry run - generate plumed.dat files only
  python3 %(prog)s --dry-run

  # Test on starting PDB files (single frame)
  python3 %(prog)s --test-pdb

  # Full run - generate and execute PLUMED driver on XTC trajectories
  python3 %(prog)s

  # Process specific folders only
  python3 %(prog)s --folders jak2_taf2_ifdc_rank1 rmsa04_r326_m2_s025
        """,
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Generate plumed.dat files only, do not run PLUMED",
    )
    parser.add_argument(
        "--test-pdb",
        action="store_true",
        help="Test on starting PDB files (single frame, no trajectory)",
    )
    parser.add_argument(
        "--folders",
        nargs="+",
        default=None,
        help="Process only these folder names (default: all folders)",
    )
    parser.add_argument(
        "--features-csv",
        type=Path,
        default=FEATURES_CSV,
        help=f"Path to features CSV (default: {FEATURES_CSV})",
    )
    parser.add_argument(
        "--pdb",
        type=Path,
        default=None,
        help="Process a single PDB file directly (outputs COLVAR next to it)",
    )
    parser.add_argument(
        "--pdb-dir",
        type=Path,
        default=None,
        help="Process all .pdb files in a directory (outputs COLVAR next to each PDB)",
    )

    args = parser.parse_args()

    # Mutually exclusive: --pdb / --pdb-dir vs folder-based modes
    standalone_mode = args.pdb is not None or args.pdb_dir is not None
    if standalone_mode and args.folders:
        print("[ERROR] --folders cannot be used with --pdb or --pdb-dir")
        sys.exit(1)
    if standalone_mode and args.test_pdb:
        print("[ERROR] --test-pdb cannot be used with --pdb or --pdb-dir")
        sys.exit(1)

    # --- Load features ---
    if not args.features_csv.exists():
        print(f"[ERROR] Features CSV not found: {args.features_csv}")
        sys.exit(1)

    feature_names = load_features(args.features_csv)
    print(f"Loaded {len(feature_names)} features from {args.features_csv.name}")

    # --- Discover work items ---
    if standalone_mode:
        # Standalone PDB mode: collect individual PDB files
        if args.pdb:
            if not args.pdb.exists():
                print(f"[ERROR] PDB file not found: {args.pdb}")
                sys.exit(1)
            pdb_files = [args.pdb]
        else:
            if not args.pdb_dir.exists():
                print(f"[ERROR] Directory not found: {args.pdb_dir}")
                sys.exit(1)
            pdb_files = sorted(args.pdb_dir.glob("*.pdb"))
            if not pdb_files:
                print(f"[ERROR] No .pdb files found in {args.pdb_dir}")
                sys.exit(1)
        print(f"Found {len(pdb_files)} PDB file(s) to process")
    else:
        # Folder-based mode (original behaviour)
        all_folders = discover_folders(BASE_DIR)
        if args.folders:
            folders = [BASE_DIR / f for f in args.folders if (BASE_DIR / f).exists()]
            missing = set(args.folders) - {f.name for f in folders}
            if missing:
                print(f"[WARN] Folders not found: {missing}")
        else:
            folders = all_folders

        print(f"Found {len(folders)} simulation folders")

    # --- Process each work item ---
    results = []

    if standalone_mode:
        for pdb_file in pdb_files:
            # Load Universe directly from the PDB
            try:
                universe = mda.Universe(str(pdb_file))
            except Exception as e:
                print(f"\n[ERROR] Failed to load {pdb_file.name}: {e}")
                results.append({
                    "folder": pdb_file.name,
                    "status": "failed",
                    "error": f"Failed to load PDB: {e}",
                })
                continue

            print(f"  Loaded {len(universe.atoms)} atoms, {universe.trajectory.n_frames} frame(s)")

            # Use the PDB's parent dir as work folder;
            # pass standalone_pdb so process_folder uses this file as topo+trajectory.
            result = process_folder(
                pdb_file.parent, feature_names, universe,
                dry_run=args.dry_run, test_pdb=False,
                standalone_pdb=pdb_file,
            )
            results.append(result)
    else:
        for folder in folders:
            # Load protein-only topology for atom-index mapping
            universe = build_atom_index_map(folder, None)
            if universe is None:
                results.append({
                    "folder": folder.name,
                    "status": "failed",
                    "error": "No valid topology found",
                })
                continue

            result = process_folder(
                folder, feature_names, universe,
                dry_run=args.dry_run, test_pdb=args.test_pdb,
            )
            results.append(result)

    # --- Summary ---
    print(f"\n{'=' * 60}")
    print("SUMMARY")
    print(f"{'=' * 60}")

    statuses = {}
    for r in results:
        s = r["status"]
        statuses[s] = statuses.get(s, 0) + 1

    for status, count in sorted(statuses.items()):
        print(f"  {status}: {count}")

    # List generated files
    plumed_files = [r for r in results if r.get("plumed_dat")]
    if plumed_files:
        print(f"\nGenerated PLUMED input files:")
        for r in plumed_files:
            print(f"  {r['plumed_dat']}")

    colvar_files = [r for r in results if r.get("colvar")]
    if colvar_files:
        print(f"\nGenerated COLVAR files:")
        for r in colvar_files:
            print(f"  {r['colvar']}")

    # Exit with error if any failures
    if any(r["status"] == "failed" for r in results):
        sys.exit(1)


if __name__ == "__main__":
    main()
