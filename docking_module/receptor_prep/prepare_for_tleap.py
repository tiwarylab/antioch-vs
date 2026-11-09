#!/usr/bin/env python3
"""Prepare reviewed AMBER-state PDB files for LEaP with CPPTRAJ prepareforleap.

This companion to prepare_amber_protonation.py preserves explicit AMBER
residue names (ASH/GLH/HID/HIE/HIP/CYX) and detects final disulfide
connectivity with CPPTRAJ.
Optionally, it runs Reduce to add and optimize hydrogens before CPPTRAJ uses
their placement to assign explicit histidine names. It writes a prepared PDB,
a disulfide report, and a minimal LEaP modifications file. When no required
disulfides are supplied for a group of AF2-style PDBs, it selects the model
with the highest residue-mean pLDDT as a bookkeeping reference and reports
models with a different detected disulfide set.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from pathlib import Path


def input_pdbs(input_path: Path, recursive: bool, pattern: str) -> list[Path]:
    if input_path.is_file():
        if input_path.suffix.lower() != ".pdb":
            raise ValueError(f"Input file must be a PDB, got {input_path}.")
        return [input_path]
    if input_path.is_dir():
        glob_pattern = f"**/{pattern}" if recursive else pattern
        pdbs = sorted(path for path in input_path.glob(glob_pattern) if path.is_file())
        if pdbs:
            return pdbs
        raise FileNotFoundError(f"No files matching {pattern!r} found in {input_path}.")
    raise FileNotFoundError(f"Input path does not exist: {input_path}")


def cpptraj_input(
    input_pdb: Path,
    output_pdb: Path,
    cpptraj_leap_file: Path,
    leap_unit: str,
    keep_water: bool,
    keep_hydrogens: bool,
    detect_histidines: bool,
    detect_disulfides: bool,
    disulfide_sg_cutoff: float,
    altloc: str,
    prepare_sugars: bool,
) -> str:
    command = (
        f"prepareforleap crdset CRD name Final pdbout {output_pdb} "
        f"out {cpptraj_leap_file} leapunitname {leap_unit} keepaltloc {altloc}"
    )
    if not keep_water:
        command += " nowat"
    if not keep_hydrogens:
        command += " noh"
    if not detect_histidines:
        command += " nohisdetect"
    if not detect_disulfides:
        command += " nodisulfides"
    else:
        command += f" disulfidecut {disulfide_sg_cutoff:.3f}"
    if not prepare_sugars:
        command += " nosugars"
    return "\n".join(
        (
            f"parm {input_pdb}",
            f"loadcrd {input_pdb} name CRD",
            command,
            "go",
            "",
        )
    )


REQUIRED_DISULFIDE_PATTERN = re.compile(r"^(?P<chain_a>[^:\s]+):(?P<resi_a>-?\d+)-(?P<chain_b>[^:\s]+):(?P<resi_b>-?\d+)$")
LEAP_BOND_PATTERN = re.compile(r"^bond\s+\S+\.(?P<resi_a>-?\d+)\.SG\s+\S+\.(?P<resi_b>-?\d+)\.SG\s*$")


def parse_required_disulfide(value: str) -> tuple[tuple[str, int], tuple[str, int]]:
    match = REQUIRED_DISULFIDE_PATTERN.fullmatch(value)
    if not match:
        raise ValueError(
            f"Invalid required disulfide {value!r}. Use CHAIN:RESIDUE-CHAIN:RESIDUE, "
            "for example A:75-A:160."
        )
    return (
        (match.group("chain_a"), int(match.group("resi_a"))),
        (match.group("chain_b"), int(match.group("resi_b"))),
    )


def sg_atoms_by_residue(pdb_path: Path) -> dict[int, list[dict[str, object]]]:
    atoms: dict[int, list[dict[str, object]]] = {}
    for line in pdb_path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or line[12:16].strip() != "SG":
            continue
        try:
            residue = int(line[22:26])
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        atoms.setdefault(residue, []).append(
            {"chain": line[21].strip() or "_", "residue": residue, "resname": line[17:20].strip(), "xyz": xyz}
        )
    return atoms


def detected_disulfides(cpptraj_bond_file: Path, prepared_pdb: Path) -> list[dict[str, object]]:
    """Map CPPTRAJ's LEaP bonds back to chain-aware PDB residues."""
    sg_atoms = sg_atoms_by_residue(prepared_pdb)
    detected: list[dict[str, object]] = []
    for line in cpptraj_bond_file.read_text(errors="replace").splitlines():
        match = LEAP_BOND_PATTERN.match(line)
        if not match:
            continue
        first_residue = int(match.group("resi_a"))
        second_residue = int(match.group("resi_b"))
        first_atoms = sg_atoms.get(first_residue, [])
        second_atoms = sg_atoms.get(second_residue, [])
        if len(first_atoms) != 1 or len(second_atoms) != 1:
            raise ValueError(
                "Could not map a CPPTRAJ disulfide bond unambiguously back to PDB chains: "
                f"{line!r}. This helper currently requires unique residue numbers per input PDB."
            )
        first, second = first_atoms[0], second_atoms[0]
        first_xyz = first["xyz"]
        second_xyz = second["xyz"]
        assert isinstance(first_xyz, tuple) and isinstance(second_xyz, tuple)
        distance = sum((a - b) ** 2 for a, b in zip(first_xyz, second_xyz)) ** 0.5
        detected.append(
            {
                "chain_1": first["chain"],
                "residue_1": first["residue"],
                "resname_1": first["resname"],
                "chain_2": second["chain"],
                "residue_2": second["residue"],
                "resname_2": second["resname"],
                "sg_distance_angstrom": distance,
                "cpptraj_bond": line,
            }
        )
    return detected


def pair_key(first: tuple[str, int], second: tuple[str, int]) -> frozenset[tuple[str, int]]:
    return frozenset((first, second))


def pair_label(pair: frozenset[tuple[str, int]]) -> str:
    first, second = sorted(pair)
    return f"{first[0]}:{first[1]}-{second[0]}:{second[1]}"


def residue_mean_plddt(pdb_path: Path) -> float | None:
    """Return the mean of per-residue B-factor means from an AF2-style PDB."""
    residue_scores: dict[tuple[str, str, str], list[float]] = {}
    for line in pdb_path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            score = float(line[60:66])
        except ValueError:
            continue
        residue_key = (line[21].strip() or "_", line[22:26].strip(), line[26].strip())
        residue_scores.setdefault(residue_key, []).append(score)
    if not residue_scores:
        return None
    return sum(sum(scores) / len(scores) for scores in residue_scores.values()) / len(residue_scores)


def detected_pair_keys(report_path: Path) -> set[frozenset[tuple[str, int]]]:
    with report_path.open(newline="") as handle:
        return {
            pair_key(
                (row["chain_1"], int(row["residue_1"])),
                (row["chain_2"], int(row["residue_2"])),
            )
            for row in csv.DictReader(handle)
        }


def write_auto_reference_disulfide_summary(
    output_root: Path,
    processed: list[tuple[Path, Path]],
    dry_run: bool,
) -> None:
    """Compare detected bonds with the highest residue-mean pLDDT input model."""
    scores = {pdb_path: residue_mean_plddt(pdb_path) for pdb_path, _ in processed}
    usable_scores = {pdb_path: score for pdb_path, score in scores.items() if score is not None}
    summary_path = output_root / "disulfide_reference.json"
    if not usable_scores:
        summary_path.write_text(
            json.dumps(
                {
                    "selection": "unavailable_no_residue_b_factors",
                    "message": "No reference was selected because the input PDBs lack residue B-factor values.",
                },
                indent=2,
            )
            + "\n"
        )
        return

    # `processed` is deterministic input order, so ties select the first path.
    reference_pdb = max(usable_scores, key=lambda path: usable_scores[path])
    reference_output = next(output_dir for pdb_path, output_dir in processed if pdb_path == reference_pdb)
    reference_score = usable_scores[reference_pdb]
    if dry_run:
        summary_path.write_text(
            json.dumps(
                {
                    "selection": "highest_residue_mean_plddt",
                    "reference_input_pdb": str(reference_pdb.resolve()),
                    "reference_residue_mean_plddt": reference_score,
                    "status": "not_compared_dry_run",
                },
                indent=2,
            )
            + "\n"
        )
        return

    reference_report = next(reference_output.glob("*_detected_disulfides.csv"))
    reference_pairs = detected_pair_keys(reference_report)
    model_rows: list[dict[str, object]] = []
    for pdb_path, output_dir in processed:
        manifest_path = next(output_dir.glob("*_prepareforleap_manifest.json"))
        report_path = next(output_dir.glob("*_detected_disulfides.csv"))
        detected_pairs = detected_pair_keys(report_path)
        missing = reference_pairs - detected_pairs
        extra = detected_pairs - reference_pairs
        comparison = {
            "selection": "highest_residue_mean_plddt",
            "reference_input_pdb": str(reference_pdb.resolve()),
            "reference_residue_mean_plddt": reference_score,
            "reference_disulfides": sorted(pair_label(pair) for pair in reference_pairs),
            "matches_reference_disulfides": not missing and not extra,
            "missing_reference_disulfides": sorted(pair_label(pair) for pair in missing),
            "extra_detected_disulfides": sorted(pair_label(pair) for pair in extra),
        }
        manifest = json.loads(manifest_path.read_text())
        manifest["disulfide_reference_comparison"] = comparison
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        model_rows.append(
            {
                "input_pdb": str(pdb_path.resolve()),
                "residue_mean_plddt": scores[pdb_path],
                "matches_reference_disulfides": comparison["matches_reference_disulfides"],
                "missing_reference_disulfides": comparison["missing_reference_disulfides"],
                "extra_detected_disulfides": comparison["extra_detected_disulfides"],
            }
        )
    summary_path.write_text(
        json.dumps(
            {
                "selection": "highest_residue_mean_plddt",
                "reference_input_pdb": str(reference_pdb.resolve()),
                "reference_residue_mean_plddt": reference_score,
                "reference_disulfides": sorted(pair_label(pair) for pair in reference_pairs),
                "models": model_rows,
            },
            indent=2,
        )
        + "\n"
    )


def write_disulfide_report(
    report_path: Path,
    detected: list[dict[str, object]],
    required: list[tuple[tuple[str, int], tuple[str, int]]],
) -> list[tuple[tuple[str, int], tuple[str, int]]]:
    detected_keys = {
        pair_key(
            (str(row["chain_1"]), int(row["residue_1"])),
            (str(row["chain_2"]), int(row["residue_2"])),
        )
        for row in detected
    }
    missing = [pair for pair in required if pair_key(*pair) not in detected_keys]
    with report_path.open("w", newline="") as handle:
        fields = (
            "chain_1",
            "residue_1",
            "resname_1",
            "chain_2",
            "residue_2",
            "resname_2",
            "sg_distance_angstrom",
            "cpptraj_bond",
            "required",
        )
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        required_keys = {pair_key(*pair) for pair in required}
        for row in detected:
            report_row = dict(row)
            report_row["required"] = pair_key(
                (str(row["chain_1"]), int(row["residue_1"])),
                (str(row["chain_2"]), int(row["residue_2"])),
            ) in required_keys
            writer.writerow(report_row)
    return missing


def write_tleap_modifications(
    output_path: Path,
    detected: list[dict[str, object]],
    leap_unit: str,
    missing_required: list[tuple[tuple[str, int], tuple[str, int]]],
    dry_run: bool,
) -> None:
    lines = ["# Paste after loading the prepared PDB into the named LEaP unit."]
    if dry_run:
        lines.append("# Dry run only: CPPTRAJ did not assess disulfides.")
    elif missing_required:
        lines.append("# Not eligible for membrane preparation: required disulfides are missing.")
        for first, second in missing_required:
            lines.append(f"# Missing: {first[0]}:{first[1]}-{second[0]}:{second[1]}")
    elif detected:
        lines.append("# Disulfides detected by CPPTRAJ prepareforleap.")
        lines.extend(
            f"bond {leap_unit}.{row['residue_1']}.SG {leap_unit}.{row['residue_2']}.SG" for row in detected
        )
    else:
        lines.append("# No disulfides detected by CPPTRAJ prepareforleap.")
    output_path.write_text("\n".join(lines) + "\n")


def run_reduce(reduce_exe: str, input_pdb: Path, output_pdb: Path, log_path: Path) -> int:
    """Add hydrogens and optimize Asn/Gln/His flips with Reduce."""
    try:
        result = subprocess.run(
            [reduce_exe, "-BUILD", "-FLIP", str(input_pdb)],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"Could not find Reduce executable {reduce_exe!r}. "
            "Install Reduce or pass --reduce-exe."
        ) from exc
    output_pdb.write_text(result.stdout)
    log_path.write_text(f"Reduce exit code: {result.returncode}\n\n{result.stderr}")
    has_atoms = any(line.startswith(("ATOM  ", "HETATM")) for line in result.stdout.splitlines())
    if not has_atoms:
        raise RuntimeError(f"Reduce completed but did not write {output_pdb}. See {log_path}.")
    return result.returncode


def remove_intermediates(paths: list[Path]) -> list[str]:
    """Remove generated working files and report the names actually removed."""
    removed: list[str] = []
    for path in paths:
        if path.is_file():
            path.unlink()
            removed.append(path.name)
    return removed


def _his_propka_calls(decisions_csv: Path) -> list[str]:
    """Read PROPKA's per-His charge decision, in residue order: 'HIP' or 'neutral'."""
    calls: list[str] = []
    with decisions_csv.open() as handle:
        for row in csv.DictReader(handle):
            if row.get("original_residue", "").strip() == "HIS" or row.get("propka_residue", "").strip() == "HIS":
                calls.append("HIP" if row.get("proposed_amber_residue", "").strip() == "HIP" else "neutral")
    return calls


def _pdb_element(line: str) -> str:
    """Element for an ATOM/HETATM line: cols 77-78 if present, else the leading char of the atom name."""
    elem = line[76:78].strip().upper()
    if elem:
        return elem
    core = line[12:16].strip().lstrip("0123456789")
    return core[0].upper() if core else ""


def reconcile_prepared_pdb(prepared_pdb: Path, decisions_csv: Path | None) -> list[str]:
    """Enforce residue-name atom-set invariants that Reduce/CPPTRAJ can violate.

    Reduce places hydrogens from its own H-bond scoring; CPPTRAJ then renames residues (CYX/HIE/HID/HIP)
    from that placement plus SG-SG geometry. When the two disagree the prepared PDB can carry a hydrogen the
    final residue name forbids, and tleap aborts on the orphan atom ("does not have a type"). Two cases seen:
      * a disulfide cysteine named CYX that still has a thiol HG (Reduce missed a short SG-SG, CPPTRAJ caught
        it) -> strip the HG (a CYX is SG-bonded and never carries HG);
      * a histidine rendered HIP although PROPKA's pKa says neutral -> demote to a neutral tautomer (keep the
        ring N whose nearest external acceptor is closer; drop the other ring N-H), or promote a neutral His
        to HIP where PROPKA proposed HIP (tleap rebuilds the missing ring proton from the template).

    This is a line-level edit of the PDB: it only deletes specific hydrogen lines and renames histidine
    residues, leaving every other record — crucially the TER records and the ACE/NME cap atoms that define
    chain termini — byte-for-byte intact. (An earlier ParmEd round-trip here reconstructed TER/bonds and
    broke cap connectivity, so tleap saw capped residues as bare termini.) PROPKA's decision table
    (decisions_csv) is matched to the structure's histidines in residue order, robust to CPPTRAJ renumbering."""
    lines = prepared_pdb.read_text(errors="replace").splitlines(keepends=True)

    # Group atom lines by residue, in file order; record atom-name -> (line_index, x, y, z).
    residues: list[tuple[str, str, str, dict]] = []
    acceptors: list[tuple[float, float, float, int]] = []  # x,y,z,line_index for every N/O atom
    cur_key = None
    for i, ln in enumerate(lines):
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        chain, resseq, resname, aname = ln[21], ln[22:26], ln[17:20].strip(), ln[12:16].strip()
        x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
        key = (chain, resseq)
        if cur_key != key:
            residues.append((chain, resseq, resname, {})); cur_key = key
        residues[-1][3][aname] = (i, x, y, z)
        if _pdb_element(ln) in ("N", "O"):
            acceptors.append((x, y, z, i))

    changes: list[str] = []
    drop: set[int] = set()        # line indices to delete
    rename: dict[int, str] = {}   # line index -> new 3-char residue name

    def _serial(line_index: int) -> str:
        return lines[line_index][6:11].strip()

    # (1) CYX must not carry a thiol hydrogen.
    for chain, resseq, resname, atoms in residues:
        if resname == "CYX" and "HG" in atoms:
            li = atoms["HG"][0]
            drop.add(li)
            changes.append(f"strip atom {_serial(li)} (HG) from CYX {chain}{resseq.strip()}")

    # (2) Histidine charge must follow PROPKA.
    if decisions_csv and decisions_csv.is_file():
        calls = _his_propka_calls(decisions_csv)
        his = [r for r in residues if r[2] in ("HIS", "HIE", "HID", "HIP")]
        if len(calls) == len(his):
            for (chain, resseq, resname, atoms), call in zip(his, calls):
                if call == "neutral" and resname == "HIP":
                    taut = _pick_neutral_tautomer(atoms, acceptors)
                    drop_atom = "HD1" if taut == "HIE" else "HE2"  # HIE keeps HE2; HID keeps HD1
                    serial = "?"
                    if drop_atom in atoms:
                        li = atoms[drop_atom][0]
                        drop.add(li); serial = _serial(li)
                    for idx, *_ in atoms.values():
                        rename[idx] = taut
                    changes.append(f"demote HIP->{taut} at His {chain}{resseq.strip()} "
                                   f"(PROPKA neutral; strip atom {serial} {drop_atom})")
                elif call == "HIP" and resname in ("HIE", "HID"):
                    for idx, *_ in atoms.values():
                        rename[idx] = "HIP"  # tleap adds the missing ring proton from the template
                    changes.append(f"promote His {chain}{resseq.strip()} -> HIP (PROPKA charged)")
        else:
            changes.append(f"[warn] His count mismatch (structure {len(his)} vs decisions {len(calls)}); "
                           "skipped His reconciliation")

    if drop or rename:
        out: list[str] = []
        for i, ln in enumerate(lines):
            if i in drop:
                continue
            if i in rename:
                ln = ln[:17] + f"{rename[i]:>3}" + ln[20:]
            out.append(ln)
        prepared_pdb.write_text("".join(out))
    return changes


def _pick_neutral_tautomer(atoms: dict, acceptors: list) -> str:
    """For a HIP being neutralized, keep the ring N whose nearest external acceptor is closer.

    `atoms` maps atom-name -> (line_index, x, y, z) for the histidine; `acceptors` is every N/O atom in the
    structure as (x, y, z, line_index). Returns 'HID' (keep Nd1-H) or 'HIE' (keep Ne2-H)."""
    import math
    nd1, ne2 = atoms.get("ND1"), atoms.get("NE2")
    if nd1 is None or ne2 is None:
        return "HIE"
    self_lines = {idx for idx, *_ in atoms.values()}

    def nearest(n):
        _, x, y, z = n
        best = float("inf")
        for ax, ay, az, ai in acceptors:
            if ai in self_lines:
                continue
            d = (x - ax) ** 2 + (y - ay) ** 2 + (z - az) ** 2
            if d < best:
                best = d
        return best

    return "HID" if nearest(nd1) < nearest(ne2) else "HIE"


def process_pdb(
    pdb_path: Path,
    output_root: Path,
    cpptraj_exe: str,
    leap_unit: str,
    keep_water: bool,
    keep_hydrogens: bool,
    detect_histidines: bool,
    detect_disulfides: bool,
    disulfide_sg_cutoff: float,
    altloc: str,
    prepare_sugars: bool,
    use_reduce: bool,
    reduce_exe: str,
    required_disulfides: list[tuple[tuple[str, int], tuple[str, int]]],
    keep_intermediates: bool,
    dry_run: bool,
    protonation_decisions: Path | None = None,
) -> Path:
    output_dir = output_root / pdb_path.stem
    output_dir.mkdir(parents=True, exist_ok=True)
    prepared_pdb = output_dir / f"{pdb_path.stem}_prepareforleap.pdb"
    input_script = output_dir / f"{pdb_path.stem}_prepareforleap.cpptraj.in"
    log_path = output_dir / f"{pdb_path.stem}_prepareforleap.cpptraj.log"
    cpptraj_leap_file = output_dir / f"{pdb_path.stem}_prepareforleap_cpptraj.leap"
    disulfide_report = output_dir / f"{pdb_path.stem}_detected_disulfides.csv"
    tleap_modifications = output_dir / f"{pdb_path.stem}_tleap_modifications.leap"
    reduced_pdb = output_dir / f"{pdb_path.stem}_reduce.pdb"
    reduce_log_path = output_dir / f"{pdb_path.stem}_reduce.log"

    # Reduce supplies hydrogen placement; CPPTRAJ then translates it into
    # HID/HIE/HIP names. Otherwise preserve the previous heavy-atom workflow.
    effective_input = reduced_pdb if use_reduce else pdb_path
    effective_keep_hydrogens = keep_hydrogens or use_reduce
    effective_detect_histidines = detect_histidines or use_reduce
    reduce_exit_code = None
    if use_reduce and not dry_run:
        reduce_exit_code = run_reduce(reduce_exe, pdb_path, reduced_pdb, reduce_log_path)

    script_text = cpptraj_input(
        input_pdb=effective_input.resolve(),
        output_pdb=prepared_pdb.resolve(),
        cpptraj_leap_file=cpptraj_leap_file.resolve(),
        leap_unit=leap_unit,
        keep_water=keep_water,
        keep_hydrogens=effective_keep_hydrogens,
        detect_histidines=effective_detect_histidines,
        detect_disulfides=detect_disulfides,
        disulfide_sg_cutoff=disulfide_sg_cutoff,
        altloc=altloc,
        prepare_sugars=prepare_sugars,
    )
    input_script.write_text(script_text)

    if not dry_run:
        try:
            result = subprocess.run(
                [cpptraj_exe, "-i", str(input_script)],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
        except FileNotFoundError as exc:
            raise RuntimeError(
                f"Could not find CPPTRAJ executable {cpptraj_exe!r}. "
                "Run inside an AmberTools environment or pass --cpptraj-exe."
            ) from exc
        log_path.write_text(result.stdout)
        if result.returncode != 0:
            raise RuntimeError(f"CPPTRAJ failed for {pdb_path}. See {log_path}.")
        if not prepared_pdb.is_file() or prepared_pdb.stat().st_size == 0:
            raise RuntimeError(f"CPPTRAJ completed but did not write {prepared_pdb}. See {log_path}.")

        # Reconcile residue-name atom-set invariants Reduce/CPPTRAJ can violate (orphan CYX HG,
        # PROPKA-contradicting HIP) before the PDB reaches tleap.
        decisions = protonation_decisions
        if decisions is None:
            name = pdb_path.name
            if name.endswith("_amber_states.pdb"):
                candidate = pdb_path.parent / name.replace("_amber_states.pdb", "_protonation_decisions.csv")
                decisions = candidate if candidate.is_file() else None
        reconciliation_changes = reconcile_prepared_pdb(prepared_pdb, decisions)
        if reconciliation_changes:
            (output_dir / f"{pdb_path.stem}_reconciliation.log").write_text(
                "\n".join(reconciliation_changes) + "\n")

    if dry_run:
        detected = []
        missing_required = []
        disulfide_validation = "not_run_dry_run"
    else:
        detected = detected_disulfides(cpptraj_leap_file, prepared_pdb)
        missing_required = write_disulfide_report(disulfide_report, detected, required_disulfides)
        disulfide_validation = "passed" if not missing_required else "failed_missing_required_disulfides"
    write_tleap_modifications(tleap_modifications, detected, leap_unit, missing_required, dry_run)

    intermediate_paths = [reduced_pdb, reduce_log_path, log_path, cpptraj_leap_file]
    # A dry run is useful specifically because it leaves the generated CPPTRAJ
    # recipe for inspection, even when normal intermediates are cleaned up.
    if not dry_run:
        intermediate_paths.append(input_script)
    removed_intermediates = [] if keep_intermediates else remove_intermediates(intermediate_paths)
    manifest = {
        "input_pdb": str(pdb_path.resolve()),
        "prepared_pdb": str(prepared_pdb.resolve()),
        "reduce_used": use_reduce,
        "reduced_pdb": str(reduced_pdb.resolve()) if reduced_pdb.is_file() else None,
        "reduce_log": str(reduce_log_path.resolve()) if reduce_log_path.is_file() else None,
        "reduce_exit_code": reduce_exit_code,
        "reduce_completed_with_warnings": bool(reduce_exit_code) if reduce_exit_code is not None else None,
        "cpptraj_input": str(input_script.resolve()) if input_script.is_file() else None,
        "cpptraj_log": str(log_path.resolve()) if log_path.is_file() else None,
        "cpptraj_leap_file": str(cpptraj_leap_file.resolve()) if cpptraj_leap_file.is_file() else None,
        "tleap_modifications": str(tleap_modifications.resolve()),
        "detected_disulfides": None if dry_run else str(disulfide_report.resolve()),
        "preserved_explicit_histidine_states": not effective_detect_histidines,
        "preserved_explicit_disulfides": not detect_disulfides,
        "disulfide_sg_cutoff_angstrom": disulfide_sg_cutoff if detect_disulfides else None,
        "removed_waters": not keep_water,
        "removed_hydrogens": not effective_keep_hydrogens,
        "required_disulfides": [
            f"{first[0]}:{first[1]}-{second[0]}:{second[1]}" for first, second in required_disulfides
        ],
        "disulfide_validation": disulfide_validation,
        "eligible_for_membrane_preparation": None if dry_run else not missing_required,
        "missing_required_disulfides": [
            f"{first[0]}:{first[1]}-{second[0]}:{second[1]}" for first, second in missing_required
        ],
        "intermediates_kept": keep_intermediates,
        "removed_intermediate_files": removed_intermediates,
    }
    (output_dir / f"{pdb_path.stem}_prepareforleap_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    return output_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input",
        type=Path,
        help="One AMBER-state PDB or a directory containing AMBER-state PDBs.",
    )
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory.")
    parser.add_argument(
        "--pattern",
        default="*_amber_states.pdb",
        help="Directory input filename pattern (default: *_amber_states.pdb).",
    )
    parser.add_argument("--recursive", action="store_true", help="Search input directories recursively.")
    parser.add_argument("--cpptraj-exe", default="cpptraj", help="CPPTRAJ executable (default: cpptraj).")
    parser.add_argument("--leap-unit", default="mol", help="LEaP unit name (default: mol).")
    parser.add_argument("--keep-water", action="store_true", help="Do not strip water molecules.")
    parser.add_argument("--keep-hydrogens", action="store_true", help="Do not strip input hydrogens.")
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help=(
            "Retain Reduce and CPPTRAJ working files and logs. "
            "By default only final handoff files are kept."
        ),
    )
    parser.add_argument(
        "--reduce",
        action="store_true",
        help=(
            "Run Reduce (-BUILD -FLIP) before CPPTRAJ. This adds/optimizes hydrogens "
            "and enables CPPTRAJ histidine-name detection."
        ),
    )
    parser.add_argument(
        "--reduce-exe",
        default="reduce",
        help="Reduce executable used with --reduce (default: reduce).",
    )
    parser.add_argument(
        "--detect-histidines",
        action="store_true",
        help="Allow prepareforleap to rename His based on existing H atoms; off by default.",
    )
    parser.add_argument(
        "--detect-disulfides",
        action="store_true",
        help="Allow prepareforleap to search for disulfides; off by default.",
    )
    parser.add_argument(
        "--disulfide-sg-cutoff",
        type=float,
        default=2.5,
        help="SG-SG cutoff in Angstrom when --detect-disulfides is set (default: 2.5).",
    )
    parser.add_argument(
        "--required-disulfide",
        action="append",
        default=[],
        metavar="CHAIN:RESIDUE-CHAIN:RESIDUE",
        help=(
            "Required CPPTRAJ-detected disulfide; repeatable. Models missing any required pair "
            "are marked ineligible for membrane preparation. Example: A:75-A:160."
        ),
    )
    parser.add_argument(
        "--altloc",
        choices=("A", "highestocc"),
        default="highestocc",
        help="Alternate-location policy (default: highestocc).",
    )
    parser.add_argument(
        "--no-sugars",
        action="store_true",
        help=(
            "Skip prepareforleap's sugar detection/fixing pass. Recommended for apo protein "
            "receptors with no modelled glycans; avoids a cpptraj sugar-prep failure on inputs "
            "whose atom elements it cannot resolve."
        ),
    )
    parser.add_argument("--dry-run", action="store_true", help="Write the CPPTRAJ/LEaP input files but do not run CPPTRAJ.")
    parser.add_argument(
        "--protonation-decisions",
        type=Path,
        default=None,
        help="PROPKA protonation_decisions.csv used to reconcile His charge states (HIP vs neutral) after "
             "CPPTRAJ. If omitted, it is auto-derived from an *_amber_states.pdb input when present.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.disulfide_sg_cutoff <= 0:
        raise ValueError("--disulfide-sg-cutoff must be positive.")
    required_disulfides = [parse_required_disulfide(value) for value in args.required_disulfide]
    if required_disulfides and not args.detect_disulfides:
        raise ValueError("--required-disulfide requires --detect-disulfides; required bonds are never forced.")
    pdbs = input_pdbs(args.input.resolve(), args.recursive, args.pattern)
    processed: list[tuple[Path, Path]] = []
    for pdb_path in pdbs:
        output_dir = process_pdb(
            pdb_path=pdb_path,
            output_root=args.out_dir.resolve(),
            cpptraj_exe=args.cpptraj_exe,
            leap_unit=args.leap_unit,
            keep_water=args.keep_water,
            keep_hydrogens=args.keep_hydrogens,
            detect_histidines=args.detect_histidines,
            detect_disulfides=args.detect_disulfides,
            disulfide_sg_cutoff=args.disulfide_sg_cutoff,
            altloc=args.altloc,
            prepare_sugars=not args.no_sugars,
            use_reduce=args.reduce,
            reduce_exe=args.reduce_exe,
            required_disulfides=required_disulfides,
            keep_intermediates=args.keep_intermediates,
            dry_run=args.dry_run,
            protonation_decisions=args.protonation_decisions.resolve() if args.protonation_decisions else None,
        )
        processed.append((pdb_path, output_dir))
        print(f"{pdb_path.name}: prepareforleap files written to {output_dir}")
    if not required_disulfides:
        write_auto_reference_disulfide_summary(args.out_dir.resolve(), processed, args.dry_run)
        print(f"Automatic disulfide reference written to {args.out_dir.resolve() / 'disulfide_reference.json'}")
    else:
        # Explicit biological requirements supersede any older automatic
        # comparison created in the same output directory.
        (args.out_dir.resolve() / "disulfide_reference.json").unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
