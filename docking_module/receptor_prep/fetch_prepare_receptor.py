#!/usr/bin/env python3
"""Fetch or read a PDB complex, prepare the receptor, and extract a ligand.

This utility is intended to create simple starting inputs for the docking
pipelines in this repository. It can download a structure by PDB ID from RCSB
or use a local complex PDB file.
"""

from __future__ import annotations

import argparse
import math
import shlex
import shutil
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

PH = 7.4
RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_CCD_SDF_URL = "https://files.rcsb.org/ligands/download/{resname}_ideal.sdf"
MISSING_RESIDUE_NEAR_LIGAND_CUTOFF = 6.0
MISSING_RESIDUE_WARNING_THRESHOLD = 5


def _import_receptor_dependencies():
    try:
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer
    except ImportError as exc:
        raise ImportError(
            "Receptor preparation requires pdbfixer and openmm. Install them in "
            "your docking environment, for example: "
            "micromamba install -n docking_pipeline -c conda-forge pdbfixer openmm"
        ) from exc

    return PDBFixer, PDBFile


def _import_minimization_dependencies():
    try:
        from openmm import CustomExternalForce, LangevinIntegrator, unit
        from openmm.app import ForceField, HBonds, NoCutoff, Simulation
    except ImportError as exc:
        raise ImportError(
            "Restrained receptor minimization requires openmm. Install it in "
            "your docking environment, for example: "
            "micromamba install -n docking_pipeline -c conda-forge openmm"
        ) from exc

    return CustomExternalForce, ForceField, HBonds, LangevinIntegrator, NoCutoff, Simulation, unit


def _import_ligand_dependencies():
    try:
        from Bio.PDB import PDBIO, PDBParser, Select
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as exc:
        raise ImportError(
            "Ligand extraction/conversion requires biopython and rdkit. Install "
            "them in your docking environment before running this script."
        ) from exc

    return PDBIO, PDBParser, Select, Chem, AllChem


def fetch_pdb(pdb_id: str, output_pdb: Path) -> Path:
    pdb_id = pdb_id.strip().upper()
    if not pdb_id:
        raise ValueError("PDB ID cannot be empty.")

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    url = RCSB_PDB_URL.format(pdb_id=pdb_id)

    try:
        with urllib.request.urlopen(url) as response:
            contents = response.read()
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"Could not fetch PDB ID {pdb_id} from RCSB: HTTP {exc.code}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Could not fetch PDB ID {pdb_id} from RCSB: {exc.reason}") from exc

    if not contents.strip():
        raise RuntimeError(f"Downloaded PDB file for {pdb_id} was empty.")

    output_pdb.write_bytes(contents)
    return output_pdb


def fetch_ccd_template_sdf(resname: str, output_sdf: Path) -> Path:
    resname = resname.strip().upper()
    if not resname:
        raise ValueError("Ligand residue name cannot be empty.")

    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    url = RCSB_CCD_SDF_URL.format(resname=resname)

    try:
        with urllib.request.urlopen(url) as response:
            contents = response.read()
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"Could not fetch CCD template for ligand {resname}: HTTP {exc.code}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Could not fetch CCD template for ligand {resname}: {exc.reason}") from exc

    if not contents.strip():
        raise RuntimeError(f"Downloaded CCD template for ligand {resname} was empty.")

    output_sdf.write_bytes(contents)
    return output_sdf


def prepare_receptor(
    input_pdb: Path,
    output_pdb: Path,
    receptor_chains: list[str] | None = None,
    remove_receptor_residues: list[str] | None = None,
    keep_water: bool = False,
    add_missing_residues: bool = False,
    add_terminal_missing_residues: bool = False,
    ph: float = PH,
    report_path: Path | None = None,
    ligand_resname: str | None = None,
    ligand_chain: str | None = None,
    ligand_resseq: int | None = None,
    missing_residue_cutoff: float = MISSING_RESIDUE_NEAR_LIGAND_CUTOFF,
    missing_residue_warning_threshold: int = MISSING_RESIDUE_WARNING_THRESHOLD,
    minimize_receptor: bool = False,
    minimize_max_iterations: int = 200,
    minimize_restraint_k: float = 10.0,
    prepared_all_chains_pdb: Path | None = None,
    keep_all_files: bool = True,
    input_complex_is_temporary: bool = False,
    prepared_complex_pdb: Path | None = None,
) -> Path:
    PDBFixer, PDBFile = _import_receptor_dependencies()

    fixer = PDBFixer(filename=str(input_pdb))

    fixer.findMissingResidues()
    topology_chain_ids = _topology_chain_ids(fixer)
    detected_missing_residues = dict(fixer.missingResidues)
    missing_residue_context = summarize_missing_residues_near_ligand(
        input_pdb=input_pdb,
        fixer=fixer,
        missing_residues=detected_missing_residues,
        ligand_resname=ligand_resname,
        ligand_chain=ligand_chain,
        ligand_resseq=ligand_resseq,
        cutoff=missing_residue_cutoff,
    )
    terminal_missing_residues = {}
    if add_missing_residues and not add_terminal_missing_residues:
        terminal_missing_residues = remove_terminal_missing_residues(fixer)
    if not add_missing_residues:
        fixer.missingResidues = {}

    fixer.findNonstandardResidues()
    detected_nonstandard_residues = list(fixer.nonstandardResidues)
    fixer.replaceNonstandardResidues()

    fixer.removeHeterogens(keepWater=keep_water)

    fixer.findMissingAtoms()
    detected_missing_atoms = dict(fixer.missingAtoms)
    detected_missing_terminals = dict(fixer.missingTerminals)
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    minimization_status = {
        "enabled": minimize_receptor,
        "applied": False,
        "max_iterations": minimize_max_iterations,
        "restraint_k": minimize_restraint_k,
        "message": "Not requested.",
    }
    post_removal_cleanup_status = {
        "requested": bool(remove_receptor_residues),
        "applied": False,
        "message": "Not requested.",
    }
    if minimize_receptor and not remove_receptor_residues:
        fixer.positions = minimize_receptor_positions(
            topology=fixer.topology,
            positions=fixer.positions,
            max_iterations=minimize_max_iterations,
            restraint_k=minimize_restraint_k,
        )
        minimization_status["applied"] = True
        minimization_status["message"] = "Restrained minimization completed."
    elif minimize_receptor and remove_receptor_residues:
        minimization_status["message"] = "Deferred until after residue removal cleanup."

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    internal_prepared_all_chains_pdb = output_pdb
    if receptor_chains:
        internal_prepared_all_chains_pdb = prepared_all_chains_pdb or output_pdb.with_name(
            f"{output_pdb.stem}_all_chains{output_pdb.suffix}"
        )

    with internal_prepared_all_chains_pdb.open("w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)

    if receptor_chains:
        filter_receptor_chains(internal_prepared_all_chains_pdb, output_pdb, receptor_chains)

    removed_receptor_residue_ranges = []
    if remove_receptor_residues:
        removed_receptor_residue_ranges = remove_receptor_residue_ranges_from_pdb(
            output_pdb,
            output_pdb,
            remove_receptor_residues,
        )
        cleanup_receptor_with_pdbfixer(
            input_pdb=output_pdb,
            output_pdb=output_pdb,
            keep_water=keep_water,
            ph=ph,
            minimize_receptor=minimize_receptor,
            minimize_max_iterations=minimize_max_iterations,
            minimize_restraint_k=minimize_restraint_k,
        )
        post_removal_cleanup_status["applied"] = True
        post_removal_cleanup_status["message"] = (
            "PDBFixer cleanup was rerun after receptor residue removal without remodeling missing residues."
        )
        if minimize_receptor:
            minimization_status["applied"] = True
            minimization_status["message"] = "Restrained minimization completed after residue removal cleanup."

    residue_numbering_gaps = check_residue_numbering_gaps(output_pdb)

    if report_path is not None:
        write_preparation_report(
            report_path=report_path,
            input_pdb=input_pdb,
            input_complex_is_temporary=input_complex_is_temporary,
            output_pdb=output_pdb,
            prepared_all_chains_pdb=internal_prepared_all_chains_pdb if receptor_chains and keep_all_files else None,
            prepared_complex_pdb=prepared_complex_pdb,
            receptor_chains=receptor_chains,
            removed_receptor_residue_ranges=removed_receptor_residue_ranges,
            post_removal_cleanup_status=post_removal_cleanup_status,
            keep_water=keep_water,
            add_missing_residues=add_missing_residues,
            add_terminal_missing_residues=add_terminal_missing_residues,
            ph=ph,
            missing_residues=detected_missing_residues,
            terminal_missing_residues=terminal_missing_residues,
            topology_chain_ids=topology_chain_ids,
            missing_residue_context=missing_residue_context,
            missing_residue_cutoff=missing_residue_cutoff,
            missing_residue_warning_threshold=missing_residue_warning_threshold,
            nonstandard_residues=detected_nonstandard_residues,
            missing_atoms=detected_missing_atoms,
            missing_terminals=detected_missing_terminals,
            residue_numbering_gaps=residue_numbering_gaps,
            minimization_status=minimization_status,
        )

    return output_pdb


def cleanup_receptor_with_pdbfixer(
    input_pdb: Path,
    output_pdb: Path,
    keep_water: bool,
    ph: float,
    minimize_receptor: bool = False,
    minimize_max_iterations: int = 200,
    minimize_restraint_k: float = 10.0,
) -> Path:
    PDBFixer, PDBFile = _import_receptor_dependencies()

    fixer = PDBFixer(filename=str(input_pdb))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=keep_water)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    if minimize_receptor:
        fixer.positions = minimize_receptor_positions(
            topology=fixer.topology,
            positions=fixer.positions,
            max_iterations=minimize_max_iterations,
            restraint_k=minimize_restraint_k,
        )

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    cleanup_pdb = output_pdb.with_name(f"{output_pdb.stem}_cleanup_tmp{output_pdb.suffix}")
    with cleanup_pdb.open("w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
    cleanup_pdb.replace(output_pdb)
    return output_pdb


def minimize_receptor_positions(
    topology: Any,
    positions: Any,
    max_iterations: int = 200,
    restraint_k: float = 10.0,
) -> Any:
    CustomExternalForce, ForceField, HBonds, LangevinIntegrator, NoCutoff, Simulation, unit = (
        _import_minimization_dependencies()
    )

    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    system = forcefield.createSystem(topology, nonbondedMethod=NoCutoff, constraints=HBonds)

    restraint = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", restraint_k * unit.kilocalories_per_mole / unit.angstroms**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    restrained_atoms = 0
    for atom in topology.atoms():
        if atom.element is not None and atom.element.symbol == "H":
            continue
        pos_nm = positions[atom.index].value_in_unit(unit.nanometer)
        restraint.addParticle(atom.index, [pos_nm.x, pos_nm.y, pos_nm.z])
        restrained_atoms += 1

    if restrained_atoms:
        system.addForce(restraint)

    integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    simulation = Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy(maxIterations=max_iterations)
    state = simulation.context.getState(getPositions=True)
    return state.getPositions()


def _topology_chain_ids(fixer: Any) -> list[str]:
    return [chain.id for chain in fixer.topology.chains()]


def _topology_residues_by_chain(fixer: Any) -> dict[str, list[Any]]:
    residues_by_chain: dict[str, list[Any]] = {}
    for chain in fixer.topology.chains():
        residues_by_chain[chain.id] = list(chain.residues())
    return residues_by_chain


def remove_terminal_missing_residues(fixer: Any) -> dict[tuple[int, int], list[str]]:
    """Remove terminal missing-residue segments from PDBFixer's build list.

    PDBFixer reports missing residues from SEQRES records. Terminal tails are
    frequently unresolved in crystal structures and can be poorly constrained
    for docking prep, so ANTIOCH-VS only builds internal gaps by default.
    """
    chains = list(fixer.topology.chains())
    terminal_segments: dict[tuple[int, int], list[str]] = {}

    for key, residue_names in list(fixer.missingResidues.items()):
        chain_index, residue_index = key
        if chain_index >= len(chains):
            continue
        chain_residue_count = len(list(chains[chain_index].residues()))
        if residue_index == 0 or residue_index == chain_residue_count:
            terminal_segments[key] = list(residue_names)
            del fixer.missingResidues[key]

    return terminal_segments


def _chain_residue_atom_coords(structure: Any) -> dict[tuple[str, int], list[Any]]:
    coords: dict[tuple[str, int], list[Any]] = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag, resseq, _icode = residue.id
                if hetflag != " ":
                    continue
                atom_coords = [atom.coord for atom in residue.get_atoms()]
                if atom_coords:
                    coords[(chain.id, resseq)] = atom_coords
        break
    return coords


def check_residue_numbering_gaps(pdb_file: Path) -> dict[str, list[tuple[int, int]]]:
    _PDBIO, PDBParser, _Select, _Chem, _AllChem = _import_ligand_dependencies()
    structure = PDBParser(QUIET=True).get_structure("prepared_receptor", str(pdb_file))
    gaps: dict[str, list[tuple[int, int]]] = {}

    for model in structure:
        for chain in model:
            residue_numbers = []
            for residue in chain:
                hetflag, resseq, _icode = residue.id
                if hetflag == " ":
                    residue_numbers.append(resseq)

            unique_numbers = sorted(set(residue_numbers))
            chain_gaps = [
                (previous, current)
                for previous, current in zip(unique_numbers, unique_numbers[1:])
                if current - previous > 1
            ]
            if chain_gaps:
                gaps[chain.id] = chain_gaps
        break

    return gaps


def _selected_ligand_atom_coords(
    structure: Any,
    ligand_resname: str,
    ligand_chain: str | None,
    ligand_resseq: int | None,
) -> list[Any]:
    coords = []
    for model in structure:
        for chain in model:
            if ligand_chain is not None and chain.id != ligand_chain:
                continue
            for residue in chain:
                _hetflag, resseq, _icode = residue.id
                if residue.resname.strip() != ligand_resname:
                    continue
                if ligand_resseq is not None and resseq != ligand_resseq:
                    continue
                coords.extend(atom.coord for atom in residue.get_atoms())
        break
    return coords


def _min_distance_to_ligand(residue_coords: list[Any], ligand_coords: list[Any]) -> float | None:
    if not residue_coords or not ligand_coords:
        return None

    min_distance = None
    for residue_coord in residue_coords:
        for ligand_coord in ligand_coords:
            distance = float(((residue_coord - ligand_coord) ** 2).sum() ** 0.5)
            if min_distance is None or distance < min_distance:
                min_distance = distance
    return min_distance


def summarize_missing_residues_near_ligand(
    input_pdb: Path,
    fixer: Any,
    missing_residues: dict[tuple[int, int], list[str]],
    ligand_resname: str | None,
    ligand_chain: str | None,
    ligand_resseq: int | None,
    cutoff: float,
) -> dict[str, Any]:
    context: dict[str, Any] = {
        "ligand_selected": False,
        "nearby_segments": [],
        "nearby_missing_residue_count": 0,
    }
    if not missing_residues or not ligand_resname:
        return context

    _PDBIO, PDBParser, _Select, _Chem, _AllChem = _import_ligand_dependencies()
    structure = PDBParser(QUIET=True).get_structure("complex", str(input_pdb))
    ligand_coords = _selected_ligand_atom_coords(
        structure=structure,
        ligand_resname=ligand_resname,
        ligand_chain=ligand_chain,
        ligand_resseq=ligand_resseq,
    )
    context["ligand_selected"] = bool(ligand_coords)
    if not ligand_coords:
        return context

    chain_ids = _topology_chain_ids(fixer)
    residues_by_chain = _topology_residues_by_chain(fixer)
    residue_coords = _chain_residue_atom_coords(structure)

    for (chain_index, residue_index), residue_names in missing_residues.items():
        if chain_index >= len(chain_ids):
            continue

        chain_id = chain_ids[chain_index]
        present_residues = residues_by_chain.get(chain_id, [])
        flank_residues = []
        if residue_index > 0 and residue_index - 1 < len(present_residues):
            flank_residues.append(present_residues[residue_index - 1])
        if residue_index < len(present_residues):
            flank_residues.append(present_residues[residue_index])

        flank_distances = []
        for residue in flank_residues:
            try:
                resseq = int(residue.id)
            except (TypeError, ValueError):
                continue
            distance = _min_distance_to_ligand(
                residue_coords=residue_coords.get((chain_id, resseq), []),
                ligand_coords=ligand_coords,
            )
            if distance is not None:
                flank_distances.append((resseq, distance))

        near_distances = [(resseq, distance) for resseq, distance in flank_distances if distance <= cutoff]
        if not near_distances:
            continue

        context["nearby_segments"].append(
            {
                "chain_id": chain_id,
                "insert_before_index": residue_index,
                "residue_names": list(residue_names),
                "flank_distances": flank_distances,
            }
        )
        context["nearby_missing_residue_count"] += len(residue_names)

    return context


def _format_missing_residue_segment(
    chain_id: str,
    residue_index: int,
    residue_names: list[str],
) -> str:
    residues = ", ".join(residue_names)
    return f"chain {chain_id}, insert_before_residue_index {residue_index}: {len(residue_names)} residues [{residues}]"


def _format_fixer_residue(residue: Any) -> str:
    chain = getattr(residue, "chain", None)
    chain_id = getattr(chain, "id", "?")
    residue_id = getattr(residue, "id", "?")
    residue_name = getattr(residue, "name", "?")
    return f"{chain_id}:{residue_id} {residue_name}"


def _format_atom_name(atom: Any) -> str:
    return str(getattr(atom, "name", atom))


def write_preparation_report(
    report_path: Path,
    input_pdb: Path,
    input_complex_is_temporary: bool,
    output_pdb: Path,
    prepared_all_chains_pdb: Path | None,
    prepared_complex_pdb: Path | None,
    receptor_chains: list[str] | None,
    removed_receptor_residue_ranges: list[dict[str, Any]],
    post_removal_cleanup_status: dict[str, Any],
    keep_water: bool,
    add_missing_residues: bool,
    add_terminal_missing_residues: bool,
    ph: float,
    missing_residues: dict[tuple[int, int], list[str]],
    terminal_missing_residues: dict[tuple[int, int], list[str]],
    topology_chain_ids: list[str],
    missing_residue_context: dict[str, Any],
    missing_residue_cutoff: float,
    missing_residue_warning_threshold: int,
    nonstandard_residues: list[Any],
    missing_atoms: dict[Any, list[Any]],
    missing_terminals: dict[Any, list[Any]],
    residue_numbering_gaps: dict[str, list[tuple[int, int]]],
    minimization_status: dict[str, Any],
) -> Path:
    report_path.parent.mkdir(parents=True, exist_ok=True)

    missing_count = sum(len(residues) for residues in missing_residues.values())
    nearby_count = int(missing_residue_context.get("nearby_missing_residue_count", 0))
    all_chain_output = prepared_all_chains_pdb if prepared_all_chains_pdb else "not kept (temporary)"
    if not receptor_chains:
        all_chain_output = "same as prepared receptor"

    lines = [
        "Receptor Preparation Report",
        "===========================",
        "",
        f"Input complex: {'not kept (downloaded temporarily)' if input_complex_is_temporary else input_pdb}",
        f"Prepared receptor: {output_pdb}",
        f"Prepared all-chain receptor: {all_chain_output}",
        f"Prepared receptor-ligand complex: {prepared_complex_pdb if prepared_complex_pdb else 'not requested'}",
        f"Selected receptor chains: {', '.join(receptor_chains) if receptor_chains else 'all protein chains'}",
        f"pH for added hydrogens: {ph}",
        f"Kept waters: {keep_water}",
        f"Modeled missing residues: {add_missing_residues}",
        f"Modeled terminal missing residues: {add_missing_residues and add_terminal_missing_residues}",
        "",
        "Removed Receptor Residue Ranges",
        "-------------------------------",
    ]

    if removed_receptor_residue_ranges:
        for residue_range in removed_receptor_residue_ranges:
            matched_chains = (
                ", ".join(residue_range["matched_chains"])
                if residue_range["matched_chains"]
                else "none"
            )
            lines.append(
                f"{residue_range['spec']}: removed {residue_range['removed_count']} "
                f"protein residues from chain(s) {matched_chains}"
            )
    else:
        lines.append("No receptor residue ranges were requested for removal.")

    lines.extend(
        [
            "",
            "Post-Removal Receptor Cleanup",
            "-----------------------------",
            f"Requested: {post_removal_cleanup_status['requested']}",
            f"Applied: {post_removal_cleanup_status['applied']}",
            f"Status: {post_removal_cleanup_status['message']}",
        ]
    )

    lines.extend(
        [
            "",
            "Missing Residues",
            "----------------",
            f"Detected missing residues: {missing_count}",
        ]
    )

    if missing_residues:
        for (chain_index, residue_index), residue_names in sorted(missing_residues.items()):
            chain_label = (
                topology_chain_ids[chain_index]
                if chain_index < len(topology_chain_ids)
                else f"chain_index_{chain_index}"
            )
            lines.append(_format_missing_residue_segment(chain_label, residue_index, list(residue_names)))
    else:
        lines.append("No missing residues were reported by PDBFixer.")

    if missing_count and not add_missing_residues:
        lines.append("")
        lines.append("Note: missing residues were detected but not modeled because --add-missing-residues was not used.")
    elif terminal_missing_residues:
        lines.append("")
        lines.append(
            "Note: terminal missing-residue segments were not modeled. ANTIOCH-VS only models internal missing "
            "residue gaps by default because unresolved terminal tails can be poorly constrained for docking prep. "
            "Use --add-terminal-missing-residues to include them."
        )
        for (chain_index, residue_index), residue_names in sorted(terminal_missing_residues.items()):
            chain_label = (
                topology_chain_ids[chain_index]
                if chain_index < len(topology_chain_ids)
                else f"chain_index_{chain_index}"
            )
            lines.append(_format_missing_residue_segment(chain_label, residue_index, list(residue_names)))

    lines.extend(
        [
            "",
            "Missing Residues Near Selected Ligand",
            "-------------------------------------",
            f"Distance cutoff: {missing_residue_cutoff:.1f} A from flanking present residues to ligand atoms",
        ]
    )

    if missing_count == 0:
        lines.append("No missing-residue segments were detected; ligand-neighborhood analysis was not needed.")
    elif not missing_residue_context.get("ligand_selected"):
        lines.append("Ligand-neighborhood analysis was skipped because no matching ligand atoms were selected.")
    elif nearby_count == 0:
        lines.append("No missing-residue segments were detected near the selected ligand.")
    else:
        lines.append(f"Missing residues near ligand: {nearby_count}")
        for segment in missing_residue_context.get("nearby_segments", []):
            lines.append(
                _format_missing_residue_segment(
                    segment["chain_id"],
                    segment["insert_before_index"],
                    segment["residue_names"],
                )
            )
            for resseq, distance in segment["flank_distances"]:
                lines.append(f"  flank residue {segment['chain_id']}:{resseq} distance_to_ligand = {distance:.2f} A")

        if nearby_count >= missing_residue_warning_threshold:
            lines.append("")
            lines.append(
                "WARNING: many missing residues are near the selected ligand. PDBFixer can build simple missing "
                "segments, but this may not be reliable for modeling binding-site loops. Consider a dedicated "
                "loop/modeling workflow and manual inspection before docking."
            )

    lines.extend(
        [
            "",
            "Nonstandard Residues",
            "--------------------",
            f"Detected nonstandard residues: {len(nonstandard_residues)}",
        ]
    )
    for residue, replacement in nonstandard_residues:
        lines.append(f"{_format_fixer_residue(residue)} -> {replacement}")

    missing_atom_count = sum(len(atoms) for atoms in missing_atoms.values())
    missing_terminal_count = sum(len(atoms) for atoms in missing_terminals.values())
    lines.extend(
        [
            "",
            "Missing Atoms",
            "-------------",
            f"Missing heavy atoms added: {missing_atom_count}",
            f"Missing terminal atoms added: {missing_terminal_count}",
        ]
    )

    for residue, atoms in missing_atoms.items():
        atom_names = ", ".join(_format_atom_name(atom) for atom in atoms)
        lines.append(f"{_format_fixer_residue(residue)}: {atom_names}")
    for residue, atoms in missing_terminals.items():
        atom_names = ", ".join(_format_atom_name(atom) for atom in atoms)
        lines.append(f"{_format_fixer_residue(residue)} terminal atoms: {atom_names}")

    lines.extend(
        [
            "",
            "Final Residue Numbering Continuity",
            "----------------------------------",
        ]
    )
    if residue_numbering_gaps:
        lines.append("Residue numbering gaps remain in the prepared receptor:")
        for chain_id, gaps in residue_numbering_gaps.items():
            for previous, current in gaps:
                lines.append(f"chain {chain_id}: numbering jumps from {previous} to {current}")
        lines.append(
            "Note: numbering gaps can be harmless if they reflect terminal truncation or intentional residue numbering, "
            "but internal pocket gaps should be inspected before docking."
        )
    else:
        lines.append("No residue numbering gaps were detected in protein residues of the prepared receptor.")

    lines.extend(
        [
            "",
            "Restrained Minimization",
            "-----------------------",
            f"Requested: {minimization_status['enabled']}",
            f"Applied: {minimization_status['applied']}",
            f"Max iterations: {minimization_status['max_iterations']}",
            f"Heavy-atom restraint strength: {minimization_status['restraint_k']} kcal/mol/A^2",
            f"Status: {minimization_status['message']}",
        ]
    )

    report_path.write_text("\n".join(lines) + "\n")
    return report_path


def _make_ligand_select_class():
    _, _, Select, _, _ = _import_ligand_dependencies()

    class LigandSelect(Select):
        def __init__(self, chain_id: str | None, resname: str, resseq: int | None = None):
            self.chain_id = chain_id
            self.resname = resname
            self.resseq = resseq

        def accept_residue(self, residue):
            _hetflag, resseq, _icode = residue.id
            if residue.resname.strip() != self.resname:
                return 0
            if self.chain_id is not None and residue.parent.id != self.chain_id:
                return 0
            if self.resseq is not None and resseq != self.resseq:
                return 0
            return 1

    return LigandSelect


def _make_receptor_chain_select_class():
    _, _, Select, _, _ = _import_ligand_dependencies()

    class ReceptorChainSelect(Select):
        def __init__(self, chain_ids: set[str]):
            self.chain_ids = chain_ids

        def accept_chain(self, chain):
            return 1 if chain.id in self.chain_ids else 0

        def accept_residue(self, residue):
            hetflag, _resseq, _icode = residue.id
            return 1 if hetflag == " " or residue.resname.strip() == "HOH" else 0

    return ReceptorChainSelect


def filter_receptor_chains(input_pdb: Path, output_pdb: Path, receptor_chains: list[str]) -> Path:
    PDBIO, PDBParser, _Select, _Chem, _AllChem = _import_ligand_dependencies()
    ReceptorChainSelect = _make_receptor_chain_select_class()

    selected_chains = {chain.strip() for chain in receptor_chains if chain.strip()}
    if not selected_chains:
        raise ValueError("--receptor-chains was provided but no valid chain IDs were found.")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", str(input_pdb))
    available_chains = {chain.id for model in structure for chain in model}
    missing_chains = sorted(selected_chains - available_chains)
    if missing_chains:
        raise ValueError(
            "Requested receptor chain(s) not found in input PDB: "
            + ", ".join(missing_chains)
            + f". Available chains: {', '.join(sorted(available_chains))}"
        )

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb), ReceptorChainSelect(selected_chains))
    return output_pdb


def parse_receptor_residue_range_specs(specs: list[str]) -> list[dict[str, Any]]:
    ranges = []
    for spec in specs:
        for item in spec.split(","):
            item = item.strip()
            if not item:
                continue

            chain_id = None
            residue_range = item
            if ":" in item:
                chain_id, residue_range = item.split(":", 1)
                chain_id = chain_id.strip()
                residue_range = residue_range.strip()
                if not chain_id:
                    raise ValueError(f"Invalid residue-removal range {item!r}: chain ID is empty.")

            if "-" in residue_range:
                start_text, end_text = residue_range.split("-", 1)
            else:
                start_text = end_text = residue_range

            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid residue-removal range {item!r}. Use START-END or CHAIN:START-END."
                ) from exc

            if start > end:
                raise ValueError(f"Invalid residue-removal range {item!r}: start is greater than end.")

            ranges.append({"chain_id": chain_id, "start": start, "end": end, "spec": item})

    if not ranges:
        raise ValueError("--remove-receptor-residues was provided but no valid residue ranges were found.")

    return ranges


def _make_remove_residue_range_select_class(ranges: list[dict[str, Any]]):
    _, _, Select, _, _ = _import_ligand_dependencies()

    class RemoveResidueRangeSelect(Select):
        def accept_residue(self, residue):
            hetflag, resseq, _icode = residue.id
            if hetflag != " ":
                return 1

            chain_id = residue.parent.id
            for residue_range in ranges:
                requested_chain = residue_range.get("chain_id")
                if requested_chain is not None and chain_id != requested_chain:
                    continue
                if residue_range["start"] <= resseq <= residue_range["end"]:
                    return 0

            return 1

    return RemoveResidueRangeSelect


def remove_receptor_residue_ranges_from_pdb(input_pdb: Path, output_pdb: Path, specs: list[str]) -> list[dict[str, Any]]:
    PDBIO, PDBParser, _Select, _Chem, _AllChem = _import_ligand_dependencies()
    ranges = parse_receptor_residue_range_specs(specs)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prepared_receptor", str(input_pdb))
    available_chains = {chain.id for model in structure for chain in model}
    requested_chains = {residue_range["chain_id"] for residue_range in ranges if residue_range.get("chain_id")}
    missing_chains = sorted(requested_chains - available_chains)
    if missing_chains:
        raise ValueError(
            "Requested residue-removal chain(s) not found in prepared receptor: "
            + ", ".join(missing_chains)
            + f". Available chains: {', '.join(sorted(available_chains))}"
        )

    removed_ranges = []
    for residue_range in ranges:
        removed_count = 0
        matched_chains = set()
        for model in structure:
            for chain in model:
                requested_chain = residue_range.get("chain_id")
                if requested_chain is not None and chain.id != requested_chain:
                    continue
                for residue in chain:
                    hetflag, resseq, _icode = residue.id
                    if hetflag == " " and residue_range["start"] <= resseq <= residue_range["end"]:
                        removed_count += 1
                        matched_chains.add(chain.id)
            break

        removed_ranges.append(
            {
                "spec": residue_range["spec"],
                "chain_id": residue_range.get("chain_id"),
                "start": residue_range["start"],
                "end": residue_range["end"],
                "matched_chains": sorted(matched_chains),
                "removed_count": removed_count,
            }
        )

    RemoveResidueRangeSelect = _make_remove_residue_range_select_class(ranges)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb), RemoveResidueRangeSelect())

    return removed_ranges


def extract_ligand(
    input_pdb: Path,
    output_pdb: Path,
    resname: str,
    chain_id: str | None = None,
    resseq: int | None = None,
) -> Path:
    PDBIO, PDBParser, _, _, _ = _import_ligand_dependencies()
    LigandSelect = _make_ligand_select_class()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", str(input_pdb))
    output_pdb.parent.mkdir(parents=True, exist_ok=True)

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb), LigandSelect(chain_id=chain_id, resname=resname, resseq=resseq))

    atom_count = sum(
        1 for line in output_pdb.read_text().splitlines() if line.startswith(("ATOM", "HETATM"))
    )
    if atom_count == 0:
        raise RuntimeError(
            "Ligand extraction produced an empty file. Check --ligand-resname, "
            "--ligand-chain, and --ligand-resseq."
        )

    return output_pdb


def _pdb_coord_lines(input_pdb: Path, record_prefixes: tuple[str, ...]) -> list[str]:
    return [
        line
        for line in input_pdb.read_text().splitlines()
        if line.startswith(record_prefixes)
    ]


def write_prepared_complex_pdb(receptor_pdb: Path, ligand_pdb: Path, output_pdb: Path) -> Path:
    receptor_lines = _pdb_coord_lines(receptor_pdb, ("ATOM", "HETATM", "TER"))
    ligand_lines = _pdb_coord_lines(ligand_pdb, ("HETATM",))

    if not receptor_lines:
        raise RuntimeError(f"No receptor coordinate records found in {receptor_pdb}")
    if not ligand_lines:
        raise RuntimeError(f"No ligand HETATM records found in {ligand_pdb}")

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    output_pdb.write_text("\n".join(receptor_lines + ligand_lines + ["END"]) + "\n")
    return output_pdb


def copy_prepared_ligand_mol2(input_mol2: Path, output_mol2: Path) -> Path:
    if not input_mol2.exists():
        raise FileNotFoundError(f"Prepared ligand MOL2 does not exist: {input_mol2}")

    output_mol2.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(input_mol2, output_mol2)
    return output_mol2


def _parse_mol2_charge_sum(input_mol2: Path) -> float | None:
    in_atom_section = False
    charge_sum = 0.0
    charge_count = 0

    for line in input_mol2.read_text().splitlines():
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom_section = True
            continue
        if line.startswith("@<TRIPOS>") and in_atom_section:
            break
        if not in_atom_section or not line.strip():
            continue

        fields = line.split()
        if len(fields) < 9:
            continue
        try:
            charge_sum += float(fields[-1])
            charge_count += 1
        except ValueError:
            continue

    return charge_sum if charge_count else None


def ligand_center_of_mass(mol: Any) -> tuple[float, float, float]:
    if mol.GetNumConformers() == 0:
        raise RuntimeError("Cannot calculate ligand center of mass because the molecule has no coordinates.")

    conformer = mol.GetConformer()
    total_mass = 0.0
    weighted_x = 0.0
    weighted_y = 0.0
    weighted_z = 0.0

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        mass = float(atom.GetMass())
        position = conformer.GetAtomPosition(atom.GetIdx())
        total_mass += mass
        weighted_x += mass * position.x
        weighted_y += mass * position.y
        weighted_z += mass * position.z

    if total_mass == 0:
        raise RuntimeError("Cannot calculate ligand center of mass because no heavy atoms were found.")

    return (weighted_x / total_mass, weighted_y / total_mass, weighted_z / total_mass)


def run_mk_prepare_ligand(input_ligand: Path, output_pdbqt: Path, command: str) -> Path:
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    cmd = [command, "-i", str(input_ligand), "-o", str(output_pdbqt)]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"Could not find {command!r}. Install Meeko or pass --mk-prepare-ligand-command."
        ) from exc
    except subprocess.CalledProcessError as exc:
        output = "\n".join(part for part in [exc.stdout, exc.stderr] if part)
        raise RuntimeError(f"mk_prepare_ligand.py failed while writing {output_pdbqt}:\n{output}") from exc

    if result.stdout.strip():
        print(result.stdout.strip())
    if result.stderr.strip():
        print(result.stderr.strip(), file=sys.stderr)

    if not output_pdbqt.exists():
        raise RuntimeError(f"mk_prepare_ligand.py finished but did not create {output_pdbqt}")

    return output_pdbqt


def _pdb_atom_is_hydrogen(line: str) -> bool:
    """Return True when a PDB ATOM/HETATM record represents hydrogen."""
    atom_name = line[12:16].strip()
    element = line[76:78].strip() if len(line) >= 78 else ""
    return element.upper() == "H" or atom_name.upper().startswith("H")


def write_receptor_pdb_for_meeko(input_receptor: Path, output_receptor: Path) -> Path:
    """Write a receptor PDB copy suitable for Meeko receptor conversion.

    Meeko builds an RDKit molecule from the receptor coordinates. Passing a
    PDBFixer-prepared receptor verbatim can leave records that make RDKit infer
    impossible valences, especially pre-existing hydrogens, alternate locations,
    or CONECT records. The staged file keeps the receptor protein coordinates
    while letting Meeko assign receptor hydrogens/charges itself.
    """
    output_receptor.parent.mkdir(parents=True, exist_ok=True)
    with input_receptor.open() as source, output_receptor.open("w") as dest:
        for line in source:
            if line.startswith("CONECT"):
                continue
            if line.startswith(("ATOM", "HETATM")):
                if _pdb_atom_is_hydrogen(line):
                    continue
                altloc = line[16].strip() if len(line) > 16 else ""
                if altloc and altloc != "A":
                    continue
                if altloc == "A":
                    line = f"{line[:16]} {line[17:]}"
                if line.startswith("HETATM"):
                    continue
            dest.write(line)
    return output_receptor


def write_crystal_receptor_pdb_for_meeko(
    input_complex: Path,
    output_receptor: Path,
    receptor_chains: list[str] | None = None,
) -> Path:
    """Stage a clean heavy-atom protein receptor directly from the input complex.

    This is useful for Vina/Meeko when PDBFixer-added atoms or terminal handling
    create geometry that RDKit cannot sanitize. HETATM records are intentionally
    omitted so waters, ions, cofactors, and the reference ligand are not carried
    into the receptor PDBQT conversion.
    """
    selected_chains = set(receptor_chains or [])
    wrote_atom = False

    output_receptor.parent.mkdir(parents=True, exist_ok=True)
    with input_complex.open() as source, output_receptor.open("w") as dest:
        for line in source:
            if not line.startswith("ATOM"):
                continue
            chain_id = line[21].strip() if len(line) > 21 else ""
            if selected_chains and chain_id not in selected_chains:
                continue
            if _pdb_atom_is_hydrogen(line):
                continue
            altloc = line[16].strip() if len(line) > 16 else ""
            if altloc and altloc != "A":
                continue
            if altloc == "A":
                line = f"{line[:16]} {line[17:]}"
            dest.write(line)
            wrote_atom = True

        if wrote_atom:
            dest.write("TER\n")
        dest.write("END\n")

    if not wrote_atom:
        raise RuntimeError(
            "Could not stage a crystal receptor for Meeko because no matching protein ATOM records were found."
        )

    return output_receptor


def run_mk_prepare_receptor_from_staged_input(
    input_receptor: Path,
    output_pdbqt: Path,
    command: str,
    extra_args: list[str] | None = None,
) -> Path:
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    cmd = [command, "-i", str(input_receptor), "--write_pdbqt", str(output_pdbqt)]
    cmd.extend(extra_args or [])

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"Could not find {command!r}. Install Meeko or pass --mk-prepare-receptor-command."
        ) from exc
    except subprocess.CalledProcessError as exc:
        output = "\n".join(part for part in [exc.stdout, exc.stderr] if part)
        raise RuntimeError(f"mk_prepare_receptor.py failed while writing {output_pdbqt}:\n{output}") from exc

    if result.stdout.strip():
        print(result.stdout.strip())
    if result.stderr.strip():
        print(result.stderr.strip(), file=sys.stderr)

    if not output_pdbqt.exists():
        raise RuntimeError(f"mk_prepare_receptor.py finished but did not create {output_pdbqt}")

    return output_pdbqt


def run_mk_prepare_receptor(
    input_receptor: Path,
    output_pdbqt: Path,
    command: str,
    extra_args: list[str] | None = None,
) -> Path:
    with tempfile.TemporaryDirectory(prefix="dockflow_meeko_receptor_") as tmpdir_name:
        meeko_input = Path(tmpdir_name) / input_receptor.name
        write_receptor_pdb_for_meeko(input_receptor, meeko_input)
        return run_mk_prepare_receptor_from_staged_input(meeko_input, output_pdbqt, command, extra_args)


def run_mk_prepare_receptor_with_fallback(
    prepared_receptor: Path,
    input_complex: Path,
    output_pdbqt: Path,
    command: str,
    receptor_chains: list[str] | None = None,
    vina_receptor_source: str = "auto",
    extra_args: list[str] | None = None,
) -> Path:
    if vina_receptor_source == "prepared":
        return run_mk_prepare_receptor(prepared_receptor, output_pdbqt, command, extra_args)

    with tempfile.TemporaryDirectory(prefix="dockflow_meeko_receptor_") as tmpdir_name:
        tmpdir = Path(tmpdir_name)

        if vina_receptor_source == "crystal":
            crystal_input = tmpdir / f"{input_complex.stem}_crystal_receptor_for_meeko.pdb"
            write_crystal_receptor_pdb_for_meeko(input_complex, crystal_input, receptor_chains)
            return run_mk_prepare_receptor_from_staged_input(crystal_input, output_pdbqt, command, extra_args)

        prepared_input = tmpdir / prepared_receptor.name
        write_receptor_pdb_for_meeko(prepared_receptor, prepared_input)
        try:
            return run_mk_prepare_receptor_from_staged_input(prepared_input, output_pdbqt, command, extra_args)
        except RuntimeError as prepared_error:
            crystal_input = tmpdir / f"{input_complex.stem}_crystal_receptor_for_meeko.pdb"
            write_crystal_receptor_pdb_for_meeko(input_complex, crystal_input, receptor_chains)
            print(
                "[!] Meeko failed on the PDBFixer-prepared receptor; retrying with a clean "
                "heavy-atom receptor staged from the input crystal complex.",
                file=sys.stderr,
            )
            try:
                return run_mk_prepare_receptor_from_staged_input(crystal_input, output_pdbqt, command, extra_args)
            except RuntimeError as crystal_error:
                raise RuntimeError(
                    "mk_prepare_receptor.py failed for both Vina receptor sources.\n\n"
                    f"Prepared receptor error:\n{prepared_error}\n\n"
                    f"Crystal receptor fallback error:\n{crystal_error}"
                ) from crystal_error


def ligand_pdb_to_sdf(
    input_pdb: Path,
    output_sdf: Path,
    resname: str | None = None,
    template_sdf: Path | None = None,
    use_template_bond_orders: bool = True,
    optimize: bool = True,
    output_reference_mol2: Path | None = None,
    output_mol2: Path | None = None,
    mol2_charge_model: str = "gasteiger",
    prepared_ligand_mol2: Path | None = None,
    report_path: Path | None = None,
    keep_all_files: bool = True,
) -> dict[str, Any]:
    _, _, _, Chem, AllChem = _import_ligand_dependencies()

    mol = Chem.MolFromPDBFile(str(input_pdb), removeHs=False, sanitize=False, proximityBonding=True)
    if mol is None:
        raise ValueError("RDKit could not parse ligand from PDB. Bonding may need manual repair.")

    template_sdf_used = template_sdf
    template_source = "none"
    if use_template_bond_orders:
        if template_sdf is None:
            if resname is None:
                raise ValueError("resname is required when assigning CCD template bond orders.")
            template_sdf = output_sdf.parent / f"{resname.strip().upper()}_ideal.sdf"
            fetch_ccd_template_sdf(resname, template_sdf)
            template_source = "RCSB CCD"
        else:
            template_source = "user-provided"
        template_sdf_used = template_sdf

        template_supplier = Chem.SDMolSupplier(str(template_sdf), removeHs=False, sanitize=True)
        template = template_supplier[0] if template_supplier and len(template_supplier) else None
        if template is None:
            raise ValueError(f"RDKit could not parse CCD template SDF: {template_sdf}")

        try:
            mol = Chem.RemoveHs(mol, sanitize=False)
            template = Chem.RemoveHs(template)
            mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
            Chem.SanitizeMol(mol)
        except Exception as exc:
            raise ValueError(
                "Could not assign ligand bond orders from the CCD template. "
                "For custom ligands, provide --ccd-template-sdf or rerun with --no-template-bond-orders."
            ) from exc
    else:
        Chem.SanitizeMol(mol)

    formal_charge = Chem.GetFormalCharge(mol)
    heavy_atom_count = mol.GetNumHeavyAtoms()
    reference_center_of_mass = ligand_center_of_mass(mol)

    if output_reference_mol2 is not None:
        reference_mol = Chem.RemoveHs(mol, sanitize=False)
        write_ligand_mol2(reference_mol, output_reference_mol2, resname=resname, charge_model="formal")

    mol = Chem.AddHs(mol, addCoords=True)
    if optimize:
        force_field = AllChem.UFFGetMoleculeForceField(mol)
        if force_field is None:
            raise RuntimeError("RDKit could not build a UFF force field for ligand optimization.")

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() > 1:
                force_field.AddFixedPoint(atom.GetIdx())

        force_field.Initialize()
        force_field.Minimize(maxIts=500)

    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_sdf))
    try:
        writer.write(mol)
    finally:
        writer.close()

    if output_mol2 is not None:
        if prepared_ligand_mol2 is not None:
            copy_prepared_ligand_mol2(prepared_ligand_mol2, output_mol2)
        else:
            write_ligand_mol2(mol, output_mol2, resname=resname, charge_model=mol2_charge_model)

    ligand_report = {
        "input_pdb": input_pdb,
        "output_sdf": output_sdf,
        "output_reference_mol2": output_reference_mol2,
        "output_mol2": output_mol2,
        "resname": resname,
        "template_sdf": template_sdf_used,
        "template_source": template_source,
        "used_template_bond_orders": use_template_bond_orders,
        "formal_charge": formal_charge,
        "reference_center_of_mass": reference_center_of_mass,
        "heavy_atom_count": heavy_atom_count,
        "atom_count_with_hydrogens": mol.GetNumAtoms(),
        "optimized": optimize,
        "mol2_charge_model": mol2_charge_model if output_mol2 is not None else None,
        "prepared_ligand_mol2": prepared_ligand_mol2,
        "keep_all_files": keep_all_files,
    }
    if output_mol2 is not None:
        if prepared_ligand_mol2 is not None:
            ligand_report["mol2_charge_sum"] = _parse_mol2_charge_sum(output_mol2)
            ligand_report["mol2_charge_model"] = "user-provided"
        else:
            ligand_report["mol2_charge_sum"] = sum(_mol2_atom_charge(atom) for atom in mol.GetAtoms())

    if report_path is not None:
        append_ligand_report(report_path, ligand_report)

    return ligand_report


def _mol2_atom_name(atom: Any, used_names: set[str]) -> str:
    residue_info = atom.GetPDBResidueInfo()
    if residue_info is not None:
        name = residue_info.GetName().strip()
    else:
        name = ""

    if not name:
        name = f"{atom.GetSymbol()}{atom.GetIdx() + 1}"

    base_name = "".join(char if char.isalnum() or char in "_'" else "_" for char in name)[:8] or "A"
    atom_name = base_name
    counter = 1
    while atom_name in used_names:
        suffix = str(counter)
        atom_name = f"{base_name[: 8 - len(suffix)]}{suffix}"
        counter += 1

    used_names.add(atom_name)
    return atom_name


def _mol2_atom_type(atom: Any) -> str:
    symbol = atom.GetSymbol()
    hybridization = str(atom.GetHybridization())

    if atom.GetIsAromatic():
        return f"{symbol}.ar"
    if symbol == "C":
        if hybridization.endswith("SP3"):
            return "C.3"
        if hybridization.endswith("SP2"):
            return "C.2"
        if hybridization.endswith("SP"):
            return "C.1"
    if symbol == "N":
        if atom.GetFormalCharge() > 0:
            return "N.4"
        if hybridization.endswith("SP3"):
            return "N.3"
        if hybridization.endswith("SP2"):
            return "N.2"
        if hybridization.endswith("SP"):
            return "N.1"
    if symbol == "O":
        if atom.GetFormalCharge() < 0:
            return "O.co2"
        if hybridization.endswith("SP3"):
            return "O.3"
        return "O.2"
    if symbol == "S":
        if hybridization.endswith("SP3"):
            return "S.3"
        return "S.2"
    if symbol == "P":
        return "P.3"

    return symbol


def _mol2_bond_type(bond: Any) -> str:
    if bond.GetIsAromatic():
        return "ar"

    bond_type = str(bond.GetBondType())
    if bond_type == "DOUBLE":
        return "2"
    if bond_type == "TRIPLE":
        return "3"
    if bond_type == "AROMATIC":
        return "ar"
    return "1"


def _mol2_atom_charge(atom: Any) -> float:
    if atom.HasProp("_TriposPartialCharge"):
        return float(atom.GetProp("_TriposPartialCharge"))
    if atom.HasProp("_GasteigerCharge"):
        charge = float(atom.GetProp("_GasteigerCharge"))
        if not math.isnan(charge) and not math.isinf(charge):
            return charge
    if atom.HasProp("molFileAlias"):
        return 0.0
    return float(atom.GetFormalCharge())


def write_ligand_mol2(
    mol: Any,
    output_mol2: Path,
    resname: str | None = None,
    charge_model: str = "gasteiger",
) -> Path:
    if mol.GetNumConformers() == 0:
        raise RuntimeError("Cannot write ligand MOL2 because the molecule has no coordinates.")

    if charge_model == "gasteiger":
        _, _, _, _Chem, AllChem = _import_ligand_dependencies()
        AllChem.ComputeGasteigerCharges(mol)
        for atom in mol.GetAtoms():
            charge = float(atom.GetProp("_GasteigerCharge")) if atom.HasProp("_GasteigerCharge") else math.nan
            if math.isnan(charge) or math.isinf(charge):
                charge = float(atom.GetFormalCharge())
            atom.SetProp("_TriposPartialCharge", f"{charge:.6f}")
    elif charge_model == "formal":
        for atom in mol.GetAtoms():
            atom.SetProp("_TriposPartialCharge", f"{float(atom.GetFormalCharge()):.6f}")
    else:
        raise ValueError(f"Unsupported MOL2 charge model: {charge_model}")

    ligand_name = resname.strip().upper() if resname else "LIG"
    conformer = mol.GetConformer()
    used_names: set[str] = set()

    lines = [
        "@<TRIPOS>MOLECULE",
        ligand_name,
        f"{mol.GetNumAtoms()} {mol.GetNumBonds()} 1 0 0",
        "SMALL",
        "USER_CHARGES",
        "",
        "@<TRIPOS>ATOM",
    ]

    for atom in mol.GetAtoms():
        position = conformer.GetAtomPosition(atom.GetIdx())
        atom_name = _mol2_atom_name(atom, used_names)
        atom_type = _mol2_atom_type(atom)
        charge = _mol2_atom_charge(atom)
        lines.append(
            f"{atom.GetIdx() + 1:>7} {atom_name:<8} "
            f"{position.x:>10.4f} {position.y:>10.4f} {position.z:>10.4f} "
            f"{atom_type:<8} 1 {ligand_name:<8} {charge:>10.4f}"
        )

    lines.extend(["@<TRIPOS>BOND"])
    for bond in mol.GetBonds():
        lines.append(
            f"{bond.GetIdx() + 1:>6} {bond.GetBeginAtomIdx() + 1:>4} "
            f"{bond.GetEndAtomIdx() + 1:>4} {_mol2_bond_type(bond)}"
        )

    lines.extend(["@<TRIPOS>SUBSTRUCTURE", f"     1 {ligand_name:<8}     1 TEMP              0 ****  ****    0 ROOT"])

    output_mol2.parent.mkdir(parents=True, exist_ok=True)
    output_mol2.write_text("\n".join(lines) + "\n")
    return output_mol2


def append_ligand_report(report_path: Path, ligand_report: dict[str, Any]) -> Path:
    keep_all_files = ligand_report.get("keep_all_files", True)
    template_sdf = ligand_report.get("template_sdf")
    template_source = ligand_report["template_source"]
    template_sdf_display = template_sdf if keep_all_files or template_source == "user-provided" else "not kept (temporary)"

    lines = [
        "",
        "Ligand Preparation",
        "------------------",
        f"Ligand residue name: {ligand_report.get('resname') or 'not specified'}",
        f"Ligand input PDB: {ligand_report['input_pdb'] if keep_all_files else 'not kept (temporary)'}",
        f"Ligand SDF: {ligand_report['output_sdf'] if keep_all_files else 'not kept (temporary)'}",
        f"Reference ligand MOL2: {ligand_report['output_reference_mol2'] if ligand_report.get('output_reference_mol2') else 'not written'}",
        f"Ligand MOL2: {ligand_report['output_mol2'] if ligand_report.get('output_mol2') else 'not requested'}",
        f"Template bond orders: {ligand_report['used_template_bond_orders']}",
        f"Template source: {template_source}",
        f"Template SDF: {template_sdf_display if template_sdf else 'not used'}",
        f"Formal charge after template assignment: {ligand_report['formal_charge']}",
        "Reference ligand center of mass: "
        f"{ligand_report['reference_center_of_mass'][0]:.3f}, "
        f"{ligand_report['reference_center_of_mass'][1]:.3f}, "
        f"{ligand_report['reference_center_of_mass'][2]:.3f}",
        "Suggested Vina box center: "
        f"center_x={ligand_report['reference_center_of_mass'][0]:.3f}, "
        f"center_y={ligand_report['reference_center_of_mass'][1]:.3f}, "
        f"center_z={ligand_report['reference_center_of_mass'][2]:.3f}",
        f"Heavy atoms: {ligand_report['heavy_atom_count']}",
        f"Atoms after hydrogen addition: {ligand_report['atom_count_with_hydrogens']}",
        f"Local ligand cleanup applied: {ligand_report['optimized']}",
    ]

    if ligand_report.get("output_mol2"):
        if ligand_report.get("prepared_ligand_mol2"):
            lines.append(f"Prepared ligand MOL2 source: {ligand_report['prepared_ligand_mol2']}")
        lines.append(f"MOL2 charge model: {ligand_report['mol2_charge_model']}")
        charge_sum = ligand_report.get("mol2_charge_sum")
        lines.append(f"MOL2 charge sum: {charge_sum:.4f}" if charge_sum is not None else "MOL2 charge sum: unavailable")

    lines.extend(
        [
            "",
            "Note: ligand preparation here is intentionally lightweight. If the reported formal charge or "
            "protonation state is not appropriate for your docking job, prepare the ligand separately or provide "
            "a guide SDF with --ccd-template-sdf so the script uses your desired ligand chemistry instead of the "
            "default CCD template.",
        ]
    )

    with report_path.open("a") as handle:
        handle.write("\n".join(lines) + "\n")

    return report_path


def append_vina_report(
    report_path: Path,
    receptor_pdb: Path,
    receptor_pdbqt: Path,
    ligand_mol2: Path,
    ligand_pdbqt: Path,
    mk_prepare_receptor_command: str,
    mk_prepare_ligand_command: str,
    vina_receptor_source: str,
    mk_prepare_receptor_args: list[str],
) -> Path:
    receptor_extra_args = " ".join(shlex.quote(arg) for arg in mk_prepare_receptor_args)
    lines = [
        "",
        "Vina Outputs",
        "------------",
        f"Input receptor PDB: {receptor_pdb}",
        f"Output receptor PDBQT: {receptor_pdbqt}",
        f"Requested Vina receptor source: {vina_receptor_source}",
        (
            f"Receptor conversion command: {mk_prepare_receptor_command} -i <temporary staged receptor> "
            f"--write_pdbqt {receptor_pdbqt}"
            f"{' ' + receptor_extra_args if receptor_extra_args else ''}"
        ),
        (
            "Receptor conversion note: hydrogens, HETATM records, alternate locations, and CONECT records are "
            "stripped from the temporary receptor PDB before Meeko conversion. In auto mode, ANTIOCH-VS retries "
            "with a clean crystal-derived receptor if the PDBFixer-prepared receptor fails RDKit sanitization."
        ),
        f"Input ligand MOL2: {ligand_mol2}",
        f"Output ligand PDBQT: {ligand_pdbqt}",
        f"Ligand conversion command: {mk_prepare_ligand_command} -i {ligand_mol2} -o {ligand_pdbqt}",
    ]

    with report_path.open("a") as handle:
        handle.write("\n".join(lines) + "\n")

    return report_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fetch or read a PDB complex, prepare the receptor, and extract ligand files."
    )

    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--pdb-id", help="PDB ID to download from RCSB, for example 1XGJ.")
    source.add_argument("--input-pdb", help="Local complex PDB file to use instead of downloading.")

    parser.add_argument(
        "--receptor-chains",
        nargs="+",
        default=None,
        help="Optional receptor chain ID(s) to keep, for example --receptor-chains A B. Defaults to all protein chains.",
    )
    parser.add_argument(
        "--remove-receptor-residues",
        nargs="+",
        default=None,
        help=(
            "Optional protein residue ranges to remove after receptor preparation and chain selection. "
            "Use START-END for all final receptor chains, or CHAIN:START-END for a specific chain. "
            "Example: --remove-receptor-residues B:100-200 B:250-260."
        ),
    )
    parser.add_argument("--ligand-resname", required=True, help="Ligand residue name in the PDB, for example LIG.")
    parser.add_argument("--ligand-chain", default=None, help="Optional ligand chain ID.")
    parser.add_argument("--ligand-resseq", type=int, default=None, help="Optional ligand residue sequence number.")

    parser.add_argument("--output-dir", default="receptor_prep_out", help="Directory for downloaded and prepared files.")
    parser.add_argument("--prefix", default=None, help="Output prefix. Defaults to PDB ID or input filename stem.")
    parser.add_argument("--report-file", default=None, help="Optional report path. Defaults to <prefix>_prep_report.txt.")

    parser.add_argument("--keep-water", action="store_true", help="Keep water molecules during receptor preparation.")
    parser.add_argument(
        "--add-missing-residues",
        action="store_true",
        help=(
            "Allow PDBFixer to add internal missing residues. Terminal missing residues are skipped unless "
            "--add-terminal-missing-residues is also used. By default only missing atoms/hydrogens are added."
        ),
    )
    parser.add_argument(
        "--add-terminal-missing-residues",
        action="store_true",
        help=(
            "Also allow PDBFixer to rebuild missing residues at chain termini when --add-missing-residues is used. "
            "This is off by default because unresolved terminal tails can create poorly constrained docking inputs."
        ),
    )
    parser.add_argument("--ph", type=float, default=PH, help=f"pH used when adding hydrogens. Default: {PH}.")
    parser.add_argument(
        "--missing-residue-near-ligand-cutoff",
        type=float,
        default=MISSING_RESIDUE_NEAR_LIGAND_CUTOFF,
        help=(
            "Distance cutoff in Angstroms for flagging missing-residue segments near the selected ligand. "
            f"Default: {MISSING_RESIDUE_NEAR_LIGAND_CUTOFF}."
        ),
    )
    parser.add_argument(
        "--missing-residue-warning-threshold",
        type=int,
        default=MISSING_RESIDUE_WARNING_THRESHOLD,
        help=(
            "Warn when at least this many missing residues are near the selected ligand. "
            f"Default: {MISSING_RESIDUE_WARNING_THRESHOLD}."
        ),
    )
    parser.add_argument(
        "--no-optimize-ligand",
        action="store_true",
        help="Write ligand SDF without local UFF cleanup of added hydrogens.",
    )
    parser.add_argument(
        "--docking-software",
        choices=["generic", "dock6", "vina"],
        default="generic",
        help=(
            "Choose workflow-specific docking outputs. dock6 writes a no-H reference MOL2; "
            "vina writes receptor and ligand PDBQT files with Meeko. Default: generic."
        ),
    )
    parser.add_argument(
        "--write-prepared-complex",
        action="store_true",
        help="Also write a prepared receptor plus original ligand-pose PDB for visual inspection.",
    )
    parser.add_argument(
        "--keep-all-files",
        action="store_true",
        help="Keep intermediate complex, all-chain receptor, ligand PDB/SDF, and fetched CCD template files.",
    )
    parser.add_argument(
        "--write-ligand-mol2",
        action="store_true",
        help="Also write the extracted ligand as MOL2 using the corrected bond orders and coordinates.",
    )
    parser.add_argument(
        "--prepared-ligand-mol2",
        default=None,
        help="Copy this user-prepared ligand MOL2 to the output instead of generating MOL2 with RDKit.",
    )
    parser.add_argument(
        "--ligand-mol2-charge-model",
        choices=["gasteiger", "formal"],
        default="gasteiger",
        help="Charge model for MOL2 output. Default: gasteiger.",
    )
    parser.add_argument(
        "--ccd-template-sdf",
        default=None,
        help="Optional ligand template SDF with correct bond orders. Defaults to RCSB CCD <resname>_ideal.sdf.",
    )
    parser.add_argument(
        "--no-template-bond-orders",
        action="store_true",
        help="Do not fetch/use a CCD template to assign ligand bond orders.",
    )
    parser.add_argument(
        "--mk-prepare-ligand-command",
        default="mk_prepare_ligand.py",
        help="Command used to convert the Vina target ligand to PDBQT. Default: mk_prepare_ligand.py.",
    )
    parser.add_argument(
        "--mk-prepare-receptor-command",
        default="mk_prepare_receptor.py",
        help="Command used to convert the Vina receptor to PDBQT. Default: mk_prepare_receptor.py.",
    )
    parser.add_argument(
        "--mk-prepare-receptor-args",
        default="",
        help=(
            "Additional arguments to pass to mk_prepare_receptor.py during Vina receptor PDBQT conversion. "
            "Quote the whole string, for example: --mk-prepare-receptor-args \"--allow_bad_res\"."
        ),
    )
    parser.add_argument(
        "--meeko-allow-bad-res",
        action="store_true",
        help=(
            "Pass --allow_bad_res to mk_prepare_receptor.py. This lets Meeko remove residues that do not match "
            "templates, which can be useful for crystal receptors with incomplete side chains."
        ),
    )
    parser.add_argument(
        "--vina-receptor-source",
        choices=["auto", "prepared", "crystal"],
        default="auto",
        help=(
            "Receptor coordinates to feed Meeko for Vina PDBQT conversion. "
            "'prepared' uses the PDBFixer-prepared receptor, 'crystal' uses a clean heavy-atom receptor "
            "staged directly from the input complex, and 'auto' tries prepared then crystal. Default: auto."
        ),
    )
    parser.add_argument(
        "--minimize-receptor",
        action="store_true",
        help="Run a short restrained OpenMM minimization on the prepared receptor.",
    )
    parser.add_argument(
        "--minimize-max-iterations",
        type=int,
        default=200,
        help="Maximum iterations for restrained receptor minimization. Default: 200.",
    )
    parser.add_argument(
        "--minimize-restraint-k",
        type=float,
        default=10.0,
        help="Heavy-atom positional restraint strength in kcal/mol/A^2. Default: 10.0.",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.add_terminal_missing_residues and not args.add_missing_residues:
        parser.error("--add-terminal-missing-residues requires --add-missing-residues.")

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="dockflow_receptor_prep_") as tmpdir_name:
        tmpdir = Path(tmpdir_name)

        if args.pdb_id:
            prefix = args.prefix or args.pdb_id.strip().lower()
            complex_pdb = fetch_pdb(
                args.pdb_id,
                (output_dir if args.keep_all_files else tmpdir) / f"{prefix}_complex.pdb",
            )
        else:
            complex_pdb = Path(args.input_pdb).resolve()
            if not complex_pdb.exists():
                parser.error(f"Input PDB does not exist: {complex_pdb}")
            prefix = args.prefix or complex_pdb.stem

        receptor_pdb = output_dir / f"{prefix}_receptor_prepared.pdb"
        docking_software = args.docking_software
        receptor_pdbqt = output_dir / f"{prefix}_receptor_prepared.pdbqt" if docking_software == "vina" else None
        prepared_complex_pdb = output_dir / f"{prefix}_complex_prepared.pdb" if args.write_prepared_complex else None
        ligand_pdb = (output_dir if args.keep_all_files else tmpdir) / f"{prefix}_{args.ligand_resname}_ligand.pdb"
        ligand_sdf = (output_dir if args.keep_all_files else tmpdir) / f"{prefix}_{args.ligand_resname}_ligand.sdf"
        reference_ligand_mol2 = (
            output_dir / f"{prefix}_{args.ligand_resname}_reference_ligand.mol2"
            if docking_software == "dock6"
            else None
        )
        prepared_ligand_mol2 = Path(args.prepared_ligand_mol2).resolve() if args.prepared_ligand_mol2 else None
        workflow_needs_ligand_mol2 = docking_software in {"dock6", "vina"}
        ligand_mol2 = (
            output_dir / f"{prefix}_{args.ligand_resname}_ligand.mol2"
            if args.write_ligand_mol2 or prepared_ligand_mol2 is not None or workflow_needs_ligand_mol2
            else None
        )
        ligand_pdbqt = (
            output_dir / "ligands" / f"{prefix}_{args.ligand_resname}_ligand.pdbqt"
            if docking_software == "vina"
            else None
        )
        report_path = Path(args.report_file).resolve() if args.report_file else output_dir / f"{prefix}_prep_report.txt"
        prepared_all_chains_pdb = (
            output_dir if args.keep_all_files else tmpdir
        ) / f"{prefix}_receptor_prepared_all_chains.pdb"
        ccd_template_sdf = Path(args.ccd_template_sdf).resolve() if args.ccd_template_sdf else None
        mk_prepare_receptor_args = shlex.split(args.mk_prepare_receptor_args)
        if args.meeko_allow_bad_res and "--allow_bad_res" not in mk_prepare_receptor_args:
            mk_prepare_receptor_args.append("--allow_bad_res")

        print(f"Complex PDB: {complex_pdb}")
        if args.receptor_chains:
            print(f"Preparing all receptor chains before final chain selection: {', '.join(args.receptor_chains)}")
        print(f"Preparing receptor: {receptor_pdb}")
        prepare_receptor(
            input_pdb=complex_pdb,
            output_pdb=receptor_pdb,
            receptor_chains=args.receptor_chains,
            remove_receptor_residues=args.remove_receptor_residues,
            keep_water=args.keep_water,
            add_missing_residues=args.add_missing_residues,
            add_terminal_missing_residues=args.add_terminal_missing_residues,
            ph=args.ph,
            report_path=report_path,
            ligand_resname=args.ligand_resname,
            ligand_chain=args.ligand_chain,
            ligand_resseq=args.ligand_resseq,
            missing_residue_cutoff=args.missing_residue_near_ligand_cutoff,
            missing_residue_warning_threshold=args.missing_residue_warning_threshold,
            minimize_receptor=args.minimize_receptor,
            minimize_max_iterations=args.minimize_max_iterations,
            minimize_restraint_k=args.minimize_restraint_k,
            prepared_all_chains_pdb=prepared_all_chains_pdb,
            keep_all_files=args.keep_all_files,
            input_complex_is_temporary=bool(args.pdb_id and not args.keep_all_files),
            prepared_complex_pdb=prepared_complex_pdb,
        )

        if receptor_pdbqt is not None:
            print(f"Writing Vina receptor PDBQT: {receptor_pdbqt}")
            run_mk_prepare_receptor_with_fallback(
                prepared_receptor=receptor_pdb,
                input_complex=complex_pdb,
                output_pdbqt=receptor_pdbqt,
                command=args.mk_prepare_receptor_command,
                receptor_chains=args.receptor_chains,
                vina_receptor_source=args.vina_receptor_source,
                extra_args=mk_prepare_receptor_args,
            )

        print(f"Extracting ligand PDB: {ligand_pdb}")
        extract_ligand(
            input_pdb=complex_pdb,
            output_pdb=ligand_pdb,
            resname=args.ligand_resname,
            chain_id=args.ligand_chain,
            resseq=args.ligand_resseq,
        )

        if prepared_complex_pdb is not None:
            print(f"Writing prepared complex PDB: {prepared_complex_pdb}")
            write_prepared_complex_pdb(receptor_pdb, ligand_pdb, prepared_complex_pdb)

        print(f"Writing ligand SDF: {ligand_sdf}")
        if reference_ligand_mol2 is not None:
            print(f"Writing DOCK6 reference ligand MOL2: {reference_ligand_mol2}")
        if ligand_mol2 is not None:
            if prepared_ligand_mol2 is not None:
                print(f"Copying prepared ligand MOL2: {prepared_ligand_mol2} -> {ligand_mol2}")
            else:
                print(f"Writing ligand MOL2: {ligand_mol2} ({args.ligand_mol2_charge_model} charges)")
        ligand_report = ligand_pdb_to_sdf(
            input_pdb=ligand_pdb,
            output_sdf=ligand_sdf,
            resname=args.ligand_resname,
            template_sdf=ccd_template_sdf,
            use_template_bond_orders=not args.no_template_bond_orders,
            optimize=not args.no_optimize_ligand,
            output_reference_mol2=reference_ligand_mol2,
            output_mol2=ligand_mol2,
            mol2_charge_model=args.ligand_mol2_charge_model,
            prepared_ligand_mol2=prepared_ligand_mol2,
            report_path=report_path,
            keep_all_files=args.keep_all_files,
        )

        if ligand_pdbqt is not None:
            if ligand_mol2 is None:
                raise RuntimeError("Internal error: Vina workflow requires a ligand MOL2 input for PDBQT conversion.")
            print(f"Writing Vina ligand PDBQT: {ligand_pdbqt}")
            run_mk_prepare_ligand(ligand_mol2, ligand_pdbqt, args.mk_prepare_ligand_command)
            if receptor_pdbqt is None:
                raise RuntimeError("Internal error: Vina workflow requires a receptor PDBQT output.")
            append_vina_report(
                report_path,
                receptor_pdb,
                receptor_pdbqt,
                ligand_mol2,
                ligand_pdbqt,
                args.mk_prepare_receptor_command,
                args.mk_prepare_ligand_command,
                args.vina_receptor_source,
                mk_prepare_receptor_args,
            )

        print("\nDone.")
        print(f"Prepared receptor: {receptor_pdb}")
        if receptor_pdbqt is not None:
            print(f"Receptor PDBQT:    {receptor_pdbqt}")
        if prepared_complex_pdb is not None:
            print(f"Prepared complex:  {prepared_complex_pdb}")
        if args.keep_all_files:
            if args.pdb_id:
                print(f"Complex PDB:       {complex_pdb}")
            if args.receptor_chains:
                print(f"All-chain receptor: {prepared_all_chains_pdb}")
            print(f"Ligand PDB:        {ligand_pdb}")
            print(f"Ligand SDF:        {ligand_sdf}")
        center = ligand_report["reference_center_of_mass"]
        print(f"Reference center:  {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}")
        if reference_ligand_mol2 is not None:
            print(f"Reference MOL2:    {reference_ligand_mol2}")
        if ligand_mol2 is not None:
            print(f"Ligand MOL2:       {ligand_mol2}")
        if ligand_pdbqt is not None:
            print(f"Ligand PDBQT:      {ligand_pdbqt}")
        print(f"Preparation report: {report_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
