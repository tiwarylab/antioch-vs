#!/usr/bin/env python3
"""Separate prepared ligand PDBQTs by Vina-GPU and Vina CPU parser support.

The helper inspects the PDBQT atom-type field; it does not alter molecular
structures or atom types.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


# This is the complete AutoDock atom-kind table in the official Vina-GPU 2.1
# source, plus its documented Se -> S alias.
VINA_GPU_2_1_LEGACY_TYPES = frozenset(
    {
        "A",
        "Br",
        "C",
        "Ca",
        "Cl",
        "F",
        "Fe",
        "H",
        "HD",
        "I",
        "Mg",
        "Mn",
        "N",
        "NA",
        "O",
        "OA",
        "P",
        "S",
        "SA",
        "Se",
        "Zn",
    }
)

# This is the complete atom-kind table in current AutoDock Vina, plus its
# documented Se -> S alias. Atom types outside this set are not routed to
# either engine by default. In particular, boron (B) is not a Vina atom type.
VINA_CPU_CURRENT_TYPES = frozenset(
    {
        "A",
        "At",
        "Br",
        "C",
        "Ca",
        "CG0",
        "CG1",
        "CG2",
        "CG3",
        "Cl",
        "F",
        "Fe",
        "G0",
        "G1",
        "G2",
        "G3",
        "H",
        "HD",
        "I",
        "Mg",
        "Mn",
        "N",
        "NA",
        "O",
        "OA",
        "P",
        "S",
        "SA",
        "Se",
        "Si",
        "W",
        "Zn",
    }
)


@dataclass(frozen=True)
class Classification:
    """One input PDBQT and its engine compatibility decision."""

    source: Path
    atom_types: tuple[str, ...]
    unsupported_gpu_types: tuple[str, ...]
    unsupported_cpu_types: tuple[str, ...]
    atom_record_count: int
    status: str
    detail: str


def parse_pdbqt_atom_types(path: Path) -> tuple[tuple[str, ...], int]:
    """Return unique PDBQT atom types and the number of ATOM/HETATM records."""
    atom_types: set[str] = set()
    atom_records = 0

    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(("ATOM", "HETATM")):
                continue

            atom_records += 1
            # Atom type is the final whitespace-delimited PDBQT token.  This
            # intentionally avoids the historical two-character fixed-width
            # assumption: newer Meeko macrocycle types such as CG0 are wider.
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: malformed PDBQT atom record.")
            atom_type = fields[-1]
            atom_types.add(atom_type)

    if atom_records == 0:
        raise ValueError(f"{path}: no ATOM or HETATM records found.")
    return tuple(sorted(atom_types)), atom_records


def classify_pdbqt(
    path: Path,
    gpu_allowed_types: frozenset[str],
    cpu_allowed_types: frozenset[str],
) -> Classification:
    """Classify one PDBQT without changing it."""
    try:
        atom_types, atom_records = parse_pdbqt_atom_types(path)
    except (OSError, ValueError) as exc:
        return Classification(
            source=path,
            atom_types=(),
            unsupported_gpu_types=(),
            unsupported_cpu_types=(),
            atom_record_count=0,
            status="invalid_pdbqt",
            detail=str(exc),
        )

    unsupported_gpu = tuple(atom_type for atom_type in atom_types if atom_type not in gpu_allowed_types)
    unsupported_cpu = tuple(atom_type for atom_type in atom_types if atom_type not in cpu_allowed_types)
    if unsupported_cpu:
        return Classification(
            source=path,
            atom_types=atom_types,
            unsupported_gpu_types=unsupported_gpu,
            unsupported_cpu_types=unsupported_cpu,
            atom_record_count=atom_records,
            status="vina_unsupported",
            detail="Vina CPU profile does not accept: " + ", ".join(unsupported_cpu),
        )
    if unsupported_gpu:
        return Classification(
            source=path,
            atom_types=atom_types,
            unsupported_gpu_types=unsupported_gpu,
            unsupported_cpu_types=(),
            atom_record_count=atom_records,
            status="cpu_only",
            detail="Vina-GPU profile does not accept: " + ", ".join(unsupported_gpu),
        )
    return Classification(
        source=path,
        atom_types=atom_types,
        unsupported_gpu_types=(),
        unsupported_cpu_types=(),
        atom_record_count=atom_records,
        status="gpu_compatible",
        detail="All atom types are accepted by the Vina-GPU profile.",
    )


def collect_inputs(input_path: Path) -> list[Path]:
    """Collect one PDBQT or direct PDBQT children of a ligand directory."""
    if input_path.is_file():
        if input_path.suffix.lower() != ".pdbqt":
            raise ValueError(f"Input file must end in .pdbqt: {input_path}")
        return [input_path]
    if not input_path.is_dir():
        raise FileNotFoundError(f"Input path does not exist: {input_path}")

    files = sorted(path for path in input_path.iterdir() if path.is_file() and path.suffix.lower() == ".pdbqt")
    if not files:
        raise ValueError(f"No direct .pdbqt files found in input directory: {input_path}")
    return files


def materialize(source: Path, destination: Path, mode: str) -> None:
    """Create a standalone output-file entry using the requested materialization mode."""
    if mode == "copy":
        shutil.copy2(source, destination)
    elif mode == "hardlink":
        os.link(source, destination)
    elif mode == "symlink":
        destination.symlink_to(source.resolve())
    else:  # argparse prevents this branch.
        raise ValueError(f"Unsupported materialization mode: {mode}")


def write_outputs(
    classifications: list[Classification],
    output_dir: Path,
    materialize_mode: str,
    cpu_ligand_set: str,
    input_path: Path,
    gpu_allowed_types: frozenset[str],
    cpu_allowed_types: frozenset[str],
) -> None:
    """Write separated input directories plus CSV and JSON provenance manifests."""
    gpu_dir = output_dir / "vina_gpu_compatible"
    cpu_dir = output_dir / "vina_cpu_compatible"
    unsupported_dir = output_dir / "vina_unsupported"
    invalid_dir = output_dir / "invalid_pdbqt"
    for directory in (gpu_dir, cpu_dir, unsupported_dir, invalid_dir):
        directory.mkdir(parents=True, exist_ok=False)

    manifest_rows: list[dict[str, object]] = []
    for item in classifications:
        filename = item.source.name
        destinations: list[Path] = []
        if item.status == "gpu_compatible":
            destinations = [gpu_dir / filename]
            if cpu_ligand_set == "all_compatible":
                destinations.append(cpu_dir / filename)
        elif item.status == "cpu_only":
            destinations = [cpu_dir / filename]
        elif item.status == "vina_unsupported":
            destinations = [unsupported_dir / filename]
        else:
            destinations = [invalid_dir / filename]

        for destination in destinations:
            materialize(item.source, destination, materialize_mode)

        manifest_rows.append(
            {
                "source_path": str(item.source),
                "filename": filename,
                "status": item.status,
                "atom_record_count": item.atom_record_count,
                "atom_types": ";".join(item.atom_types),
                "unsupported_gpu_types": ";".join(item.unsupported_gpu_types),
                "unsupported_cpu_types": ";".join(item.unsupported_cpu_types),
                "engine_destinations": ";".join(destination.parent.name for destination in destinations),
                "detail": item.detail,
            }
        )

    manifest_csv = output_dir / "vina_ligand_compatibility.csv"
    with manifest_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]))
        writer.writeheader()
        writer.writerows(manifest_rows)

    counts = Counter(item.status for item in classifications)
    digest = hashlib.sha256(
        "\n".join(str(item.source) for item in classifications).encode("utf-8")
    ).hexdigest()
    summary = {
        "schema_version": 2,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_path": str(input_path),
        "input_file_order_sha256": digest,
        "gpu_profile": "vina_gpu_2_1_legacy",
        "gpu_allowed_types": sorted(gpu_allowed_types),
        "cpu_profile": "autodock_vina_current",
        "cpu_allowed_types": sorted(cpu_allowed_types),
        "materialize_mode": materialize_mode,
        "cpu_ligand_set": cpu_ligand_set,
        "counts": {
            "input_pdbqt_files": len(classifications),
            "vina_gpu_compatible": counts["gpu_compatible"],
            "vina_cpu_compatible": counts["gpu_compatible"] + counts["cpu_only"],
            "cpu_only": counts["cpu_only"],
            "vina_unsupported": counts["vina_unsupported"],
            "invalid_pdbqt": counts["invalid_pdbqt"],
            "vina_gpu_input": counts["gpu_compatible"],
            "vina_cpu_input": counts["gpu_compatible"] + counts["cpu_only"]
            if cpu_ligand_set == "all_compatible"
            else counts["cpu_only"],
            "gpu_cpu_overlap": counts["gpu_compatible"] if cpu_ligand_set == "all_compatible" else 0,
        },
        "outputs": {
            "vina_gpu_compatible_dir": str(gpu_dir),
            "vina_cpu_compatible_dir": str(cpu_dir),
            "vina_unsupported_dir": str(unsupported_dir),
            "invalid_pdbqt_dir": str(invalid_dir),
            "manifest_csv": str(manifest_csv),
        },
    }
    (output_dir / "vina_ligand_compatibility_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=False) + "\n", encoding="utf-8"
    )


def parse_allowed_types(value: str) -> frozenset[str]:
    """Parse a comma-separated explicit Vina atom-type allowlist."""
    parsed = frozenset(part.strip() for part in value.split(",") if part.strip())
    if not parsed:
        raise argparse.ArgumentTypeError("An atom-type allowlist must contain at least one atom type.")
    return parsed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Prepared ligand .pdbqt file or directory.")
    parser.add_argument("--output-dir", required=True, type=Path, help="New directory for separated ligand inputs and manifests.")
    parser.add_argument(
        "--materialize",
        choices=("copy", "hardlink", "symlink"),
        default="copy",
        help="How to place input files in output directories (default: copy).",
    )
    parser.add_argument(
        "--cpu-allowed-types",
        type=parse_allowed_types,
        default=VINA_CPU_CURRENT_TYPES,
        help=(
            "Comma-separated atom types accepted by the installed Vina CPU build. "
            "Default: current AutoDock Vina profile."
        ),
    )
    parser.add_argument(
        "--gpu-allowed-types",
        type=parse_allowed_types,
        default=VINA_GPU_2_1_LEGACY_TYPES,
        help=(
            "Comma-separated atom types accepted by the installed Vina-GPU build. "
            "Default: legacy Vina-GPU 2.1 profile."
        ),
    )
    parser.add_argument(
        "--cpu-ligand-set",
        choices=("all_compatible", "cpu_only"),
        default="all_compatible",
        help=(
            "Whether Vina CPU also receives GPU-compatible PDBQTs "
            "(all_compatible, the default) or only CPU-only PDBQTs "
            "(cpu_only, no GPU/CPU overlap)."
        ),
    )
    parser.add_argument("--overwrite", action="store_true", help="Remove an existing --output-dir before writing.")
    args = parser.parse_args(argv)

    input_path = args.input.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if input_path.is_dir() and output_dir.is_relative_to(input_path):
        parser.error("--output-dir must not be inside an input ligand directory.")
    if output_dir.exists():
        if not args.overwrite:
            parser.error(f"Output directory already exists: {output_dir}. Use --overwrite to replace it.")
        shutil.rmtree(output_dir)

    inputs = collect_inputs(input_path)
    classifications = [
        classify_pdbqt(path, args.gpu_allowed_types, args.cpu_allowed_types) for path in inputs
    ]
    write_outputs(
        classifications,
        output_dir,
        args.materialize,
        args.cpu_ligand_set,
        input_path,
        args.gpu_allowed_types,
        args.cpu_allowed_types,
    )

    counts = Counter(item.status for item in classifications)
    print(f"Input PDBQTs: {len(classifications)}")
    print(f"Vina-GPU compatible: {counts['gpu_compatible']}")
    print(f"CPU only: {counts['cpu_only']}")
    print(f"Unsupported by configured Vina profiles: {counts['vina_unsupported']}")
    print(f"Invalid PDBQT: {counts['invalid_pdbqt']}")
    print(f"Vina-GPU ligand directory: {output_dir / 'vina_gpu_compatible'}")
    print(f"Vina CPU ligand directory: {output_dir / 'vina_cpu_compatible'}")
    print(f"CPU ligand set: {args.cpu_ligand_set}")
    print(f"Manifest: {output_dir / 'vina_ligand_compatibility.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
