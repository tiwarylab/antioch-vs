#!/usr/bin/env python3
"""Run PROPKA and make a reviewed AMBER residue-state proposal.

The program writes a renamed PDB, its PROPKA report, and a per-residue
decision table. It deliberately does not guess neutral histidine tautomers
or emit LEaP bond commands; CPPTRAJ prepareforleap is the single downstream
source of disulfide connectivity.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


ACIDIC_RENAMES = {"ASP": "ASH", "GLU": "GLH"}
TITRATABLE_RESIDUES = {"ASP", "GLU", "HIS", "CYS", "LYS", "ARG", "TYR"}
SUMMARY_PATTERN = re.compile(
    r"^\s*(ASP|GLU|HIS|CYS|LYS|ARG|TYR)\s+(-?\d+)\s+(\S)\s+(-?\d+(?:\.\d+)?)\b"
)


@dataclass(frozen=True)
class ResidueKey:
    chain: str
    resseq: int
    icode: str = ""

    @property
    def label(self) -> str:
        return f"{self.chain}:{self.resseq}{self.icode}" if self.chain else f":{self.resseq}{self.icode}"


@dataclass(frozen=True)
class PkaValue:
    residue_name: str
    pka: float


def parse_residue_label(text: str) -> ResidueKey:
    """Parse CHAIN:RESIDUE or CHAIN:RESIDUEINSERTION, for example A:48."""
    match = re.fullmatch(r"([^:]*):(-?\d+)([A-Za-z]?)", text.strip())
    if not match:
        raise argparse.ArgumentTypeError(
            f"Invalid residue label {text!r}; expected CHAIN:RESIDUE, for example A:48."
        )
    return ResidueKey(match.group(1), int(match.group(2)), match.group(3))


def parse_override(text: str) -> tuple[ResidueKey, str]:
    try:
        label, residue_name = text.split("=", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid override {text!r}; expected CHAIN:RESIDUE=AMBER_NAME, for example A:48=ASP."
        ) from exc
    residue_name = residue_name.upper()
    if not re.fullmatch(r"[A-Z0-9]{3}", residue_name):
        raise argparse.ArgumentTypeError(f"Invalid AMBER residue name {residue_name!r}.")
    return parse_residue_label(label), residue_name


def pka_output_path(pdb_path: Path) -> Path:
    return pdb_path.with_suffix(".pka")


def run_propka(pdb_path: Path, executable: str, reuse_existing: bool) -> Path:
    report = pka_output_path(pdb_path)
    if reuse_existing:
        if not report.is_file():
            raise FileNotFoundError(f"--reuse-existing-pka was requested, but {report} does not exist.")
        return report

    try:
        result = subprocess.run(
            [executable, pdb_path.name],
            cwd=pdb_path.parent,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"Could not find {executable!r}. Activate the PROPKA environment or pass --propka-exe."
        ) from exc

    if result.returncode != 0:
        raise RuntimeError(f"PROPKA failed for {pdb_path}:\n{result.stdout}")
    if not report.is_file():
        raise RuntimeError(
            f"PROPKA returned success for {pdb_path}, but did not create the expected report {report}.\n"
            f"Output:\n{result.stdout}"
        )
    return report


def parse_pka_summary(report_path: Path) -> dict[ResidueKey, PkaValue]:
    """Read only the stable SUMMARY OF THIS PREDICTION section."""
    values: dict[ResidueKey, PkaValue] = {}
    in_summary = False
    for line in report_path.read_text(errors="replace").splitlines():
        if line.strip() == "SUMMARY OF THIS PREDICTION":
            in_summary = True
            continue
        if not in_summary:
            continue
        match = SUMMARY_PATTERN.match(line)
        if match:
            residue_name, resseq, chain, value = match.groups()
            values[ResidueKey(chain, int(resseq))] = PkaValue(residue_name, float(value))
    if not values:
        raise ValueError(f"No protein pKa values found in {report_path}.")
    return values


def pdb_residues(pdb_path: Path) -> dict[ResidueKey, str]:
    residues: dict[ResidueKey, str] = {}
    for line in pdb_path.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 27:
            key = ResidueKey(line[21].strip(), int(line[22:26]), line[26].strip())
            residues.setdefault(key, line[17:20].strip().upper())
    if not residues:
        raise ValueError(f"No ATOM/HETATM residue records found in {pdb_path}.")
    return residues


def sg_coordinates(pdb_path: Path) -> dict[ResidueKey, tuple[float, float, float]]:
    coordinates: dict[ResidueKey, tuple[float, float, float]] = {}
    for line in pdb_path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        if line[12:16].strip() != "SG":
            continue
        key = ResidueKey(line[21].strip(), int(line[22:26]), line[26].strip())
        coordinates[key] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return coordinates


def find_disulfide_pairs(
    coordinates: dict[ResidueKey, tuple[float, float, float]],
    candidates: set[ResidueKey],
    max_distance: float,
) -> tuple[list[tuple[ResidueKey, ResidueKey, float]], list[str]]:
    """Pair only isolated close SG-SG candidates; ambiguous contacts are withheld."""
    close: dict[ResidueKey, list[tuple[ResidueKey, float]]] = {key: [] for key in candidates if key in coordinates}
    keys = sorted(close, key=lambda key: (key.chain, key.resseq, key.icode))
    for index, first in enumerate(keys):
        for second in keys[index + 1 :]:
            distance = math.dist(coordinates[first], coordinates[second])
            if distance <= max_distance:
                close[first].append((second, distance))
                close[second].append((first, distance))

    pairs: list[tuple[ResidueKey, ResidueKey, float]] = []
    warnings: list[str] = []
    seen: set[ResidueKey] = set()
    for first in keys:
        if first in seen:
            continue
        neighbours = close[first]
        if len(neighbours) == 0:
            warnings.append(f"{first.label}: PROPKA disulfide candidate has no SG partner within {max_distance:.2f} A.")
            continue
        if len(neighbours) != 1:
            warnings.append(f"{first.label}: ambiguous SG partners; no CYX rename was applied.")
            continue
        second, distance = neighbours[0]
        if len(close[second]) != 1:
            warnings.append(f"{first.label}/{second.label}: ambiguous SG pairing; no CYX rename was applied.")
            continue
        pairs.append((first, second, distance))
        seen.update({first, second})
    return pairs, warnings


def decide_states(
    residues: dict[ResidueKey, str],
    pka_values: dict[ResidueKey, PkaValue],
    ph: float,
    cutoff: float,
    cysteine_mode: str,
    disulfide_pka: float,
    disulfide_distance: float,
    coordinates: dict[ResidueKey, tuple[float, float, float]],
    protected: set[ResidueKey],
    overrides: dict[ResidueKey, str],
) -> tuple[dict[ResidueKey, str], list[dict[str, object]], list[tuple[ResidueKey, ResidueKey, float]], list[str]]:
    decisions: dict[ResidueKey, str] = {}
    rows: list[dict[str, object]] = []
    warnings: list[str] = []

    cysteine_candidates = {
        key
        for key, value in pka_values.items()
        if residues.get(key) == "CYS" and value.pka >= disulfide_pka and key not in protected and key not in overrides
    }
    disulfides: list[tuple[ResidueKey, ResidueKey, float]] = []
    if cysteine_mode == "propka":
        disulfides, warnings = find_disulfide_pairs(coordinates, cysteine_candidates, disulfide_distance)
        for first, second, _ in disulfides:
            decisions[first] = "CYX"
            decisions[second] = "CYX"

    for key, original_name in residues.items():
        value = pka_values.get(key)
        proposed_name = original_name
        reason = "No PROPKA pKa available; unchanged."
        action = "keep"

        if key in overrides:
            proposed_name = overrides[key]
            reason = "Explicit user override."
            action = "override"
        elif key in protected:
            reason = "Protected from automatic renaming."
            action = "protected"
        elif key in decisions:
            proposed_name = decisions[key]
            reason = f"PROPKA disulfide candidate with a unique SG partner within {disulfide_distance:.2f} A."
            action = "rename"
        elif value and original_name in ACIDIC_RENAMES:
            if value.pka >= ph + cutoff:
                proposed_name = ACIDIC_RENAMES[original_name]
                reason = f"pKa {value.pka:.2f} >= pH {ph:.2f} + cutoff {cutoff:.2f}."
                action = "rename"
            else:
                reason = f"pKa {value.pka:.2f} below protonation threshold {ph + cutoff:.2f}."
        elif value and original_name == "HIS":
            if value.pka >= ph + cutoff:
                proposed_name = "HIP"
                reason = f"pKa {value.pka:.2f} >= pH {ph:.2f} + cutoff {cutoff:.2f}."
                action = "rename"
            else:
                reason = (
                    f"pKa {value.pka:.2f}; neutral His tautomer cannot be derived from PROPKA. "
                    "Review as HID or HIE explicitly."
                )
                action = "review"
        elif value and original_name == "CYS":
            if cysteine_mode == "keep":
                reason = "Cysteine mode is keep; unchanged."
            elif value.pka >= disulfide_pka:
                reason = "PROPKA disulfide candidate was not paired unambiguously; unchanged."
                action = "review"
            elif value.pka <= ph - cutoff:
                proposed_name = "CYM"
                reason = f"pKa {value.pka:.2f} <= pH {ph:.2f} - cutoff {cutoff:.2f}."
                action = "rename"
            else:
                reason = f"pKa {value.pka:.2f} does not cross the deprotonation threshold {ph - cutoff:.2f}."
        elif value and original_name in {"LYS", "ARG", "TYR"}:
            reason = (
                f"pKa {value.pka:.2f}; automatic nonstandard {original_name} states are intentionally disabled."
            )
            action = "review"

        if proposed_name != original_name:
            decisions[key] = proposed_name
        rows.append(
            {
                "chain": key.chain,
                "resseq": key.resseq,
                "insertion_code": key.icode,
                "original_residue": original_name,
                "propka_residue": value.residue_name if value else "",
                "pka": f"{value.pka:.2f}" if value else "",
                "proposed_amber_residue": proposed_name,
                "action": action,
                "reason": reason,
            }
        )
    return decisions, rows, disulfides, warnings


def write_renamed_pdb(input_path: Path, output_path: Path, decisions: dict[ResidueKey, str]) -> None:
    output_lines: list[str] = []
    for line in input_path.read_text(errors="replace").splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 27:
            key = ResidueKey(line[21].strip(), int(line[22:26]), line[26].strip())
            if key in decisions:
                line = f"{line[:17]}{decisions[key]:>3s}{line[20:]}"
        output_lines.append(line)
    output_path.write_text("".join(output_lines))


def write_outputs(
    output_dir: Path,
    pdb_path: Path,
    report_path: Path,
    decisions: dict[ResidueKey, str],
    rows: list[dict[str, object]],
    disulfides: list[tuple[ResidueKey, ResidueKey, float]],
    warnings: list[str],
    dry_run: bool,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = pdb_path.stem
    if not dry_run:
        write_renamed_pdb(pdb_path, output_dir / f"{stem}_amber_states.pdb", decisions)
        shutil.copy2(report_path, output_dir / report_path.name)

    with (output_dir / f"{stem}_protonation_decisions.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    # The legacy PROPKA bond fragment is intentionally retired. CPPTRAJ
    # prepares and reports the final disulfide connectivity downstream.
    (output_dir / f"{stem}_disulfide_bonds.leap").unlink(missing_ok=True)

    manifest = {
        "input_pdb": str(pdb_path.resolve()),
        "propka_report": str(report_path.resolve()),
        "renamed_pdb": None if dry_run else str((output_dir / f"{stem}_amber_states.pdb").resolve()),
        "renamed_residues": {key.label: value for key, value in sorted(decisions.items(), key=lambda item: item[0].label)},
        "disulfides": [
            {"first": first.label, "second": second.label, "sg_distance_angstrom": round(distance, 3)}
            for first, second, distance in disulfides
        ],
        "warnings": warnings,
    }
    (output_dir / f"{stem}_protonation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def input_pdbs(input_path: Path, recursive: bool) -> list[Path]:
    if input_path.is_file():
        if input_path.suffix.lower() != ".pdb":
            raise ValueError(f"Input file must be a .pdb file, got {input_path}.")
        return [input_path]
    if input_path.is_dir():
        pattern = "**/*.pdb" if recursive else "*.pdb"
        files = sorted(path for path in input_path.glob(pattern) if path.is_file())
        if files:
            return files
        raise FileNotFoundError(f"No PDB files found in {input_path}.")
    raise FileNotFoundError(f"Input path does not exist: {input_path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="One PDB file or a directory containing PDB files.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Directory for prepared PDBs and reports.")
    parser.add_argument("--ph", type=float, default=7.0, help="Target pH (default: 7.0).")
    parser.add_argument(
        "--pka-cutoff",
        type=float,
        default=2.0,
        help="Required pKa-pH separation for an automatic state change (default: 2.0).",
    )
    parser.add_argument("--propka-exe", default="propka3", help="PROPKA executable (default: propka3).")
    parser.add_argument("--reuse-existing-pka", action="store_true", help="Do not rerun PROPKA; require a sibling .pka file.")
    parser.add_argument("--recursive", action="store_true", help="Search a directory recursively for PDB files.")
    parser.add_argument("--dry-run", action="store_true", help="Write reports but not the renamed PDB.")
    parser.add_argument(
        "--cysteine-mode",
        choices=("keep", "propka"),
        default="keep",
        help="keep leaves CYS unchanged; propka can create CYM/CYX states (default: keep).",
    )
    parser.add_argument(
        "--disulfide-pka",
        type=float,
        default=99.0,
        help="PROPKA pKa at or above which Cys is treated as a disulfide candidate (default: 99.0).",
    )
    parser.add_argument(
        "--disulfide-sg-cutoff",
        type=float,
        default=2.3,
        help="Maximum SG-SG distance for automatic CYX pairing in Angstrom (default: 2.3).",
    )
    parser.add_argument(
        "--protect",
        action="append",
        type=parse_residue_label,
        default=[],
        metavar="CHAIN:RESIDUE",
        help="Keep one residue unchanged despite its pKa; repeatable.",
    )
    parser.add_argument(
        "--override",
        action="append",
        type=parse_override,
        default=[],
        metavar="CHAIN:RESIDUE=NAME",
        help="Explicit AMBER residue name; repeatable, and takes precedence over --protect.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.pka_cutoff < 0:
        raise ValueError("--pka-cutoff must be non-negative.")
    if args.disulfide_sg_cutoff <= 0:
        raise ValueError("--disulfide-sg-cutoff must be positive.")

    overrides = dict(args.override)
    for pdb_path in input_pdbs(args.input.resolve(), args.recursive):
        report_path = run_propka(pdb_path, args.propka_exe, args.reuse_existing_pka)
        residues = pdb_residues(pdb_path)
        decisions, rows, disulfides, warnings = decide_states(
            residues=residues,
            pka_values=parse_pka_summary(report_path),
            ph=args.ph,
            cutoff=args.pka_cutoff,
            cysteine_mode=args.cysteine_mode,
            disulfide_pka=args.disulfide_pka,
            disulfide_distance=args.disulfide_sg_cutoff,
            coordinates=sg_coordinates(pdb_path),
            protected=set(args.protect),
            overrides=overrides,
        )
        target_dir = args.out_dir.resolve() / pdb_path.stem
        write_outputs(
            output_dir=target_dir,
            pdb_path=pdb_path,
            report_path=report_path,
            decisions=decisions,
            rows=rows,
            disulfides=disulfides,
            warnings=warnings,
            dry_run=args.dry_run,
        )
        renamed = ", ".join(f"{key.label}->{name}" for key, name in decisions.items()) or "none"
        print(f"{pdb_path.name}: proposed renames: {renamed}")
        print(f"  Output: {target_dir}")
        for warning in warnings:
            print(f"  WARNING: {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
