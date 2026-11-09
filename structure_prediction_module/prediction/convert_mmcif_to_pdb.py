#!/usr/bin/env python3
"""Convert AlphaFold-style mmCIF files to PDB with Biopython."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def safe_stem(path: Path, root: Path) -> str:
    relative = path.relative_to(root).with_suffix("")
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", "__".join(relative.parts))


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert all mmCIF files below a directory to PDB.")
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    try:
        from Bio.PDB import PDBIO, MMCIFParser
    except ImportError as exc:
        raise ImportError(
            "mmCIF conversion requires Biopython. Install biopython in the Python environment used for this runner."
        ) from exc

    input_dir = args.input_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cif_files = sorted(input_dir.rglob("*.cif"))
    if not cif_files:
        raise FileNotFoundError(f"No mmCIF files found in {input_dir}")

    parser_cif = MMCIFParser(QUIET=True)
    for index, cif_file in enumerate(cif_files, start=1):
        structure = parser_cif.get_structure(f"model_{index}", str(cif_file))
        target = output_dir / f"{safe_stem(cif_file, input_dir)}.pdb"
        writer = PDBIO()
        writer.set_structure(structure)
        writer.save(str(target))
        print(f"Converted: {cif_file} -> {target}")


if __name__ == "__main__":
    main()
