#!/usr/bin/env python3
"""
relabel_only_protein.py — restore original residue numbering on a truncated system's topology PDB.

When a receptor is truncated before MD (a loop excised + ACE/NME caps), tleap renumbers the topology
sequentially 1..N, collapsing the gap. The per-folder only_protein.pdb therefore no longer matches the
AF2RAVE feature CSV's original numbering, so the COLVAR feature names fail to map. This relabels it back
(atoms/coords untouched, so the matching XTC still aligns) and renames AMBER CYX -> standard CYS
(the COLVAR mapper expands HIS variants but not CYS).

Mapping for a single excised range: residues below the cut keep their number; the two caps sit in the
gap; residues at/after `--split` get `+offset`.

    python3 relabel_only_protein.py only_protein.pdb --split 202 --offset 30 [--keep-cyx]
"""
import argparse
import MDAnalysis as mda

a = argparse.ArgumentParser()
a.add_argument("pdb")
a.add_argument("--split", type=int, required=True, help="first sequential resid of the post-cut segment")
a.add_argument("--offset", type=int, required=True, help="added to resid >= split")
a.add_argument("--keep-cyx", action="store_true", help="do not rename CYX -> CYS")
a.add_argument("--max-before", type=int, default=100000, help="guard: skip if a non-cap resid already exceeds this")
args = a.parse_args()

u = mda.Universe(args.pdb)
caps = {"ACE", "NME", "NHE"}
noncap = [int(r.resid) for r in u.residues if r.resname not in caps]
# Guard against a double-apply: once relabeled, post-cut residues exceed the pre-relabel max. Set
# --max-before to the pre-relabel residue count to make this strict; left high it is effectively off.
if max(noncap) > args.max_before:
    print(f"[skip] {args.pdb}: max non-cap resid {max(noncap)} > {args.max_before} — already relabeled")
    raise SystemExit(0)

resids = u.residues.resids.copy()
for i, r in enumerate(u.residues):
    if r.resname not in caps and int(resids[i]) >= args.split:
        resids[i] = int(resids[i]) + args.offset
u.residues.resids = resids

n_cyx = 0
if not args.keep_cyx:
    cyx = u.select_atoms("resname CYX")
    n_cyx = len(cyx.residues)
    if n_cyx:
        cyx.residues.resnames = ["CYS"] * n_cyx

u.atoms.write(args.pdb)
print(f"relabeled {args.pdb}: +{args.offset} for resid>={args.split}, {n_cyx} CYX->CYS")
