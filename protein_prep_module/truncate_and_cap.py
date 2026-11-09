#!/usr/bin/env python3
"""truncate_and_cap.py — excise residue range(s) from a protein and cap the new ends with ACE/NME.

All structure manipulation and I/O go through ParmEd (AmberTools' structure library), so there is no
hand-formatted PDB parsing/writing. The only modelling assumption is that each cap reuses the flanking
removed residue's own backbone atoms, which places the cap at exact peptide geometry relative to its
neighbour; the ACE/NME force-field parameters and the methyl/amide hydrogens are then supplied by ff14SB
when tleap builds the topology downstream.

Use case: a predicted structure (e.g. AF2) may model a long disordered loop (such as a GPCR ICL3) that
swings out into solvent and inflates the solvation/membrane box. When that loop is dispensable it can be
removed to shrink the box and lower the simulation cost.

For each removed range [s, e]:
  * residues s+1 .. e-1 are stripped;
  * residue s becomes NME (keep backbone N; rename CA -> CH3)   -> caps the C-terminus of residue s-1;
  * residue e becomes ACE (rename CA -> CH3; keep backbone C, O) -> caps the N-terminus of residue e+1.
The second segment is put on a fresh chain id so a TER is written between the caps (tleap then treats the
two pieces as separate molecules). Residue numbering is preserved, so partners outside the removed range
(e.g. disulfide cysteines) keep their numbers.

Usage:
    truncate_and_cap.py IN.pdb OUT.pdb --remove 200-231 [230-240 ...] [--chain A]
"""
import argparse
import string
import sys

import parmed as pmd


def parse_ranges(tokens):
    ranges = []
    for tok in tokens:
        for part in tok.replace(",", " ").split():
            if "-" not in part:
                sys.exit(f"[truncate_and_cap] bad range '{part}' (expected START-END)")
            a, b = part.split("-", 1)
            s, e = int(a), int(b)
            if e <= s + 1:
                sys.exit(f"[truncate_and_cap] range {s}-{e} too small to cap "
                         f"(need >=1 interior residue to remove)")
            ranges.append((s, e))
    return ranges


def main():
    ap = argparse.ArgumentParser(description="Excise residue range(s) and cap with ACE/NME (ParmEd).")
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("--remove", nargs="+", required=True,
                    help="inclusive START-END residue ranges to excise (e.g. 200-231)")
    ap.add_argument("--chain", default=None, help="restrict to this chain id (default: all chains)")
    args = ap.parse_args()
    ranges = parse_ranges(args.remove)

    st = pmd.load_file(args.inp)

    def eligible(res):
        return args.chain is None or res.chain == args.chain

    # Decide, per residue, what to do — set names on the kept cap atoms, collect atoms to strip.
    # (Names are set on the live objects before stripping; stripping reindexes but keeps the names.)
    del_atoms = []            # ParmEd Atom objects to remove
    report = []
    used_chains = {a.residue.chain for a in st.atoms}
    free_chains = [c for c in string.ascii_uppercase if c not in used_chains]
    seg2_chain = {}           # original chain -> new chain id for the post-cut segment

    for res in st.residues:
        if not eligible(res):
            continue
        rn = res.number
        for (s, e) in ranges:
            if s < rn < e:                                  # interior -> remove whole residue
                del_atoms.extend(res.atoms)
            elif rn == s:                                   # -> NME: keep ONLY the backbone N
                # tleap builds the NME methyl carbon (template name "C", not "CH3") and its hydrogens
                # from the ff template, placed off this N; providing our own methyl carbon would be an
                # orphan atom ("not in residue templates") that fails type assignment.
                res.name = "NME"
                for a in list(res.atoms):
                    if a.name != "N":
                        del_atoms.append(a)
            elif rn == e:                                   # -> ACE (CA->CH3, keep C, O)
                res.name = "ACE"
                for a in list(res.atoms):
                    if a.name == "CA":
                        a.name = "CH3"
                    elif a.name not in ("C", "O"):
                        del_atoms.append(a)
                report.append((res.chain, s, e))

    if not report:
        print("[truncate_and_cap] WARNING: no ranges matched any chain; writing an unchanged copy")

    # Put each post-cut segment on a fresh chain id so ParmEd writes a TER between the caps.
    for res in st.residues:
        if not eligible(res):
            continue
        for (s, e) in ranges:
            if res.number >= e:
                if res.chain not in seg2_chain:
                    if not free_chains:
                        break
                    seg2_chain[res.chain] = free_chains.pop(0)
                res.chain = seg2_chain[res.chain]
                break

    # Strip the removed atoms in one pass, by atom serial mask (1-based, pre-strip indices are stable here).
    if del_atoms:
        serials = sorted({a.idx + 1 for a in del_atoms})
        st.strip("@" + ",".join(map(str, serials)))

    # renumber=False keeps original residue numbers (so the numbering gap marks the removed loop and
    # partners outside it keep their numbers); the fresh chain id on the second segment gives the TER.
    st.save(args.out, overwrite=True, renumber=False)

    for ch, s, e in report:
        print(f"[truncate_and_cap] chain {ch}: removed {s+1}..{e-1}, cap NME@{s} + ACE@{e} "
              f"(kept flanks {s-1}, {e+1})")
    print(f"[truncate_and_cap] wrote {args.out}")


if __name__ == "__main__":
    main()
