#!/usr/bin/env python3
"""Build p53 DBD hotspot mutants (Table 6 panel) from the WT structure."""
import os
import sys
from pymol import cmd

args = sys.argv[1:]
if "--" in sys.argv:
    args = sys.argv[sys.argv.index("--") + 1:]
if len(args) < 2:
    sys.exit("Usage: pymol -cq build_mutants.py -- <wt.pdb> <outdir> [chain]")

WT_PDB = args[0]
OUTDIR = args[1]
CHAIN  = args[2] if len(args) > 2 else "A"
os.makedirs(OUTDIR, exist_ok=True)

MUTATIONS = [
    (175, "ALA", "R175A"),
    (134, "LEU", "F134L"),
    (249, "SER", "R249S"),
    (237, "ILE", "M237I"),
    (242, "SER", "C242S"),
    (143, "ALA", "V143A"),
    (145, "GLN", "L145Q"),
    (151, "SER", "P151S"),
    (157, "PHE", "V157F"),
    (195, "THR", "I195T"),
    (232, "THR", "I232T"),
    (255, "PHE", "I255F"),
    (270, "CYS", "F270C"),
]


def build_one(resi, target3, label):
    cmd.reinitialize()
    cmd.load(WT_PDB, "m")
    sel = f"chain {CHAIN} and resi {resi}"
    if cmd.count_atoms(sel) == 0:
        print(f"[SKIP] {label}: no atoms at 'chain {CHAIN} and resi {resi}'")
        return
    cmd.wizard("mutagenesis")
    cmd.refresh_wizard()
    cmd.get_wizard().set_mode(target3)
    cmd.get_wizard().do_select(sel)
    cmd.get_wizard().apply()
    cmd.set_wizard()
    out = os.path.join(OUTDIR, f"{label}.pdb")
    cmd.save(out, "polymer or resn ZN")
    print(f"[ok] {label:<6} at {resi:>3} -> {out}")

for resi, target3, label in MUTATIONS:
    build_one(resi, target3, label)
print("\nDone.")
