#!/usr/bin/env python3
import re, sys, numpy as np, pandas as pd, MDAnalysis as mda
OCC = "/home/raza/p53_thermal_xai/analysis/qc_residue_contacts/residue_pair_occupancy_last20ns_350K.tsv"
TOP = "/home/raza/p53_thermal_xai/analysis/qc_native_contacts/top/WT_t0.pdb"
OUT = "/home/raza/p53_thermal_xai/analysis/regional_contacts"
MIN_THR = 0.02
LOOPS = {"L1": range(112, 125), "L2/Zn": list(range(164, 195)) + [238, 242], "L3": range(237, 251), "H2": range(278, 287)}
import os; os.makedirs(OUT, exist_ok=True)

u = mda.Universe(TOP)
try:
    from MDAnalysis.analysis.dssp import DSSP
    ss = DSSP(u).run().results.dssp[0]
    resids = [int(r.resid) for r in u.select_atoms("protein").residues]
    strand = {r for r, s in zip(resids, ss) if s == "E"}
except Exception as e:
    sys.exit(f"DSSP failed ({e}); install MDAnalysis>=2.8")

def region_of(r):
    if r in LOOPS["L2/Zn"]: return "L2/Zn"
    for k in ("L1", "L3", "H2"):
        if r in LOOPS[k]: return k
    return "beta-core" if r in strand else "other-DBD"

allres = sorted(set(resids))
reg = {r: region_of(r) for r in allres}
print("# Region definitions (WT DSSP for beta-core):")
for k in ["beta-core", "L1", "L2/Zn", "L3", "H2", "other-DBD"]:
    rr = [r for r in allres if reg[r] == k]; print(f"#   {k:10s} n={len(rr):3d}  {rr}")

occ = pd.read_csv(OCC, sep="\t")
occ["reg_i"] = occ.res_i.map(reg); occ["reg_j"] = occ.res_j.map(reg)
runs = [c for c in occ.columns if re.match(r".+_r[123]$", c)]
systems = sorted({c.rsplit("_", 1)[0] for c in runs} - {"WT"})
regions = ["beta-core", "L1", "L2/Zn", "L3", "H2", "other-DBD"]

def mask_region(k): return (occ.reg_i == k) | (occ.reg_j == k)

rows = []
for s in systems:
    mres = int(re.findall(r"\d+", s)[0])
    reps = sorted(c for c in runs if c.rsplit("_", 1)[0] == s)
    regs = regions + ["MutSite"]
    for k in regs:
        m = ((occ.res_i == mres) | (occ.res_j == mres)) if k == "MutSite" else mask_region(k)
        n = int(m.sum())
        if n == 0: continue
        wt1, wt2 = occ.loc[m, "WT_r1"].mean(), occ.loc[m, "WT_r2"].mean()
        wt = (wt1 + wt2) / 2; thr = max(MIN_THR, abs(wt1 - wt2))
        row = {"Mutation": s, "Region": k, "n_pairs": n, "WT_r1": wt1, "WT_r2": wt2,
               "WT_contact_occupancy": wt, "noise_thr": thr}
        deltas = []
        for c in reps:
            v = occ.loc[m, c].mean(); row["Mut_" + c.rsplit("_", 1)[1]] = v; deltas.append(v - wt)
            row["d_" + c.rsplit("_", 1)[1]] = v - wt
        if all(d < -thr for d in deltas): call = "DECREASE"
        elif all(d > thr for d in deltas): call = "INCREASE"
        elif all(abs(d) <= thr for d in deltas): call = "no_change"
        else: call = "inconsistent"
        row["Reproducible"] = call; row["n_reps"] = len(reps)
        rows.append(row)

df = pd.DataFrame(rows)
num = [c for c in df.columns if c.startswith(("WT", "Mut_", "d_", "noise"))]
df[num] = df[num].round(3)
df.to_csv(f"{OUT}/regional_contact_occupancy_350K.tsv", sep="\t", index=False)
cols = ["Mutation", "Region", "n_pairs", "WT_contact_occupancy", "noise_thr"] + \
       sorted(c for c in df.columns if c.startswith("Mut_")) + sorted(c for c in df.columns if c.startswith("d_")) + ["Reproducible"]
print(df[cols].to_string(index=False))
sym = {"DECREASE": "DOWN", "INCREASE": "UP", "no_change": ".", "inconsistent": "~"}
piv = df.pivot(index="Mutation", columns="Region", values="Reproducible").reindex(columns=regions + ["MutSite"]).replace(sym).fillna("n/a")
piv.to_csv(f"{OUT}/regional_reproducibility_matrix_350K.tsv", sep="\t")
print("\nDOWN/UP = both replicas beyond noise_thr (max(0.02, |WT_r1-WT_r2|)) in same direction; ~ = replicas disagree; . = within noise")
print(piv.to_string())
print(f"\nSaved: {OUT}/regional_contact_occupancy_350K.tsv\n       {OUT}/regional_reproducibility_matrix_350K.tsv")
