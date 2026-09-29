#!/usr/bin/env python3
# Residue-level WT native contact retention (last 20 ns), 350 K, WT + 20 mutants x 2 replicas
import os, sys, glob, numpy as np, pandas as pd, MDAnalysis as mda
from multiprocessing import Pool
from MDAnalysis.coordinates.XTC import XTCReader
from MDAnalysis.lib.distances import distance_array, calc_bonds

B   = "/home/raza/p53_thermal_xai/data/structures"
TOP = "/home/raza/p53_thermal_xai/analysis/qc_native_contacts/top"
OUT = "/home/raza/p53_thermal_xai/analysis/qc_residue_contacts"
CUT, SEQSEP, WIN_PS = 4.5, 3, 20000.0

systems = sorted(os.path.basename(f)[:-7] for f in glob.glob(f"{TOP}/*_t0.pdb"))
if "WT" not in systems: sys.exit("WT topology PDB missing")

def heavy_by_res(u):
    ag = u.select_atoms("protein and not name H*")
    d = {}
    for a in ag: d.setdefault(int(a.resid), []).append(int(a.index))
    rn = {int(r.resid): r.resname for r in ag.residues}
    return {k: np.array(v) for k, v in d.items()}, rn

# ---- WT reference residue-pair network ----
uwt = mda.Universe(f"{TOP}/WT_t0.pdb")
wt_atoms, wt_rn = heavy_by_res(uwt)
resids = sorted(wt_atoms)
pos = uwt.atoms.positions
pairs = []
for a, i in enumerate(resids):
    for j in resids[a+1:]:
        if abs(j - i) <= SEQSEP: continue
        dmin = distance_array(pos[wt_atoms[i]], pos[wt_atoms[j]]).min()
        if dmin < CUT: pairs.append((i, j, round(float(dmin), 3)))
P = pd.DataFrame(pairs, columns=["res_i", "res_j", "WT_t0_dmin_A"])
P.insert(1, "WT_resname_i", [wt_rn[i] for i in P.res_i])
P.insert(3, "WT_resname_j", [wt_rn[j] for j in P.res_j])
print(f"# WT reference: {len(resids)} residues (resid {resids[0]}-{resids[-1]}), "
      f"{len(P)} residue-pair contacts (min heavy-atom distance < {CUT} Å, |i-j|>{SEQSEP})")

def run(args):
    s, r = args
    u = mda.Universe(f"{TOP}/{s}_t0.pdb")
    at, rn = heavy_by_res(u)
    if sorted(at) != resids: return (s, r, None, f"residue numbering differs from WT")
    A, Bb, starts, n = [], [], [], 0
    for i, j in zip(P.res_i, P.res_j):
        ai, aj = at[i], at[j]
        g = np.array(np.meshgrid(ai, aj, indexing="ij")).reshape(2, -1)
        starts.append(n); A.append(g[0]); Bb.append(g[1]); n += g.shape[1]
    A, Bb, starts = np.concatenate(A), np.concatenate(Bb), np.array(starts)
    xtc = f"{B}/c36_pbcfix/{s}_{r}.xtc"
    if not os.path.exists(xtc): return (s, r, None, "missing xtc")
    rd = XTCReader(xtc)
    if rd.n_atoms < len(u.atoms): return (s, r, None, f"atom mismatch xtc {rd.n_atoms} < pdb {len(u.atoms)}")
    start = rd.n_frames - int(round(WIN_PS / rd.dt)) - 1
    frac, occ, times = [], np.zeros(len(P)), []
    for ts in rd[start:]:
        d = calc_bonds(ts.positions[A], ts.positions[Bb], box=ts.dimensions)
        c = np.minimum.reduceat(d, starts) < CUT
        occ += c; frac.append(c.mean()); times.append(ts.time / 1000.0)
    rd.close()
    frac = np.array(frac)
    mut = {k: rn[k] for k in rn if rn[k] != wt_rn.get(k)}
    return (s, r, dict(frac=frac, occ=occ / len(frac), t0=times[0], t1=times[-1], mut=mut), None)

if __name__ == "__main__":
    jobs = [(s, r) for s in systems for r in ("r1", "r2")]
    with Pool(min(len(jobs), os.cpu_count() or 4)) as pool:
        res = pool.map(run, jobs)
    rows, occ_cols = [], {}
    for s, r, d, err in res:
        if err: print(f"ERROR {s} {r}: {err}"); continue
        occ_cols[f"{s}_{r}"] = np.round(d["occ"], 4)
        rows.append({"System": s, "Rep": r,
                     "Mutation_check": ",".join(f"{wt_rn[k]}{k}{v}" for k, v in sorted(d["mut"].items())) or "none",
                     "N_frames": len(d["frac"]), "Window_ns": f"{d['t0']:.2f}-{d['t1']:.2f}",
                     "WT_contacts_retained_last20ns": f"{d['frac'].mean():.3f} ± {d['frac'].std(ddof=1):.3f}"})
    order = ["WT"] + [x for x in systems if x != "WT"]
    df = pd.DataFrame(rows)
    df["System"] = pd.Categorical(df.System, order); df = df.sort_values(["System", "Rep"])
    df.to_csv(f"{OUT}/residue_contact_retention_last20ns_350K.tsv", sep="\t", index=False)
    occ = pd.concat([P, pd.DataFrame(occ_cols)[[f"{s}_{r}" for s in order for r in ("r1", "r2") if f"{s}_{r}" in occ_cols]]], axis=1)
    occ.to_csv(f"{OUT}/residue_pair_occupancy_last20ns_350K.tsv", sep="\t", index=False)
    print("\nWT residue-pair contacts retained, last 20 ns, 350 K (fraction per frame, mean ± SD)")
    print(df.to_string(index=False))
    print(f"\nSaved: {OUT}/residue_contact_retention_last20ns_350K.tsv\n"
          f"       {OUT}/residue_pair_occupancy_last20ns_350K.tsv  ({len(P)} pairs x {len(occ_cols)} runs)")
