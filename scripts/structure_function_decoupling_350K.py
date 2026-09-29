#!/usr/bin/env python3
import os, sys, glob, numpy as np, pandas as pd, MDAnalysis as mda
from multiprocessing import Pool
from MDAnalysis.coordinates.XTC import XTCReader
from MDAnalysis.lib.distances import distance_array, calc_bonds
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

B    = "/home/raza/p53_thermal_xai/data/structures"
TOP  = "/home/raza/p53_thermal_xai/analysis/qc_native_contacts/top"
FUNK = "/home/raza/p53_thermal_xai/data/external/funk2025/funk2025_scores_20mutants.tsv"
OUT  = "/home/raza/p53_thermal_xai/analysis/structure_function"
CUT, SEQSEP, T0, T1 = 4.5, 3, 20.0, 100.0
os.makedirs(OUT, exist_ok=True)
systems = sorted(os.path.basename(f)[:-7] for f in glob.glob(f"{TOP}/*_t0.pdb"))

def heavy_by_res(u):
    d = {}
    for a in u.select_atoms("protein and not name H*"): d.setdefault(int(a.resid), []).append(int(a.index))
    return {k: np.array(v) for k, v in d.items()}

uwt = mda.Universe(f"{TOP}/WT_t0.pdb"); wt_atoms = heavy_by_res(uwt); resids = sorted(wt_atoms); pos = uwt.atoms.positions
PAIRS = [(i, j) for a, i in enumerate(resids) for j in resids[a+1:]
         if abs(j-i) > SEQSEP and distance_array(pos[wt_atoms[i]], pos[wt_atoms[j]]).min() < CUT]
print(f"# WT residue-pair contacts: {len(PAIRS)} (min heavy-atom < {CUT} Å, |i-j|>{SEQSEP}); window {T0}-{T1} ns")

def run(args):
    s, r = args
    u = mda.Universe(f"{TOP}/{s}_t0.pdb"); at = heavy_by_res(u)
    if sorted(at) != resids: return (s, r, None, "residue numbering differs")
    A, Bb, st, n = [], [], [], 0
    for i, j in PAIRS:
        g = np.array(np.meshgrid(at[i], at[j], indexing="ij")).reshape(2, -1)
        st.append(n); A.append(g[0]); Bb.append(g[1]); n += g.shape[1]
    A, Bb, st = np.concatenate(A), np.concatenate(Bb), np.array(st)
    xtc = f"{B}/c36_pbcfix/{s}_{r}.xtc"
    if not os.path.exists(xtc): return (s, r, None, "missing xtc")
    rd = XTCReader(xtc); t, f = [], []
    for ts in rd:
        d = calc_bonds(ts.positions[A], ts.positions[Bb], box=ts.dimensions)
        f.append((np.minimum.reduceat(d, st) < CUT).mean()); t.append(ts.time/1000.0)
    rd.close()
    return (s, r, (np.array(t), np.array(f)), None)

if __name__ == "__main__":
    jobs = [(s, r) for s in systems for r in ("r1", "r2")]
    with Pool(min(len(jobs), os.cpu_count() or 4)) as p: res = p.map(run, jobs)
    ts = {}
    for s, r, d, err in res:
        if err: print(f"ERROR {s} {r}: {err}"); continue
        ts[(s, r)] = d
    np.savez_compressed(f"{OUT}/residue_contact_retention_timeseries_350K.npz",
                        **{f"{s}_{r}_t": v[0] for (s, r), v in ts.items()}, **{f"{s}_{r}_f": v[1] for (s, r), v in ts.items()})
    def win(s, r):
        t, f = ts[(s, r)]; m = (t >= T0) & (t <= T1); return f[m].mean(), m.sum()
    wt = np.mean([win("WT", "r1")[0], win("WT", "r2")[0]])
    print(f"# WT: r1={win('WT','r1')[0]:.4f} r2={win('WT','r2')[0]:.4f} mean={wt:.4f}")
    funk = pd.read_csv(FUNK, sep="\t").set_index("Mutation")
    rows = []
    for s in [x for x in systems if x != "WT"]:
        q1, n1 = win(s, "r1"); q2, n2 = win(s, "r2"); qm = (q1 + q2) / 2
        rows.append({"Mutation": s, "Q_r1_20-100": q1, "Q_r2_20-100": q2, "Q_mean": qm,
                     "Delta_Q_vs_WT": qm - wt, "dQ_r1": q1 - wt, "dQ_r2": q2 - wt,
                     "CRISPR_RFS": funk.loc[s, "rfs_crispr"] if s in funk.index else np.nan, "n_frames": f"{n1}/{n2}"})
    df = pd.DataFrame(rows)
    df.round(4).to_csv(f"{OUT}/structure_function_table_350K.tsv", sep="\t", index=False)
    print(df[["Mutation", "Q_r1_20-100", "Q_r2_20-100", "Q_mean", "Delta_Q_vs_WT", "CRISPR_RFS", "n_frames"]].round(4).to_string(index=False))
    ok = df.dropna(subset=["CRISPR_RFS"])
    rho, p = spearmanr(ok.Delta_Q_vs_WT, ok.CRISPR_RFS)
    print(f"\nSpearman rho(Delta_Q_vs_WT, CRISPR_RFS) = {rho:.3f}, p = {p:.3g}, n = {len(ok)}")
    fig, ax = plt.subplots(figsize=(7, 5))
    lo = ok[["dQ_r1", "dQ_r2"]].min(axis=1); hi = ok[["dQ_r1", "dQ_r2"]].max(axis=1)
    ax.errorbar(ok.Delta_Q_vs_WT, ok.CRISPR_RFS, xerr=[ok.Delta_Q_vs_WT - lo, hi - ok.Delta_Q_vs_WT],
                fmt="o", color="tab:blue", ecolor="grey", elinewidth=1, capsize=3, ms=5)
    for _, r in ok.iterrows(): ax.annotate(r.Mutation, (r.Delta_Q_vs_WT, r.CRISPR_RFS), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.axvline(0, ls="--", lw=0.7, color="k"); ax.axhline(0, ls=":", lw=0.7, color="k")
    ax.set_xlabel("ΔQ vs WT (WT residue-pair contact retention, 20–100 ns, 350 K)\nbars = r1–r2 range")
    ax.set_ylabel("CRISPR RFS (Funk et al. 2025; +1 nonsense-like, −1 synonymous-like)")
    ax.set_title(f"Structure–function: Spearman ρ = {rho:.2f} (p = {p:.2g}, n = {len(ok)})")
    fig.tight_layout(); fig.savefig(f"{OUT}/deltaQ_vs_CRISPR_RFS_350K.png", dpi=300)
    print(f"\nSaved: {OUT}/structure_function_table_350K.tsv\n       {OUT}/deltaQ_vs_CRISPR_RFS_350K.png\n       {OUT}/residue_contact_retention_timeseries_350K.npz")
