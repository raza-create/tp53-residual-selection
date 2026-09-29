#!/usr/bin/env python3
# Reproducibility of MD-derived structural metrics across independent replicas.
# 21 genotypes x 2 seeds. Trajectory summaries (20-100 ns) are the observations, never frames.
import os, sys, glob, numpy as np, pandas as pd, MDAnalysis as mda
from multiprocessing import Pool
from MDAnalysis.coordinates.XTC import XTCReader
from MDAnalysis.lib.distances import distance_array, calc_bonds
from scipy.stats import f as fdist

B   = "/home/raza/p53_thermal_xai/data/structures"
TOP = "/home/raza/p53_thermal_xai/analysis/qc_native_contacts/top"
RMS = "/home/raza/p53_thermal_xai/analysis/qc_rmsd/xvg"
OUT = "/home/raza/p53_thermal_xai/analysis/replica_icc"; os.makedirs(OUT, exist_ok=True)
CUT, SEQSEP, T0, T1 = 4.5, 3, 20.0, 100.0
LOOPS = {"L1": list(range(112,125)), "L2_Zn": list(range(164,195))+[238,242],
         "L3": list(range(237,251)), "H2": list(range(278,287))}
REGIONS = ["beta_core","L1","L2_Zn","L3","H2","other_DBD"]

systems = sorted(os.path.basename(f)[:-7] for f in glob.glob(TOP+"/*_t0.pdb"))
if "WT" not in systems: sys.exit("WT topology missing")

def heavy_by_res(u):
    d = {}
    for a in u.select_atoms("protein and not name H*"): d.setdefault(int(a.resid), []).append(int(a.index))
    return {k: np.array(v) for k, v in d.items()}

uwt = mda.Universe(TOP+"/WT_t0.pdb"); wt_at = heavy_by_res(uwt)
resids = sorted(wt_at); pos = uwt.atoms.positions
try:
    from MDAnalysis.analysis.dssp import DSSP
    ss = DSSP(uwt).run().results.dssp[0]
    strand = {r for r, c in zip([int(x.resid) for x in uwt.select_atoms("protein").residues], ss) if c == "E"}
except Exception as e:
    sys.exit("DSSP failed: %s" % e)

def region_of(r):
    if r in LOOPS["L2_Zn"]: return "L2_Zn"
    for k in ("L1","L3","H2"):
        if r in LOOPS[k]: return k
    return "beta_core" if r in strand else "other_DBD"
reg = {r: region_of(r) for r in resids}

PAIRS = [(i,j) for a,i in enumerate(resids) for j in resids[a+1:]
         if abs(j-i) > SEQSEP and distance_array(pos[wt_at[i]], pos[wt_at[j]]).min() < CUT]
pair_reg = [{reg[i], reg[j]} for i, j in PAIRS]
masks = {g: np.array([g in s for s in pair_reg]) for g in REGIONS}
print("# WT native contacts: %d | region pair counts: %s" % (len(PAIRS), {g:int(m.sum()) for g,m in masks.items()}))
print("# window %.0f-%.0f ns; observations = trajectory summaries (n=2 per genotype)" % (T0, T1))

def run(args):
    s, r = args
    u = mda.Universe(TOP+"/%s_t0.pdb" % s); at = heavy_by_res(u)
    if sorted(at) != resids: return (s, r, None, "residue numbering differs")
    A, Bb, st, n = [], [], [], 0
    for i, j in PAIRS:
        g = np.array(np.meshgrid(at[i], at[j], indexing="ij")).reshape(2, -1)
        st.append(n); A.append(g[0]); Bb.append(g[1]); n += g.shape[1]
    A, Bb, st = np.concatenate(A), np.concatenate(Bb), np.array(st)
    xtc = B+"/c36_pbcfix/%s_%s.xtc" % (s, r)
    if not os.path.exists(xtc): return (s, r, None, "missing xtc")
    rd = XTCReader(xtc); acc, t = [], []
    for ts in rd:
        d = calc_bonds(ts.positions[A], ts.positions[Bb], box=ts.dimensions)
        acc.append(np.minimum.reduceat(d, st) < CUT); t.append(ts.time/1000.0)
    rd.close()
    acc = np.array(acc); t = np.array(t); w = (t >= T0) & (t <= T1)
    res = {"Q": acc[w].mean()}
    for g, m in masks.items(): res[g] = acc[w][:, m].mean()
    rf = RMS+"/%s_%s_rmsd.xvg" % (s, r)
    if os.path.exists(rf):
        a = np.loadtxt(rf); mm = (a[:,0] >= T0) & (a[:,0] <= T1); res["RMSD"] = a[mm,1].mean()*10.0
    else:
        res["RMSD"] = np.nan
    res["n_frames"] = int(w.sum())
    return (s, r, res, None)

def icc1(y, k=2):
    """One-way random effects ICC(1,1) on groups x k; y is (n_groups, k)."""
    n = y.shape[0]; gm = y.mean()
    msb = k*((y.mean(1)-gm)**2).sum()/(n-1)
    msw = ((y - y.mean(1, keepdims=True))**2).sum()/(n*(k-1))
    icc = (msb-msw)/(msb+(k-1)*msw)
    F = msb/msw; df1, df2 = n-1, n*(k-1)
    fu = F/fdist.ppf(0.975, df1, df2); fl = F/fdist.ppf(0.025, df1, df2)
    lo = (fu-1)/(fu+k-1); hi = (fl-1)/(fl+k-1)
    s_b = max((msb-msw)/k, 0.0)
    return icc, lo, hi, s_b, msw, F, fdist.sf(F, df1, df2)

if __name__ == "__main__":
    jobs = [(s, r) for s in systems for r in ("r1","r2")]
    with Pool(min(len(jobs), os.cpu_count() or 4)) as p: res = p.map(run, jobs)
    rows = []
    for s, r, d, err in res:
        if err: print("ERROR %s %s: %s" % (s, r, err)); continue
        rows.append(dict(Mutation=s, Replica=r, RMSD=d["RMSD"], Q=d["Q"],
                         beta_core=d["beta_core"], L1=d["L1"], L2_Zn=d["L2_Zn"],
                         L3=d["L3"], H2=d["H2"], other_DBD=d["other_DBD"], n_frames=d["n_frames"]))
    df = pd.DataFrame(rows)
    order = ["WT"]+[x for x in systems if x != "WT"]
    df["Mutation"] = pd.Categorical(df.Mutation, order); df = df.sort_values(["Mutation","Replica"])
    df.round(4).to_csv(OUT+"/trajectory_summaries_20_100ns.tsv", sep="\t", index=False)
    print(df[["Mutation","Replica","RMSD","Q","beta_core","L1","L2_Zn","L3","H2","other_DBD"]].round(4).to_string(index=False))

    metrics = ["RMSD","Q","beta_core","L1","L2_Zn","L3","H2","other_DBD"]
    w = df.pivot(index="Mutation", columns="Replica", values=metrics)
    out = []
    for m in metrics:
        y = w[m][["r1","r2"]].dropna().values
        icc, lo, hi, s_b, s_w, F, p = icc1(y)
        out.append(dict(Metric=m, n_genotypes=len(y), ICC=icc, CI_low=lo, CI_high=hi,
                        var_between=s_b, var_within=s_w, sd_between=np.sqrt(s_b), sd_within=np.sqrt(s_w),
                        F=F, p=p, mean=y.mean(), range_of_genotype_means=y.mean(1).ptp()))
    t = pd.DataFrame(out)
    t.to_csv(OUT+"/icc_results.tsv", sep="\t", index=False)
    print("\nICC(1,1), one-way random effects, %d genotypes x 2 replicas" % len(w))
    print(t[["Metric","ICC","CI_low","CI_high","var_between","var_within","sd_between","sd_within","F","p"]]
          .round(4).to_string(index=False))
    print("\nScale context (genotype means):")
    print(t[["Metric","mean","range_of_genotype_means","sd_between","sd_within"]].round(4).to_string(index=False))
    print("\nSaved: "+OUT+"/trajectory_summaries_20_100ns.tsv")
    print("       "+OUT+"/icc_results.tsv")
