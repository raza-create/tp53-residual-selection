#!/usr/bin/env python3
# Step 5: distribution-level comparison. Pairwise 1-D Wasserstein between cancer types,
# for mutation-supply (expected) vs observed tumour RFS distributions. No p-values.
import itertools, numpy as np, pandas as pd
from scipy.stats import wasserstein_distance

P = "/home/raza/p53_thermal_xai/data/external/pcawg_signatures"
D = "/home/raza/p53_thermal_xai/data/external/tp53_selection"
ORDER = ["LUSC","OV","HNSC","BRCA","COADREAD","LUAD","STAD","BLCA","UCEC","GBM","LIHC","PAAD"]

uni = pd.read_csv(D + "/tp53_universe_126_292.tsv", sep="\t")
uni = uni[uni.CRISPR_RFS.notna()].reset_index(drop=True)
obs = pd.read_csv(D + "/tcga_observed_events_126_292.tsv", sep="\t")
rfs = uni.CRISPR_RFS.values.astype(float)

chan = pd.read_csv(P + "/COSMIC_v3_SBS_GRCh37_WESrenorm.txt", sep="\t", index_col=0)
sigs = list(chan.columns)
Ex = pd.read_csv(P + "/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
Ex["patient"] = Ex["Sample Names"].str.slice(0, 12)
Ex = Ex.drop_duplicates("patient").set_index("patient")
O = pd.read_csv(P + "/context_counts_GRCh37_96_exome.csv", index_col=0).sum(1)
tri = chan.index.str[0] + chan.index.str[2] + chan.index.str[6]
Ovec = O.reindex(tri.values).values.astype(float)
ci = np.array([list(chan.index).index(c) for c in uni.sbs96])

obs = obs[obs.patient.isin(Ex.index)].reset_index(drop=True)
pl = obs.patient.unique()
q = Ex.loc[pl, sigs].astype(float).values @ chan.values.T
W = (q / Ovec[None, :])[:, ci]
Pm = W / W.sum(1, keepdims=True)
pidx = {p: k for k, p in enumerate(pl)}

Wt, Ot = {}, {}
for ct in ORDER:
    m = obs.cancer.values == ct
    rows = np.array([pidx[p] for p in obs.patient.values[m]])
    Wt[ct] = Pm[rows].mean(0)
    Ot[ct] = obs.CRISPR_RFS.values[m]
    assert abs(Wt[ct].sum() - 1) < 1e-9

pairs = list(itertools.combinations(ORDER, 2))
rows = []
for a, b in pairs:
    we = wasserstein_distance(rfs, rfs, Wt[a], Wt[b])
    wo = wasserstein_distance(Ot[a], Ot[b])
    rows.append(dict(A=a, B=b, W_expected=we, W_observed=wo))
T = pd.DataFrame(rows)
T.round(5).to_csv(D + "/wasserstein_pairs.tsv", sep="\t", index=False)

def blk(x): return x.mean(), x.median(), x.min(), x.max()
me = blk(T.W_expected); mo = blk(T.W_observed)
print("EXPECTED mutation-supply distributions (%d pairs)" % len(T))
print("Mean pairwise Wasserstein distance = %.5f" % me[0])
print("Median                             = %.5f" % me[1])
print("Min                                = %.5f" % me[2])
print("Max                                = %.5f" % me[3])
print("\nOBSERVED tumour distributions")
print("Mean pairwise Wasserstein distance = %.5f" % mo[0])
print("Median                             = %.5f" % mo[1])
print("Min                                = %.5f" % mo[2])
print("Max                                = %.5f" % mo[3])
print("\nObserved / expected mean-distance ratio = %.4f" % (mo[0]/me[0]))

print("\nTop 10 pairs by expected distance")
print("Cancer_A  Cancer_B  W_expected  W_observed")
for r in T.nlargest(10, "W_expected").itertuples():
    print("%-9s %-9s %11.5f %11.5f" % (r.A, r.B, r.W_expected, r.W_observed))

print("\nContext: spread of each distribution's own RFS")
print("%-9s %10s %10s %10s %10s" % ("Cancer", "exp_mean", "exp_SD", "obs_mean", "obs_SD"))
for ct in ORDER:
    em = (Wt[ct]*rfs).sum(); esd = np.sqrt((Wt[ct]*(rfs-em)**2).sum())
    print("%-9s %10.4f %10.4f %10.4f %10.4f" % (ct, em, esd, Ot[ct].mean(), Ot[ct].std(ddof=1)))
print("\nuniverse RFS: mean %.4f SD %.4f range %.4f..%.4f" % (rfs.mean(), rfs.std(ddof=1), rfs.min(), rfs.max()))
print("Saved: %s/wasserstein_pairs.tsv" % D)
