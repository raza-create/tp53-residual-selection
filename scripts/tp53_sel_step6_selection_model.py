#!/usr/bin/env python3
# Step 6: exponential-tilt selection model on top of sample-specific mutation supply.
# P_im(beta) propto p_im * exp(beta * RFS_m). MLE per cancer + global; LRT; distribution fit.
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, wasserstein_distance, norm

P = "/home/raza/p53_thermal_xai/data/external/pcawg_signatures"
D = "/home/raza/p53_thermal_xai/data/external/tp53_selection"
ORDER = ["LUSC","OV","HNSC","BRCA","COADREAD","LUAD","STAD","BLCA","UCEC","GBM","LIHC","PAAD"]

uni = pd.read_csv(D + "/tp53_universe_126_292.tsv", sep="\t")
uni = uni[uni.CRISPR_RFS.notna()].reset_index(drop=True)
obs = pd.read_csv(D + "/tcga_observed_events_126_292.tsv", sep="\t")

chan = pd.read_csv(P + "/COSMIC_v3_SBS_GRCh37_WESrenorm.txt", sep="\t", index_col=0)
sigs = list(chan.columns)
Ex = pd.read_csv(P + "/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
Ex["patient"] = Ex["Sample Names"].str.slice(0, 12); Ex = Ex.drop_duplicates("patient").set_index("patient")
O = pd.read_csv(P + "/context_counts_GRCh37_96_exome.csv", index_col=0).sum(1)
tri = chan.index.str[0] + chan.index.str[2] + chan.index.str[6]
Ovec = O.reindex(tri.values).values.astype(float)
ci = np.array([list(chan.index).index(c) for c in uni.sbs96])

obs = obs[obs.patient.isin(Ex.index)].reset_index(drop=True)
pl = obs.patient.unique()
q = Ex.loc[pl, sigs].astype(float).values @ chan.values.T
Wq = (q / Ovec[None, :])[:, ci]
Pm = Wq / Wq.sum(1, keepdims=True)
pidx = {p: k for k, p in enumerate(pl)}

# collapse nucleotide routes to distinct protein substitutions (RFS is per protein change)
key = uni.wt_aa + uni.residue.astype(str) + uni.alt_aa
codes, uq = pd.factorize(key)
G = np.zeros((Pm.shape[0], len(uq)))
np.add.at(G.T, codes, Pm.T)
rfsG = pd.Series(uni.CRISPR_RFS.values).groupby(codes).first().values
print("candidates %d -> distinct protein substitutions %d" % (len(uni), len(uq)))

okey = obs.wt_aa + obs.residue.astype(str) + obs.alt_aa
gmap = {k: i for i, k in enumerate(uq)}
missing = sorted(set(okey) - set(gmap))
if missing: print("observed substitutions not in universe (dropped): %s" % missing)
keep = okey.isin(gmap).values
obs, okey = obs[keep].reset_index(drop=True), okey[keep].reset_index(drop=True)
gi = okey.map(gmap).values
pi = np.array([pidx[p] for p in obs.patient])
print("events used = %d" % len(obs))

def nll(b, rows, gs):
    Z = G[rows] * np.exp(b * rfsG)[None, :]
    return -(np.log(Z[np.arange(len(gs)), gs]) - np.log(Z.sum(1))).sum()
def fit(rows, gs):
    r = minimize_scalar(nll, bounds=(-20, 20), args=(rows, gs), method="bounded",
                        options={"xatol": 1e-8})
    b = r.x
    Z = G[rows] * np.exp(b * rfsG)[None, :]; Pb = Z / Z.sum(1, keepdims=True)
    m1 = Pb @ rfsG; m2 = Pb @ (rfsG**2)
    info = (m2 - m1**2).sum()
    return b, -r.fun, 1/np.sqrt(info)

print("\nCancer      n     beta_t       SE            95% CI")
rows_out = []
ll_t = 0.0
for ct in ORDER:
    m = obs.cancer.values == ct
    b, ll, se = fit(pi[m], gi[m]); ll_t += ll
    rows_out.append(dict(cancer=ct, n=int(m.sum()), beta=b, se=se,
                         lo=b-1.96*se, hi=b+1.96*se, ll=ll))
    print("%-9s %4d %10.4f %8.4f   %8.4f - %.4f" % (ct, m.sum(), b, se, b-1.96*se, b+1.96*se))
bg, ll_g, se_g = fit(pi, gi)
print("\nGlobal beta        = %.4f" % bg)
print("Global beta 95%% CI = %.4f - %.4f" % (bg-1.96*se_g, bg+1.96*se_g))
k = len(ORDER)
lr = 2*(ll_t - ll_g); df = k-1
print("\nLogLik global model          = %.3f" % ll_g)
print("LogLik tissue-specific model = %.3f" % ll_t)
print("AIC global                   = %.3f" % (-2*ll_g + 2*1))
print("AIC tissue-specific          = %.3f" % (-2*ll_t + 2*k))
print("Likelihood-ratio statistic   = %.3f" % lr)
print("df                           = %d" % df)
print("p                            = %.4g" % chi2.sf(lr, df))

print("\nCancer     W_observed_vs_supply   W_observed_vs_fitted_selection")
for r in rows_out:
    ct = r["cancer"]; m = obs.cancer.values == ct
    sup = G[pi[m]].mean(0)
    Z = G[pi[m]] * np.exp(r["beta"] * rfsG)[None, :]; fit_d = (Z/Z.sum(1, keepdims=True)).mean(0)
    ov = obs.CRISPR_RFS.values[m]
    w1 = wasserstein_distance(rfsG, ov, sup, np.ones(len(ov)))
    w2 = wasserstein_distance(rfsG, ov, fit_d, np.ones(len(ov)))
    r["W_supply"], r["W_fitted"] = w1, w2
    print("%-9s %19.5f %31.5f" % (ct, w1, w2))
pd.DataFrame(rows_out).round(5).to_csv(D + "/selection_model_beta.tsv", sep="\t", index=False)
print("\nSaved: %s/selection_model_beta.tsv" % D)
