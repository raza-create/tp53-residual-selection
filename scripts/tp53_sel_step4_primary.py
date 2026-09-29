#!/usr/bin/env python3
# Step 4 PRIMARY TEST: delta-RFS per cancer + cross-cancer functional convergence.
# Null = sample-matched mutational opportunity, conditional on a TP53 126-292 missense SNV occurring.
import io, zipfile, numpy as np, pandas as pd

P  = "/home/raza/p53_thermal_xai/data/external/pcawg_signatures"
D  = "/home/raza/p53_thermal_xai/data/external/tp53_selection"
NB, SEED = 100000, 20260929
ORDER = ["LUSC","OV","HNSC","BRCA","COADREAD","LUAD","STAD","BLCA","UCEC","GBM","LIHC","PAAD"]
rng = np.random.default_rng(SEED)

uni = pd.read_csv(D + "/tp53_universe_126_292.tsv", sep="\t")
uni = uni[uni.CRISPR_RFS.notna()].reset_index(drop=True)
obs = pd.read_csv(D + "/tcga_observed_events_126_292.tsv", sep="\t")
print("universe = %d candidates | observed events = %d | cancers = %d"
      % (len(uni), len(obs), obs.cancer.nunique()))

chan = pd.read_csv(P + "/COSMIC_v3_SBS_GRCh37_WESrenorm.txt", sep="\t", index_col=0)
sigs = list(chan.columns)
Eexp = pd.read_csv(P + "/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
Eexp["patient"] = Eexp["Sample Names"].str.slice(0, 12)
Eexp = Eexp.drop_duplicates("patient").set_index("patient")

G = pd.read_csv(P + "/context_counts_GRCh37_96.csv", index_col=0).sum(1)   # unused, kept for record
O = pd.read_csv(P + "/context_counts_GRCh37_96_exome.csv", index_col=0).sum(1)
tri = chan.index.str[0] + chan.index.str[2] + chan.index.str[6]
Ovec = O.reindex(tri.values).values.astype(float)

ci = np.array([list(chan.index).index(c) for c in uni.sbs96])   # candidate -> channel row
rfs = uni.CRISPR_RFS.values.astype(float)

pats = obs.patient.unique()
miss = [p for p in pats if p not in Eexp.index]
if miss: print("patients without exposures (dropped): %d" % len(miss))
obs = obs[obs.patient.isin(Eexp.index)].reset_index(drop=True)

Em = Eexp.loc[pats[np.isin(pats, Eexp.index)], sigs].astype(float).values
pidx = {p: k for k, p in enumerate(Eexp.loc[pats[np.isin(pats, Eexp.index)]].index)}
q = Em @ chan.values.T                    # patients x 96 expected WES spectrum
r = q / Ovec[None, :]                     # per-site rate
W = r[:, ci]                              # patients x candidates
Pm = W / W.sum(1, keepdims=True)
mu_pat = Pm @ rfs
neff = 1.0 / (Pm**2).sum(1)

row = np.array([pidx[p] for p in obs.patient])
obs["mu"] = mu_pat[row]
obs["neff"] = neff[row]
cum = np.cumsum(Pm, axis=1)[row]          # events x candidates, per-event CDF

res = []
for ct in ORDER:
    m = obs.cancer.values == ct
    res.append(dict(cancer=ct, n=int(m.sum()),
                    obs=obs.CRISPR_RFS.values[m].mean(), exact=obs.mu.values[m].mean()))
R = pd.DataFrame(res); R["delta"] = R.obs - R.exact

idx = {ct: np.where(obs.cancer.values == ct)[0] for ct in ORDER}
null = np.empty((NB, len(ORDER))); V = np.empty(NB)
for b in range(NB):
    u = rng.random(len(obs))
    pick = (cum < u[:, None]).sum(1)
    v = rfs[np.minimum(pick, len(rfs)-1)]
    means = np.array([v[idx[ct]].mean() for ct in ORDER])
    null[b] = means; V[b] = means.var(ddof=1)

R["null_lo"] = np.percentile(null, 2.5, axis=0)
R["null_hi"] = np.percentile(null, 97.5, axis=0)
R["P"] = [(1 + (null[:, k] >= R.obs.values[k]).sum()) / (NB + 1) for k in range(len(ORDER))]

print("\nCancer      n  Observed_mean_RFS  Exact_null_mean_RFS  Delta_RFS  Null_2.5%  Null_97.5%  empirical_P")
for r_ in R.itertuples():
    print("%-9s %4d %18.4f %20.4f %10.4f %10.4f %11.4f %12.5f"
          % (r_.cancer, r_.n, r_.obs, r_.exact, r_.delta, r_.null_lo, r_.null_hi, r_.P))

Vobs = R.obs.values.var(ddof=1)
print("\nGLOBAL CONVERGENCE")
print("Observed between-cancer variance = %.6g" % Vobs)
print("Median null variance             = %.6g" % np.median(V))
print("95%% null interval                = %.6g - %.6g" % (np.percentile(V,2.5), np.percentile(V,97.5)))
print("Observed/null variance ratio     = %.4f" % (Vobs/np.median(V)))
print("Empirical convergence P          = %.5f" % ((1 + (V <= Vobs).sum())/(NB+1)))

print("\nQC")
print("Maximum |Monte-Carlo null mean - exact null mean| = %.6f"
      % np.abs(null.mean(0) - R.exact.values).max())
print("Minimum effective number of candidate mutations across samples = %.1f" % obs.neff.min())
print("  (median N_eff = %.1f, 5th pct = %.1f)" % (obs.neff.median(), obs.neff.quantile(.05)))
R.round(6).to_csv(D + "/primary_deltaRFS_results.tsv", sep="\t", index=False)
np.savez_compressed(D + "/null_distributions.npz", null=null, V=V, cancers=np.array(ORDER))
print("\nSaved: %s/primary_deltaRFS_results.tsv and null_distributions.npz" % D)
