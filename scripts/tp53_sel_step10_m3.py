#!/usr/bin/env python3
# Step 10: M3 adds the two p53-null Giacomelli screens. Full fit, cross-fit, Poisson-binomial.
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import chi2

P = "/home/raza/p53_thermal_xai/data/external/pcawg_signatures"
D = "/home/raza/p53_thermal_xai/data/external/tp53_selection"
GIA = "/home/raza/p53_thermal_xai/data/external/giacomelli2018/41588_2018_204_MOESM6_ESM.xlsx"
SEED, NFOLD = 20260929, 5
ORDER = ["LUSC","OV","HNSC","BRCA","COADREAD","LUAD","STAD","BLCA","UCEC","GBM","LIHC","PAAD"]
FOCUS = ["Y220C","V157F","I195T","R248Q","R273H","R273L","R273S","G266V",
         "H179R","R158L","R175H","R248W","R249S","R282W","G245S"]
rng = np.random.default_rng(SEED)

uni = pd.read_csv(D + "/tp53_universe_126_292.tsv", sep="\t")
uni = uni[uni.CRISPR_RFS.notna()].reset_index(drop=True)
obs = pd.read_csv(D + "/tcga_observed_events_126_292.tsv", sep="\t")
chan = pd.read_csv(P + "/COSMIC_v3_SBS_GRCh37_WESrenorm.txt", sep="\t", index_col=0)
sigs = list(chan.columns)
Ex = pd.read_csv(P + "/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
Ex["patient"] = Ex["Sample Names"].str.slice(0,12); Ex = Ex.drop_duplicates("patient").set_index("patient")
O = pd.read_csv(P + "/context_counts_GRCh37_96_exome.csv", index_col=0).sum(1)
tri = chan.index.str[0] + chan.index.str[2] + chan.index.str[6]
Ovec = O.reindex(tri.values).values.astype(float)
ci = np.array([list(chan.index).index(c) for c in uni.sbs96])

obs = obs[obs.patient.isin(Ex.index)].reset_index(drop=True)
pl = obs.patient.unique()
q = Ex.loc[pl, sigs].astype(float).values @ chan.values.T
Wq = (q / Ovec[None,:])[:, ci]; Pm = Wq / Wq.sum(1, keepdims=True)
pidx = {p:k for k,p in enumerate(pl)}
key = (uni.wt_aa + uni.residue.astype(str) + uni.alt_aa).values
codes, uq = pd.factorize(key)
G = np.zeros((Pm.shape[0], len(uq))); np.add.at(G.T, codes, Pm.T)
G = G / G.sum(1, keepdims=True)
rfsG = pd.Series(uni.CRISPR_RFS.values).groupby(codes).first().values

g = pd.read_excel(GIA, header=1)
g["Mutation"] = g.AA_wt.astype(str) + g.Position.astype(int).astype(str) + g.AA_variant.astype(str)
cols = {"D":"A549_p53WT_Nutlin-3_Z-score", "N":"A549_p53NULL_Nutlin-3_Z-score",
        "E":"A549_p53NULL_Etoposide_Z-score"}
gg = g.assign(**{k: pd.to_numeric(g[v], errors="coerce") for k, v in cols.items()}) \
      .groupby("Mutation")[list(cols)].mean()
M = pd.Series(uq).map(lambda m: True).values
V = gg.reindex(uq)
print("alleles with all three Giacomelli scores: %d / %d" % (V.notna().all(1).sum(), len(uq)))
z = lambda a: (a - np.nanmean(a))/np.nanstd(a)
R  = z(rfsG)
Dv = z(V["D"].values)
Nv = z(V["N"].values)
Ev = z(-V["E"].values)          # sign-flipped: higher = greater disruption
print("correlations among predictors:")
print(pd.DataFrame(np.corrcoef(np.vstack([R, Dv, Nv, Ev])),
                   index=["R","D","N","E"], columns=["R","D","N","E"]).round(3).to_string())

okey = (obs.wt_aa + obs.residue.astype(str) + obs.alt_aa).values
gmap = {k:i for i,k in enumerate(uq)}
gi = np.array([gmap[k] for k in okey]); pi = np.array([pidx[p] for p in obs.patient])
X2 = np.column_stack([R, Dv]); X3 = np.column_stack([R, Dv, Nv, Ev])

def mknll(X, mask=None):
    m = np.ones(len(gi), bool) if mask is None else mask
    def f(b):
        Z = G[pi[m]] * np.exp(X @ np.atleast_1d(b))[None,:]
        return -(np.log(Z[np.arange(m.sum()), gi[m]]) - np.log(Z.sum(1))).sum()
    return f
def fit(X, x0, mask=None):
    r = minimize(mknll(X, mask), x0, method="L-BFGS-B")
    H = r.hess_inv.todense() if hasattr(r.hess_inv, "todense") else np.asarray(r.hess_inv)
    return r.x, -r.fun, np.sqrt(np.diag(H))

b2, ll2, _  = fit(X2, [1.1, 0.9])
b3, ll3, s3 = fit(X3, [1.0, 0.8, 0.2, 0.2])
print("\nM2:\nbeta_R = %.4f\nbeta_D = %.4f\nlogLik = %.3f\nAIC    = %.3f" % (b2[0], b2[1], ll2, -2*ll2+4))
print("\nM3:")
for nm, v, s in zip(["beta_R","beta_D","beta_NullNutlin","beta_NullEtoposide"], b3, s3):
    print("%-20s %8.4f   SE %.4f   95%% CI %8.4f - %.4f" % (nm, v, s, v-1.96*s, v+1.96*s))
print("logLik = %.3f\nAIC    = %.3f" % (ll3, -2*ll3+8))
lr = 2*(ll3-ll2)
print("\nM2 vs M3:\nLR = %.3f\ndf = 2\np  = %.4g\nDelta AIC = %.3f"
      % (lr, chi2.sf(lr,2), (-2*ll3+8)-(-2*ll2+4)))

up = obs.patient.unique().copy(); rng.shuffle(up)
fold = {p: i % NFOLD for i, p in enumerate(up)}
fo = np.array([fold[p] for p in obs.patient])
Pcf2 = np.zeros((len(obs), len(uq))); Pcf3 = np.zeros_like(Pcf2)
h2 = h3 = 0.0; wins = 0
for f in range(NFOLD):
    tr, te = fo != f, fo == f
    a2 = minimize(mknll(X2, tr), [1.1,0.9], method="L-BFGS-B").x
    a3 = minimize(mknll(X3, tr), [1.0,0.8,0.2,0.2], method="L-BFGS-B").x
    l2 = -mknll(X2, te)(a2)/te.sum(); l3 = -mknll(X3, te)(a3)/te.sum()
    h2 += l2/NFOLD; h3 += l3/NFOLD; wins += (l3 > l2)
    for X, a, out in ((X2,a2,Pcf2), (X3,a3,Pcf3)):
        Z = G[pi[te]] * np.exp(X @ a)[None,:]; out[te] = Z/Z.sum(1,keepdims=True)
print("\nMean held-out log-likelihood per event:\nM2 = %.5f\nM3 = %.5f\nDelta = %.5f" % (h2, h3, h3-h2))
print("M3 beats M2 in %d / %d folds" % (wins, NFOLD))

def pb_upper(p, k):
    n = len(p); N = n+1
    w = np.exp(2j*np.pi*np.arange(N)/N)
    pmf = np.real(np.fft.fft(np.prod(1-p[None,:]+p[None,:]*w[:,None], axis=1))/N)
    pmf = np.clip(pmf, 0, None); pmf /= pmf.sum()
    return pmf[k:].sum()
def alleles(Pcf):
    E = Pcf.sum(0); Oc = np.bincount(gi, minlength=len(uq))
    Pv = np.array([pb_upper(Pcf[:,m], Oc[m]) if Oc[m]>0 else 1.0 for m in range(len(uq))])
    o = np.argsort(Pv); r = np.empty(len(Pv),int); r[o]=np.arange(1,len(Pv)+1)
    F = np.empty(len(Pv)); F[o] = np.minimum(np.minimum.accumulate((Pv[o]*len(Pv)/r[o])[::-1])[::-1],1)
    return Oc, E, Pv, F
O2,E2,p2,f2 = alleles(Pcf2); O3,E3,p3,f3 = alleles(Pcf3)
s2 = set(np.array(uq)[(f2<0.05)&(O2>E2)]); s3 = set(np.array(uq)[(f3<0.05)&(O3>E3)])
print("\nFDR < 0.05:\ncross-fitted M2 = %d\ncross-fitted M3 = %d" % (len(s2), len(s3)))
print("FDR < 0.10:\nM2 = %d\nM3 = %d" % (((f2<0.10)&(O2>E2)).sum(), ((f3<0.10)&(O3>E3)).sum()))
print("\nOf the M2 %d: remain under M3 = %d ; lost = %d" % (len(s2), len(s2&s3), len(s2-s3)))
if s3-s2: print("new under M3 = %d: %s" % (len(s3-s2), sorted(s3-s2)))

T = pd.DataFrame(dict(Mutation=uq, Observed=O2, M2_expected=E2, M2_FDR=f2,
                      M3_expected=E3, M3_O_over_E=O3/E3, M3_P=p3, M3_FDR=f3))
T.round(6).to_csv(D + "/residual_crossfit_M3.tsv", sep="\t", index=False)
print("\nMutation  Observed  M2_expected  M2_FDR  M3_expected  M3_O/E  M3_FDR")
print(T.set_index("Mutation").reindex(FOCUS)[["Observed","M2_expected","M2_FDR","M3_expected","M3_O_over_E","M3_FDR"]]
        .to_string(float_format=lambda v: "%.4g" % v))
print("\nTop 20 under M3")
print(T.sort_values("M3_P").head(20)[["Mutation","Observed","M3_expected","M3_O_over_E","M3_P","M3_FDR"]]
        .to_string(index=False, float_format=lambda v: "%.4g" % v))
print("\nSaved: %s/residual_crossfit_M3.tsv" % D)
