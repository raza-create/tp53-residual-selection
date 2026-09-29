#!/usr/bin/env python3
# Step 9: cross-fitted M2 + exact Poisson-binomial allele inference.
import numpy as np, pandas as pd
from scipy.optimize import minimize

P = "/home/raza/p53_thermal_xai/data/external/pcawg_signatures"
D = "/home/raza/p53_thermal_xai/data/external/tp53_selection"
GIA = "/home/raza/p53_thermal_xai/data/external/giacomelli2018/41588_2018_204_MOESM6_ESM.xlsx"
SEED, NFOLD = 20260929, 5
ORDER = ["LUSC","OV","HNSC","BRCA","COADREAD","LUAD","STAD","BLCA","UCEC","GBM","LIHC","PAAD"]
FOCUS = ["Y220C","V157F","I195T","R248Q","R273H","R273L","R273S","G266V",
         "R249S","R175H","R248W","R282W","G245S","R175C","E171K","P278A"]
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
dne = g.assign(DNE=pd.to_numeric(g["A549_p53WT_Nutlin-3_Z-score"], errors="coerce")) \
       .dropna(subset=["DNE"]).groupby("Mutation").DNE.mean()
DNE = pd.Series(uq).map(dne).values
R = (rfsG - rfsG.mean())/rfsG.std(ddof=0)
Dv = (DNE - DNE.mean())/DNE.std(ddof=0)
X = np.column_stack([R, Dv])

okey = (obs.wt_aa + obs.residue.astype(str) + obs.alt_aa).values
gmap = {k:i for i,k in enumerate(uq)}
gi = np.array([gmap[k] for k in okey]); pi = np.array([pidx[p] for p in obs.patient])

up = obs.patient.unique().copy(); rng.shuffle(up)
fold = {p: i % NFOLD for i, p in enumerate(up)}
fo = np.array([fold[p] for p in obs.patient])

def nll_on(mask, b):
    Z = G[pi[mask]] * np.exp(X @ np.atleast_1d(b))[None,:]
    return -(np.log(Z[np.arange(mask.sum()), gi[mask]]) - np.log(Z.sum(1))).sum()

Pcf = np.zeros((len(obs), len(uq)))
print("fold   n_train  n_test   beta_R   beta_D")
for f in range(NFOLD):
    tr, te = fo != f, fo == f
    b = minimize(lambda z: nll_on(tr, z), [1.1, 0.9], method="L-BFGS-B").x
    Z = G[pi[te]] * np.exp(X @ b)[None,:]
    Pcf[te] = Z / Z.sum(1, keepdims=True)
    print("%4d %8d %7d %8.4f %8.4f" % (f, tr.sum(), te.sum(), b[0], b[1]))

E = Pcf.sum(0); Oc = np.bincount(gi, minlength=len(uq))
print("\nCROSS-FITTED M2\n")
print("Events = %d" % len(obs)); print("Protein alleles = %d" % len(uq))

def pb_upper(p, k):
    """P(X >= k) for Poisson-binomial via FFT characteristic function."""
    n = len(p); N = n + 1
    w = np.exp(2j*np.pi*np.arange(N)/N)
    cf = np.prod(1 - p[None,:] + p[None,:]*w[:,None], axis=1)
    pmf = np.real(np.fft.fft(cf)/N)
    pmf = np.clip(pmf, 0, None); pmf /= pmf.sum()
    return pmf[k:].sum()

Pv = np.array([pb_upper(Pcf[:, m], Oc[m]) if Oc[m] > 0 else 1.0 for m in range(len(uq))])
o = np.argsort(Pv); r = np.empty(len(Pv), int); r[o] = np.arange(1, len(Pv)+1)
FDR = np.empty(len(Pv)); FDR[o] = np.minimum(np.minimum.accumulate((Pv[o]*len(Pv)/r[o])[::-1])[::-1], 1)

ct_by = {}
for k, c in zip(gi, obs.cancer.values): ct_by.setdefault(k, set()).add(c)
T = pd.DataFrame(dict(Mutation=uq, Observed=Oc, Crossfit_expected=E, PoissonBinomial_P=Pv, FDR=FDR))
T["O_over_E"] = T.Observed/T.Crossfit_expected
T["N_cancer_types"] = [len(ct_by.get(i,())) for i in range(len(uq))]
T["Cancer_types"] = [",".join(sorted(ct_by.get(i,()), key=ORDER.index)) for i in range(len(uq))]
enr = (T.FDR < 0.05) & (T.Observed > T.Crossfit_expected)
print("\nFDR < 0.05 enriched alleles = %d" % enr.sum())
print("FDR < 0.10 enriched alleles = %d" % (((T.FDR<0.10)&(T.Observed>T.Crossfit_expected)).sum()))

prev = pd.read_csv(D + "/residual_M1_vs_M2.tsv", sep="\t")
old = set(prev.loc[(prev.M2_FDR<0.05)&(prev.M2_OE>1), "Mutation"])
new = set(T.loc[enr, "Mutation"])
print("\nOverlap with original full-fit M2 %d = %d" % (len(old), len(old & new)))
print("Original M2 hits lost = %d" % len(old - new))
print("New hits = %d" % len(new - old))

cols = ["Mutation","Observed","Crossfit_expected","O_over_E","PoissonBinomial_P","FDR","Cancer_types"]
print("\nTop 30 by Poisson-binomial P")
print(T.sort_values(["PoissonBinomial_P","O_over_E"], ascending=[True,False]).head(30)[cols]
        .to_string(index=False, float_format=lambda v: "%.4g" % v))
print("\nPrespecified alleles")
print(T.set_index("Mutation").reindex(FOCUS).reset_index()[cols]
        .to_string(index=False, float_format=lambda v: "%.4g" % v))
T.round(6).to_csv(D + "/residual_crossfit_poissonbinomial.tsv", sep="\t", index=False)
print("\nSaved: %s/residual_crossfit_poissonbinomial.tsv" % D)
