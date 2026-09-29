#!/usr/bin/env python3
# Step 8: does dominant-negative activity (Giacomelli WT-Nutlin) explain the residual selection?
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import chi2

P = "/home/raza/p53_thermal_xai/data/external/pcawg_signatures"
D = "/home/raza/p53_thermal_xai/data/external/tp53_selection"
GIA = "/home/raza/p53_thermal_xai/data/external/giacomelli2018/41588_2018_204_MOESM6_ESM.xlsx"
NB, SEED, NFOLD = 100000, 20260929, 5
FOCUS = ["P278A","R273S","G266V","R273L","I195T","Y220C","V157F","R249S",
         "R248Q","R175H","R248W","R282W","R273H","G245S","R175C","E171K"]
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
Gsup = np.zeros((Pm.shape[0], len(uq))); np.add.at(Gsup.T, codes, Pm.T)
rfsG = pd.Series(uni.CRISPR_RFS.values).groupby(codes).first().values

g = pd.read_excel(GIA, header=1)
g["Mutation"] = g.AA_wt.astype(str) + g.Position.astype(int).astype(str) + g.AA_variant.astype(str)
g["DNE"] = pd.to_numeric(g["A549_p53WT_Nutlin-3_Z-score"], errors="coerce")
dne_map = g.dropna(subset=["DNE"]).groupby("Mutation").DNE.mean()
DNE = pd.Series(uq).map(dne_map).values
ok = ~np.isnan(DNE)
print("Candidate protein alleles with both scores = %d (of %d)" % (ok.sum(), len(uq)))

okey = (obs.wt_aa + obs.residue.astype(str) + obs.alt_aa).values
gmap = {k:i for i,k in enumerate(uq)}
gi_all = np.array([gmap.get(k, -1) for k in okey])
keepev = (gi_all >= 0) & ok[np.maximum(gi_all,0)]
obs2 = obs[keepev].reset_index(drop=True)
print("Observed TCGA events retained              = %d (of %d)" % (len(obs2), len(obs)))

sub = np.where(ok)[0]
Gs = Gsup[:, sub]; Gs = Gs / Gs.sum(1, keepdims=True)
remap = {v:i for i,v in enumerate(sub)}
gi = np.array([remap[g_] for g_ in gi_all[keepev]])
pi = np.array([pidx[p] for p in obs2.patient])
R = (rfsG[sub] - rfsG[sub].mean())/rfsG[sub].std(ddof=0)
Dv = (DNE[sub] - np.nanmean(DNE[sub]))/np.nanstd(DNE[sub])
# orient DNE: higher = stronger dominant-negative (Giacomelli WT-Nutlin enrichment = DN)
print("DNE orientation: corr(DNE_z, RFS_z) = %.3f" % np.corrcoef(Dv, R)[0,1])

def mk(X):
    def nll(b):
        Z = Gs[pi] * np.exp(X @ np.atleast_1d(b))[None,:]
        return -(np.log(Z[np.arange(len(gi)), gi]) - np.log(Z.sum(1))).sum()
    return nll
def fit(X, x0):
    r = minimize(mk(X), x0, method="L-BFGS-B")
    return r.x, -r.fun, r
X1 = R[:,None]; X2 = np.column_stack([R, Dv])
b1, ll1, _ = fit(X1, [2.0]); b2, ll2, r2 = fit(X2, [2.0, 0.0])
H = r2.hess_inv.todense() if hasattr(r2.hess_inv,"todense") else np.linalg.inv(r2.hess_inv)
se = np.sqrt(np.diag(H))
print("\nM1: mutation supply + CRISPR RFS")
print("beta_R = %.4f\nlogLik = %.3f\nAIC    = %.3f" % (b1[0], ll1, -2*ll1+2))
print("\nM2: mutation supply + CRISPR RFS + DNE")
print("beta_R      = %.4f\nbeta_D      = %.4f\nSE beta_D   = %.4f\n95%% CI beta_D = %.4f - %.4f"
      % (b2[0], b2[1], se[1], b2[1]-1.96*se[1], b2[1]+1.96*se[1]))
print("logLik = %.3f\nAIC    = %.3f" % (ll2, -2*ll2+4))
lr = 2*(ll2-ll1)
print("\nLikelihood-ratio test M1 vs M2:\nLR = %.3f\ndf = 1\np  = %.4g" % (lr, chi2.sf(lr,1)))
print("Delta AIC = %.3f" % ((-2*ll2+4)-(-2*ll1+2)))

up = obs2.patient.unique(); rng.shuffle(up)
fold = {p: i % NFOLD for i, p in enumerate(up)}
fo = np.array([fold[p] for p in obs2.patient])
h1 = h2 = 0; wins = 0
for f in range(NFOLD):
    tr, te = fo != f, fo == f
    def sub_nll(X, mask):
        def nll(b):
            Z = Gs[pi[mask]] * np.exp(X @ np.atleast_1d(b))[None,:]
            return -(np.log(Z[np.arange(mask.sum()), gi[mask]]) - np.log(Z.sum(1))).sum()
        return nll
    a1 = minimize(sub_nll(X1, tr), [2.0], method="L-BFGS-B").x
    a2 = minimize(sub_nll(X2, tr), [2.0, 0.0], method="L-BFGS-B").x
    l1 = -sub_nll(X1, te)(a1)/te.sum(); l2 = -sub_nll(X2, te)(a2)/te.sum()
    h1 += l1/NFOLD; h2 += l2/NFOLD; wins += (l2 > l1)
print("\nMean held-out log-likelihood (per event):\nM1 = %.5f\nM2 = %.5f" % (h1, h2))
print("Delta held-out log-likelihood = %.5f" % (h2-h1))
print("Number of folds where M2 beats M1 = %d / %d" % (wins, NFOLD))

def residual(X, b):
    Z = Gs[pi] * np.exp(X @ np.atleast_1d(b))[None,:]
    Pf = Z/Z.sum(1,keepdims=True); E = Pf.sum(0)
    Oc = np.bincount(gi, minlength=len(sub))
    cum = np.cumsum(Pf,axis=1); ge = np.zeros(len(sub), np.int64)
    for _ in range(NB):
        u = rng.random(len(gi))
        ge += (np.bincount(np.minimum((cum<u[:,None]).sum(1), len(sub)-1), minlength=len(sub)) >= Oc)
    p = (1+ge)/(NB+1)
    o = np.argsort(p); r = np.empty(len(p),int); r[o]=np.arange(1,len(p)+1)
    fdr = np.empty(len(p)); fdr[o] = np.minimum(np.minimum.accumulate((p[o]*len(p)/r[o])[::-1])[::-1],1)
    return Oc, E, p, fdr
O1,E1,p1,f1 = residual(X1,b1); O2,E2,p2,f2 = residual(X2,b2)
n1 = ((f1<0.05)&(O1>E1)).sum(); n2 = ((f2<0.05)&(O2>E2)).sum()
print("\nFDR<0.05 residual-enriched alleles under M1 = %d" % n1)
print("FDR<0.05 residual-enriched alleles under M2 = %d" % n2)
s1 = set(np.where((f1<0.05)&(O1>E1))[0]); s2 = set(np.where((f2<0.05)&(O2>E2))[0])
print("Of the M1 set: remain significant under M2 = %d ; no longer significant = %d"
      % (len(s1&s2), len(s1-s2)))
nm = np.array(uq)[sub]
T = pd.DataFrame(dict(Mutation=nm, M1_OE=O1/E1, M1_FDR=f1, M2_OE=O2/E2, M2_FDR=f2, Observed=O1))
T.round(5).to_csv(D + "/residual_M1_vs_M2.tsv", sep="\t", index=False)
print("\nMutation   M1_O/E   M1_FDR   M2_O/E   M2_FDR")
print(T.set_index("Mutation").reindex(FOCUS)[["M1_OE","M1_FDR","M2_OE","M2_FDR"]]
        .to_string(float_format=lambda v: "%.4g" % v))
print("\nSaved: %s/residual_M1_vs_M2.tsv" % D)
