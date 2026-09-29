#!/usr/bin/env python3
# Step 7: residual allele-specific selection after sample-specific mutation supply + global LOF selection.
import numpy as np, pandas as pd

P = "/home/raza/p53_thermal_xai/data/external/pcawg_signatures"
D = "/home/raza/p53_thermal_xai/data/external/tp53_selection"
BETA, NB, SEED = 2.3416, 100000, 20260929
ORDER = ["LUSC","OV","HNSC","BRCA","COADREAD","LUAD","STAD","BLCA","UCEC","GBM","LIHC","PAAD"]
FOCUS = ["R175H","G245S","R248Q","R248W","R249S","R273H","R282W","Y220C"]
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
rfsG = pd.Series(uni.CRISPR_RFS.values).groupby(codes).first().values

okey = (obs.wt_aa + obs.residue.astype(str) + obs.alt_aa).values
gmap = {k:i for i,k in enumerate(uq)}
keep = np.array([k in gmap for k in okey]); obs, okey = obs[keep].reset_index(drop=True), okey[keep]
pi = np.array([pidx[p] for p in obs.patient]); gi = np.array([gmap[k] for k in okey])

Z = G[pi] * np.exp(BETA*rfsG)[None,:]; Pfit = Z / Z.sum(1, keepdims=True)
E = Pfit.sum(0)
Ocount = np.bincount(gi, minlength=len(uq))
print("events = %d | alleles tested = %d | sum(E) = %.1f" % (len(obs), len(uq), E.sum()))

cum = np.cumsum(Pfit, axis=1)
ge = np.zeros(len(uq), dtype=np.int64)
for b in range(NB):
    u = rng.random(len(obs))
    pick = np.minimum((cum < u[:,None]).sum(1), len(uq)-1)
    c = np.bincount(pick, minlength=len(uq))
    ge += (c >= Ocount)
Pval = (1+ge)/(NB+1)

o = np.argsort(Pval); r = np.empty(len(uq), int); r[o] = np.arange(1, len(uq)+1)
fdr = np.minimum.accumulate((Pval[o]*len(uq)/r[o])[::-1])[::-1]
FDR = np.empty(len(uq)); FDR[o] = np.minimum(fdr, 1)

ct_by = {}
for k, c in zip(gi, obs.cancer.values): ct_by.setdefault(k, set()).add(c)
T = pd.DataFrame(dict(Mutation=uq, Observed=Ocount, Expected=E, RFS=rfsG,
                      empirical_P=Pval, FDR=FDR))
T["Observed_over_Expected"] = T.Observed/T.Expected
T["N_cancer_types"] = [len(ct_by.get(i,())) for i in range(len(uq))]
T["Cancer_types"] = [",".join(sorted(ct_by.get(i,()), key=ORDER.index)) for i in range(len(uq))]
T.round(5).to_csv(D + "/residual_selection_alleles.tsv", sep="\t", index=False)

print("\nNumber of alleles tested          = %d" % len(T))
for t in (2,5,10): print("Alleles with observed count >= %-2d = %d" % (t, (T.Observed>=t).sum()))
print("\nFDR < 0.05 residual-enriched alleles = %d" % ((T.FDR<0.05)&(T.Observed>T.Expected)).sum())
print("FDR < 0.10 residual-enriched alleles = %d" % ((T.FDR<0.10)&(T.Observed>T.Expected)).sum())

cols = ["Mutation","Observed","Expected","Observed_over_Expected","empirical_P","FDR","N_cancer_types","RFS","Cancer_types"]
print("\nTop 30 by empirical P (ties broken by O/E)")
print(T.sort_values(["empirical_P","Observed_over_Expected"], ascending=[True,False]).head(30)[cols]
        .to_string(index=False, float_format=lambda v: "%.4g" % v))
print("\nPrespecified hotspot alleles")
print(T[T.Mutation.isin(FOCUS)].set_index("Mutation").reindex(FOCUS).reset_index()[cols]
        .to_string(index=False, float_format=lambda v: "%.4g" % v))
print("\nSaved: %s/residual_selection_alleles.tsv" % D)
