#!/usr/bin/env python3
# Step 13: M4 = M3 + HSF1 potentiation, on the locked SBS1536 supply model and the same five patient folds.
import io, zipfile, numpy as np, pandas as pd, urllib.request
from scipy.optimize import minimize
from scipy.stats import chi2

P   = "/home/raza/p53_thermal_xai/data/external/pcawg_signatures"
D   = "/home/raza/p53_thermal_xai/data/external/tp53_selection"
GIA = "/home/raza/p53_thermal_xai/data/external/giacomelli2018/41588_2018_204_MOESM6_ESM.xlsx"
HSF = "/home/raza/p53_thermal_xai/data/external/halim2026_hsf1/NIHMS2136175-supplement-4.xlsx"
SEED, NFOLD, CHR = 20260929, 5, "17"
ORDER = ["LUSC","OV","HNSC","BRCA","COADREAD","LUAD","STAD","BLCA","UCEC","GBM","LIHC","PAAD"]
SIG = {"LUSC":"Lung-SCC","OV":"Ovary-AdenoCa","HNSC":"Head-SCC","BRCA":"Breast-cancer",
       "COADREAD":"ColoRect-AdenoCa","LUAD":"Lung-AdenoCa","STAD":"Stomach-AdenoCa",
       "BLCA":"Transitional-cell-carcinoma","UCEC":"Uterus-AdenoCa","GBM":"CNS-GBM",
       "LIHC":"Liver-HCC","PAAD":"Panc-AdenoCa"}
COMP = {"A":"T","C":"G","G":"C","T":"A"}
rng = np.random.default_rng(SEED)

uni = pd.read_csv(D + "/tp53_universe_126_292.tsv", sep="\t")
uni = uni[uni.CRISPR_RFS.notna()].reset_index(drop=True)
obs = pd.read_csv(D + "/tcga_observed_events_126_292.tsv", sep="\t")
obs = obs[obs.cancer.isin(ORDER)].reset_index(drop=True)

lo, hi = uni.gpos.min()-3, uni.gpos.max()+3
gseq = urllib.request.urlopen("https://grch37.rest.ensembl.org/sequence/region/human/%s:%d..%d:1?content-type=text/plain"
                              % (CHR, lo, hi), timeout=60).read().decode().strip().upper()
def pent(gp, rp, ap):
    w = gseq[gp-2-lo: gp+3-lo]
    if rp in "CT": return w, "%s>%s" % (rp, ap)
    return "".join(COMP[b] for b in reversed(w)), "%s>%s" % (COMP[rp], COMP[ap])
pv = [pent(r.gpos, r.ref_plus, r.alt_plus) for r in uni.itertuples()]
uni["k1536"] = [b + "|" + a for a, b in pv]

z_ = zipfile.ZipFile(P + "/WES_TCGA_2018_03_09.zip")
c15 = pd.read_csv(io.BytesIO(z_.read("WES_TCGA.1536.csv")))
c15["k1536"] = c15["Mutation type"] + "|" + c15["Pentanucleotide"]
samp = [c for c in c15.columns if "::" in c]
lab = pd.Series([s.split("::")[0] for s in samp], index=samp)
c15 = c15.set_index("k1536")[samp]
O15 = pd.read_csv(P + "/context_counts_GRCh37_1536_exome.csv", index_col=0).sum(1)
mt_of = c15.index.str.split("|").str[0]; pn_of = c15.index.str.split("|").str[1]
Ok = pd.Series(O15.reindex(pn_of).values, index=c15.index).astype(float)
parent = pd.Series((pn_of.str[1] + "[" + mt_of + "]" + pn_of.str[3]).values, index=c15.index)
obs_k = obs.merge(uni[["residue","wt_aa","alt_aa","k1536"]].drop_duplicates(["residue","wt_aa","alt_aa"]),
                  on=["residue","wt_aa","alt_aa"], how="left")
F = {}
for ct in ORDER:
    N = c15[lab[lab == SIG[ct]].index].sum(1).astype(float)
    sub = obs_k[obs_k.cancer == ct].k1536.value_counts()
    N = N.subtract(sub.reindex(N.index).fillna(0), fill_value=0).clip(lower=0)
    df = pd.DataFrame({"N": N, "O": Ok, "par": parent})
    F[ct] = ((df.N + 0.5)/df.groupby("par").N.transform("sum").add(8.0)) / (df.O/df.groupby("par").O.transform("sum"))
Fm = pd.DataFrame(F)

chan = pd.read_csv(P + "/COSMIC_v3_SBS_GRCh37_WESrenorm.txt", sep="\t", index_col=0)
sigs = list(chan.columns)
Ex = pd.read_csv(P + "/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
Ex["patient"] = Ex["Sample Names"].str.slice(0,12); Ex = Ex.drop_duplicates("patient").set_index("patient")
Oex = pd.read_csv(P + "/context_counts_GRCh37_96_exome.csv", index_col=0).sum(1)
tri = chan.index.str[0] + chan.index.str[2] + chan.index.str[6]
Ovec = Oex.reindex(tri.values).values.astype(float)
ci = np.array([list(chan.index).index(c) for c in uni.sbs96])
obs = obs[obs.patient.isin(Ex.index)].reset_index(drop=True)
pl = obs.patient.unique()
r_ic = (Ex.loc[pl, sigs].astype(float).values @ chan.values.T / Ovec[None,:])[:, ci]
pidx = {p:k for k,p in enumerate(pl)}
pat_ct = obs.drop_duplicates("patient").set_index("patient").cancer.to_dict()
cidx = {ct:i for i,ct in enumerate(ORDER)}
Fk = np.column_stack([Fm[ct].reindex(uni.k1536).values for ct in ORDER])
W = r_ic * np.array([Fk[:, cidx[pat_ct[p]]] for p in pl])
Pn = W / W.sum(1, keepdims=True)
key = (uni.wt_aa + uni.residue.astype(str) + uni.alt_aa).values
codes, uq = pd.factorize(key)
G = np.zeros((len(pl), len(uq))); np.add.at(G.T, codes, Pn.T); G = G / G.sum(1, keepdims=True)
rfsG = pd.Series(uni.CRISPR_RFS.values).groupby(codes).first().values

g = pd.read_excel(GIA, header=1)
g["Mutation"] = g.AA_wt.astype(str) + g.Position.astype(int).astype(str) + g.AA_variant.astype(str)
c3 = {"D":"A549_p53WT_Nutlin-3_Z-score","N":"A549_p53NULL_Nutlin-3_Z-score","E":"A549_p53NULL_Etoposide_Z-score"}
V = g.assign(**{k: pd.to_numeric(g[c], errors="coerce") for k,c in c3.items()}).groupby("Mutation")[list(c3)].mean().reindex(uq)
h = pd.read_excel(HSF, sheet_name="mut_log2FC_summary", header=1)
h = h[h.Mut_type == "Missense"]
HC = "mean_nutlin3-dox-37_vs_nutlin3-37"
REPS = ["nutlin3-dox-37_vs_nutlin3-37_rep%d" % i for i in (1,2,3)]
Hm = {c: pd.Series(uq).map(h.groupby("allele")[c].mean()).values for c in [HC]+REPS}
zz = lambda a: (a-np.nanmean(a))/np.nanstd(a)
X3 = np.column_stack([zz(rfsG), zz(V["D"].values), zz(V["N"].values), zz(-V["E"].values)])

okey = (obs.wt_aa + obs.residue.astype(str) + obs.alt_aa).values
gmap = {k:i for i,k in enumerate(uq)}
gi = np.array([gmap[k] for k in okey]); pi = np.array([pidx[p] for p in obs.patient])
def nll(X, mask):
    def f(b):
        Z = G[pi[mask]] * np.exp(X @ np.atleast_1d(b))[None,:]
        return -(np.log(Z[np.arange(mask.sum()), gi[mask]]) - np.log(Z.sum(1))).sum()
    return f
ALL = np.ones(len(gi), bool)
X4 = np.column_stack([X3, zz(Hm[HC])])
r3 = minimize(nll(X3, ALL), [1,.8,-.1,.4], method="L-BFGS-B")
r4 = minimize(nll(X4, ALL), [1,.8,-.1,.4,.1], method="L-BFGS-B")
ll3, ll4 = -r3.fun, -r4.fun
H4 = r4.hess_inv.todense() if hasattr(r4.hess_inv,"todense") else np.asarray(r4.hess_inv)
se = np.sqrt(np.diag(H4))
print("M3\nlogLik = %.3f\nAIC    = %.3f" % (ll3, -2*ll3+8))
print("\nM4")
for n_, v_ in zip(["beta_R","beta_D","beta_N","beta_E","beta_HSF1"], r4.x): print("%-10s = %.4f" % (n_, v_))
print("SE beta_HSF1     = %.4f" % se[4])
print("95%% CI beta_HSF1 = %.4f - %.4f" % (r4.x[4]-1.96*se[4], r4.x[4]+1.96*se[4]))
print("logLik = %.3f\nAIC    = %.3f" % (ll4, -2*ll4+10))
lr = 2*(ll4-ll3)
print("\nM3 vs M4\nLR = %.3f\ndf = 1\np  = %.4g\nDelta AIC = %.3f" % (lr, chi2.sf(lr,1), (-2*ll4+10)-(-2*ll3+8)))

up = obs.patient.unique().copy(); rng.shuffle(up)
fold = {p:i%NFOLD for i,p in enumerate(up)}; fo = np.array([fold[p] for p in obs.patient])
def cf(X, x0):
    Pcf = np.zeros((len(obs), len(uq))); per = []
    for f_ in range(NFOLD):
        tr, te = fo != f_, fo == f_
        a = minimize(nll(X, tr), x0, method="L-BFGS-B").x
        per.append((a, -nll(X, te)(a)/te.sum()))
        Z = G[pi[te]] * np.exp(X @ a)[None,:]; Pcf[te] = Z/Z.sum(1,keepdims=True)
    return Pcf, per
P3, p3f = cf(X3, [1,.8,-.1,.4]); P4, p4f = cf(X4, [1,.8,-.1,.4,.1])
print("\nfold   beta_HSF1   heldout_LL_M3   heldout_LL_M4   delta_LL")
for i,((_,l3),(a4,l4)) in enumerate(zip(p3f,p4f)):
    print("%4d %11.4f %15.5f %15.5f %11.5f" % (i, a4[4], l3, l4, l4-l3))
m3 = np.mean([x[1] for x in p3f]); m4 = np.mean([x[1] for x in p4f])
print("\nMean held-out LL/event:\nM3 = %.5f\nM4 = %.5f\nDelta = %.5f" % (m3, m4, m4-m3))
print("M4 better than M3 in %d / %d folds" % (sum(l4>l3 for (_,l3),(_,l4) in zip(p3f,p4f)), NFOLD))

def pb(p,k):
    N=len(p)+1; w=np.exp(2j*np.pi*np.arange(N)/N)
    pmf=np.real(np.fft.fft(np.prod(1-p[None,:]+p[None,:]*w[:,None],axis=1))/N)
    pmf=np.clip(pmf,0,None); pmf/=pmf.sum(); return pmf[k:].sum()
def alle(Pcf):
    E=Pcf.sum(0); Oc=np.bincount(gi,minlength=len(uq))
    pv=np.array([pb(Pcf[:,m],Oc[m]) if Oc[m]>0 else 1.0 for m in range(len(uq))])
    o=np.argsort(pv); r=np.empty(len(pv),int); r[o]=np.arange(1,len(pv)+1)
    F_=np.empty(len(pv)); F_[o]=np.minimum(np.minimum.accumulate((pv[o]*len(pv)/r[o])[::-1])[::-1],1)
    return Oc,E,pv,F_
O3,E3,q3,f3 = alle(P3); O4,E4,q4,f4 = alle(P4)
s3=set(np.array(uq)[(f3<0.05)&(O3>E3)]); s4=set(np.array(uq)[(f4<0.05)&(O4>E4)])
print("\nSBS1536-M3 FDR<0.05 = %d\nSBS1536-M4 FDR<0.05 = %d" % (len(s3), len(s4)))
print("\nOf the M3 %d: remain = %d ; lost = %d ; new = %d" % (len(s3), len(s3&s4), len(s3-s4), len(s4-s3)))
if s3-s4: print("  lost: %s" % sorted(s3-s4))
if s4-s3: print("  new : %s" % sorted(s4-s3))
ctb = {}
for k,c in zip(gi, obs.cancer.values): ctb.setdefault(k,set()).add(c)
T = pd.DataFrame(dict(Mutation=uq, Observed=O4, Expected_M4=E4, OE_M4=O4/E4, P_M4=q4, FDR_M4=f4,
                      Expected_M3=E3, FDR_M3=f3, HSF1=Hm[HC],
                      N_ct=[len(ctb.get(i,())) for i in range(len(uq))]))
T.round(6).to_csv(D + "/residual_M4.tsv", sep="\t", index=False)
print("\nTop 30 M4 residual alleles")
print(T.sort_values("P_M4").head(30)[["Mutation","Observed","Expected_M4","OE_M4","P_M4","FDR_M4","N_ct","HSF1"]]
        .to_string(index=False, float_format=lambda v: "%.4g" % v))
print("\nHSF1 representation       beta_H    mean held-out LL/event")
for nm, col in [("Mean of 3",HC),("Replicate 1",REPS[0]),("Replicate 2",REPS[1]),("Replicate 3",REPS[2])]:
    Xr = np.column_stack([X3, zz(Hm[col])])
    _, pr = cf(Xr, [1,.8,-.1,.4,.1])
    bfull = minimize(nll(Xr, ALL), [1,.8,-.1,.4,.1], method="L-BFGS-B").x
    print("%-22s %8.4f %19.5f" % (nm, bfull[4], np.mean([x[1] for x in pr])))
print("\nSaved: %s/residual_M4.tsv" % D)
