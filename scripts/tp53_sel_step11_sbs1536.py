#!/usr/bin/env python3
# Step 11: add cancer-specific pentanucleotide (1536) mutability modifier to the supply model; refit cross-fitted M3.
import io, zipfile, numpy as np, pandas as pd, urllib.request, json
from scipy.optimize import minimize

P = "/home/raza/p53_thermal_xai/data/external/pcawg_signatures"
D = "/home/raza/p53_thermal_xai/data/external/tp53_selection"
GIA = "/home/raza/p53_thermal_xai/data/external/giacomelli2018/41588_2018_204_MOESM6_ESM.xlsx"
SEED, NFOLD, CHR = 20260929, 5, "17"
ORDER = ["LUSC","OV","HNSC","BRCA","COADREAD","LUAD","STAD","BLCA","UCEC","GBM","LIHC","PAAD"]
SIG = {"LUSC":"Lung-SCC","OV":"Ovary-AdenoCa","HNSC":"Head-SCC","BRCA":"Breast-cancer",
       "COADREAD":"ColoRect-AdenoCa","LUAD":"Lung-AdenoCa","STAD":"Stomach-AdenoCa",
       "BLCA":"Transitional-cell-carcinoma","UCEC":"Uterus-AdenoCa","GBM":"CNS-GBM",
       "LIHC":"Liver-HCC","PAAD":"Panc-AdenoCa"}
FOCUS = ["Y220C","V157F","R158L","I195T","R248Q","R273L","R273H","R273S","G266V",
         "H179R","R249S","R175H","G245V","G245D","C242F"]
COMP = {"A":"T","C":"G","G":"C","T":"A"}
rng = np.random.default_rng(SEED)

uni = pd.read_csv(D + "/tp53_universe_126_292.tsv", sep="\t")
uni = uni[uni.CRISPR_RFS.notna()].reset_index(drop=True)
obs = pd.read_csv(D + "/tcga_observed_events_126_292.tsv", sep="\t")
obs = obs[obs.cancer.isin(ORDER)].reset_index(drop=True)

lo, hi = uni.gpos.min()-3, uni.gpos.max()+3
gseq = urllib.request.urlopen("https://grch37.rest.ensembl.org/sequence/region/human/%s:%d..%d:1?content-type=text/plain"
                              % (CHR, lo, hi), timeout=60).read().decode().strip().upper()
def pent(gp, ref_plus, alt_plus):
    w = gseq[gp-2-lo: gp+3-lo]
    if ref_plus in "CT": return w, "%s>%s" % (ref_plus, alt_plus)
    return "".join(COMP[b] for b in reversed(w)), "%s>%s" % (COMP[ref_plus], COMP[alt_plus])
pv = [pent(r.gpos, r.ref_plus, r.alt_plus) for r in uni.itertuples()]
uni["pent"] = [a for a, _ in pv]; uni["mt"] = [b for _, b in pv]
bad = uni[uni.pent.str[2] != uni.mt.str[0]]
print("pentanucleotide centre-base check failures: %d" % len(bad))
uni["k1536"] = uni.mt + "|" + uni.pent

z = zipfile.ZipFile(P + "/WES_TCGA_2018_03_09.zip")
c15 = pd.read_csv(io.BytesIO(z.read("WES_TCGA.1536.csv")))
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
print("\nBackground mutations used per cancer (TP53 events of that cancer subtracted):")
F = {}
for ct in ORDER:
    cols = lab[lab == SIG[ct]].index
    N = c15[cols].sum(1).astype(float)
    sub = obs_k[obs_k.cancer == ct].k1536.value_counts()
    N = N.subtract(sub.reindex(N.index).fillna(0), fill_value=0).clip(lower=0)
    print("%-9s = %s" % (ct, f"{int(N.sum()):,}"))
    df = pd.DataFrame({"N": N, "O": Ok, "par": parent})
    df["Phat"] = (df.N + 0.5) / df.groupby("par").N.transform("sum").add(8.0)
    df["q"] = df.O / df.groupby("par").O.transform("sum")
    F[ct] = (df.Phat / df.q)
Fm = pd.DataFrame(F)
print("\nF_{t,k} across cancers: median %.3f, 5-95%% %.3f-%.3f, max %.3f"
      % (np.median(Fm.values), np.percentile(Fm.values,5), np.percentile(Fm.values,95), Fm.values.max()))

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
q96 = Ex.loc[pl, sigs].astype(float).values @ chan.values.T
r_ic = (q96 / Ovec[None,:])[:, ci]
pidx = {p:k for k,p in enumerate(pl)}

key = (uni.wt_aa + uni.residue.astype(str) + uni.alt_aa).values
codes, uq = pd.factorize(key)
rfsG = pd.Series(uni.CRISPR_RFS.values).groupby(codes).first().values
g = pd.read_excel(GIA, header=1)
g["Mutation"] = g.AA_wt.astype(str) + g.Position.astype(int).astype(str) + g.AA_variant.astype(str)
cols3 = {"D":"A549_p53WT_Nutlin-3_Z-score","N":"A549_p53NULL_Nutlin-3_Z-score","E":"A549_p53NULL_Etoposide_Z-score"}
V = g.assign(**{k: pd.to_numeric(g[v], errors="coerce") for k,v in cols3.items()}).groupby("Mutation")[list(cols3)].mean().reindex(uq)
zs = lambda a: (a-np.nanmean(a))/np.nanstd(a)
X3 = np.column_stack([zs(rfsG), zs(V["D"].values), zs(V["N"].values), zs(-V["E"].values)])

pat_ct = obs.drop_duplicates("patient").set_index("patient").cancer.to_dict()
Fk = np.column_stack([Fm[ct].reindex(uni.k1536).values for ct in ORDER])
cidx = {ct:i for i,ct in enumerate(ORDER)}
def build(use1536):
    W = r_ic.copy()
    if use1536:
        fr = np.array([Fk[:, cidx[pat_ct[p]]] for p in pl])
        W = W * fr
    Pn = W / W.sum(1, keepdims=True)
    G = np.zeros((len(pl), len(uq))); np.add.at(G.T, codes, Pn.T)
    return G / G.sum(1, keepdims=True)
G96, G15 = build(False), build(True)

okey = (obs.wt_aa + obs.residue.astype(str) + obs.alt_aa).values
gmap = {k:i for i,k in enumerate(uq)}
gi = np.array([gmap[k] for k in okey]); pi = np.array([pidx[p] for p in obs.patient])
up = obs.patient.unique().copy(); rng.shuffle(up)
fold = {p:i%NFOLD for i,p in enumerate(up)}; fo = np.array([fold[p] for p in obs.patient])

def crossfit(G):
    Pcf = np.zeros((len(obs), len(uq))); h = 0.0
    def nll(mask):
        def f(b):
            Z = G[pi[mask]] * np.exp(X3 @ b)[None,:]
            return -(np.log(Z[np.arange(mask.sum()), gi[mask]]) - np.log(Z.sum(1))).sum()
        return f
    for f_ in range(NFOLD):
        tr, te = fo != f_, fo == f_
        a = minimize(nll(tr), [1.0,0.8,-0.1,0.4], method="L-BFGS-B").x
        h += (-nll(te)(a)/te.sum())/NFOLD
        Z = G[pi[te]] * np.exp(X3 @ a)[None,:]; Pcf[te] = Z/Z.sum(1,keepdims=True)
    return Pcf, h
def pb(p, k):
    N = len(p)+1; w = np.exp(2j*np.pi*np.arange(N)/N)
    pmf = np.real(np.fft.fft(np.prod(1-p[None,:]+p[None,:]*w[:,None],axis=1))/N)
    pmf = np.clip(pmf,0,None); pmf/=pmf.sum(); return pmf[k:].sum()
def alle(Pcf):
    E = Pcf.sum(0); Oc = np.bincount(gi, minlength=len(uq))
    pv = np.array([pb(Pcf[:,m], Oc[m]) if Oc[m]>0 else 1.0 for m in range(len(uq))])
    o = np.argsort(pv); r = np.empty(len(pv),int); r[o]=np.arange(1,len(pv)+1)
    F_ = np.empty(len(pv)); F_[o] = np.minimum(np.minimum.accumulate((pv[o]*len(pv)/r[o])[::-1])[::-1],1)
    return Oc, E, pv, F_
Pa, ha = crossfit(G96); Pb_, hb = crossfit(G15)
print("\nSBS1536 MUTABILITY ROBUSTNESS\n")
print("Cross-fitted held-out log-likelihood/event:")
print("SBS96 + M3   = %.5f" % ha); print("SBS1536 + M3 = %.5f" % hb); print("Delta        = %.5f" % (hb-ha))
O1,E1,p1,f1 = alle(Pa); O2,E2,p2,f2 = alle(Pb_)
s1 = set(np.array(uq)[(f1<0.05)&(O1>E1)]); s2 = set(np.array(uq)[(f2<0.05)&(O2>E2)])
print("\nResidual enrichment:")
print("SBS96-M3 FDR<0.05   = %d" % len(s1)); print("SBS1536-M3 FDR<0.05 = %d" % len(s2))
print("\nOf original %d: remain = %d ; lost = %d ; new = %d" % (len(s1), len(s1&s2), len(s1-s2), len(s2-s1)))
if s1-s2: print("  lost: %s" % sorted(s1-s2))
if s2-s1: print("  new : %s" % sorted(s2-s1))
T = pd.DataFrame(dict(Mutation=uq, Observed=O1, Expected_96=E1, OE_96=O1/E1, FDR_96=f1,
                      Expected_1536=E2, OE_1536=O2/E2, FDR_1536=f2))
T.round(6).to_csv(D + "/residual_sbs1536.tsv", sep="\t", index=False)
print("\n" + T.set_index("Mutation").reindex(FOCUS).to_string(float_format=lambda v: "%.4g" % v))
print("\nTop 20 under SBS1536-M3")
print(T.sort_values("FDR_1536").head(20).to_string(index=False, float_format=lambda v: "%.4g" % v))
print("\nSaved: %s/residual_sbs1536.tsv" % D)
