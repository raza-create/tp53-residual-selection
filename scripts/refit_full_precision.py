#!/usr/bin/env python3
# Refit M4 (SBS1536) to tight tolerance, store full-precision coefficients + exact SEs.
import io, zipfile, json, numpy as np, pandas as pd, urllib.request
from pathlib import Path
from scipy.optimize import minimize
P=Path("/home/raza/p53_thermal_xai/data/external/pcawg_signatures")
B=Path("/home/raza/p53_thermal_xai/data/external/tp53_selection")
GIA="/home/raza/p53_thermal_xai/data/external/giacomelli2018/41588_2018_204_MOESM6_ESM.xlsx"
HSF="/home/raza/p53_thermal_xai/data/external/halim2026_hsf1/NIHMS2136175-supplement-4.xlsx"
CHR="17"; COMP={"A":"T","C":"G","G":"C","T":"A"}
ORDER=["LUSC","OV","HNSC","BRCA","COADREAD","LUAD","STAD","BLCA","UCEC","GBM","LIHC","PAAD"]
SIG={"LUSC":"Lung-SCC","OV":"Ovary-AdenoCa","HNSC":"Head-SCC","BRCA":"Breast-cancer",
     "COADREAD":"ColoRect-AdenoCa","LUAD":"Lung-AdenoCa","STAD":"Stomach-AdenoCa",
     "BLCA":"Transitional-cell-carcinoma","UCEC":"Uterus-AdenoCa","GBM":"CNS-GBM",
     "LIHC":"Liver-HCC","PAAD":"Panc-AdenoCa"}
uni=pd.read_csv(B/"tp53_universe_126_292.tsv",sep="\t"); uni=uni[uni.CRISPR_RFS.notna()].reset_index(drop=True)
lo,hi=uni.gpos.min()-3,uni.gpos.max()+3
gs=urllib.request.urlopen("https://grch37.rest.ensembl.org/sequence/region/human/%s:%d..%d:1?content-type=text/plain"%(CHR,lo,hi),timeout=60).read().decode().strip().upper()
def pent(gp,rp,ap):
    w=gs[gp-2-lo:gp+3-lo]
    return (w,"%s>%s"%(rp,ap)) if rp in "CT" else ("".join(COMP[b] for b in reversed(w)),"%s>%s"%(COMP[rp],COMP[ap]))
uni["k1536"]=[b+"|"+a for a,b in (pent(r.gpos,r.ref_plus,r.alt_plus) for r in uni.itertuples())]
uni["Mutation"]=uni.wt_aa+uni.residue.astype(str)+uni.alt_aa
obs=pd.read_csv(B/"tcga_observed_events_126_292.tsv",sep="\t"); obs=obs[obs.cancer.isin(ORDER)].reset_index(drop=True)
zf=zipfile.ZipFile(P/"WES_TCGA_2018_03_09.zip")
c15=pd.read_csv(io.BytesIO(zf.read("WES_TCGA.1536.csv"))); c15["k"]=c15["Mutation type"]+"|"+c15["Pentanucleotide"]
sT=[c for c in c15.columns if "::" in c]; lab=pd.Series([s.split("::")[0] for s in sT],index=sT); c15=c15.set_index("k")[sT]
O15=pd.read_csv(P/"context_counts_GRCh37_1536_exome.csv",index_col=0).sum(1)
pn=c15.index.str.split("|").str[1]; mt=c15.index.str.split("|").str[0]
Ok=pd.Series(O15.reindex(pn).values,index=c15.index).astype(float)
par=pd.Series((pn.str[1]+"["+mt+"]"+pn.str[3]).values,index=c15.index)
ok2=obs.merge(uni[["residue","wt_aa","alt_aa","k1536"]].drop_duplicates(["residue","wt_aa","alt_aa"]),on=["residue","wt_aa","alt_aa"],how="left")
F={}
for ct in ORDER:
    N=c15[lab[lab==SIG[ct]].index].sum(1).astype(float)
    N=N.subtract(ok2[ok2.cancer==ct].k1536.value_counts().reindex(N.index).fillna(0),fill_value=0).clip(lower=0)
    d=pd.DataFrame({"N":N,"O":Ok,"par":par})
    F[ct]=((d.N+0.5)/d.groupby("par").N.transform("sum").add(8.0))/(d.O/d.groupby("par").O.transform("sum"))
Fm=pd.DataFrame(F)
chan=pd.read_csv(P/"COSMIC_v3_SBS_GRCh37_WESrenorm.txt",sep="\t",index_col=0); sigs=list(chan.columns)
Ex=pd.read_csv(P/"TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
Ex["patient"]=Ex["Sample Names"].str.slice(0,12); Ex=Ex.drop_duplicates("patient").set_index("patient")
Oex=pd.read_csv(P/"context_counts_GRCh37_96_exome.csv",index_col=0).sum(1)
tri=chan.index.str[0]+chan.index.str[2]+chan.index.str[6]; Ovec=Oex.reindex(tri.values).values.astype(float)
ci=np.array([list(chan.index).index(c) for c in uni.sbs96])
obs=obs[obs.patient.isin(Ex.index)].reset_index(drop=True); pl=obs.patient.unique()
r_ic=(Ex.loc[pl,sigs].astype(float).values@chan.values.T/Ovec[None,:])[:,ci]
pidx={p:i for i,p in enumerate(pl)}; pct=obs.drop_duplicates("patient").set_index("patient").cancer.to_dict()
cidx={c:i for i,c in enumerate(ORDER)}; Fk=np.column_stack([Fm[c].reindex(uni.k1536).values for c in ORDER])
codes,uq=pd.factorize(uni.Mutation.values)
W=r_ic*np.array([Fk[:,cidx[pct[p]]] for p in pl]); Pn=W/W.sum(1,keepdims=True)
G=np.zeros((len(pl),len(uq))); np.add.at(G.T,codes,Pn.T); G=G/G.sum(1,keepdims=True)
rfsG=pd.Series(uni.CRISPR_RFS.values).groupby(codes).first().values
g=pd.read_excel(GIA,header=1); g["Mutation"]=g.AA_wt.astype(str)+g.Position.astype(int).astype(str)+g.AA_variant.astype(str)
c3={"D":"A549_p53WT_Nutlin-3_Z-score","N":"A549_p53NULL_Nutlin-3_Z-score","E":"A549_p53NULL_Etoposide_Z-score"}
V=g.assign(**{k:pd.to_numeric(g[c],errors="coerce") for k,c in c3.items()}).groupby("Mutation")[list(c3)].mean().reindex(uq)
h=pd.read_excel(HSF,sheet_name="mut_log2FC_summary",header=1); h=h[h.Mut_type=="Missense"]
Hraw=pd.Series(uq).map(h.groupby("allele")["mean_nutlin3-dox-37_vs_nutlin3-37"].mean()).values
z=lambda a:(a-np.nanmean(a))/np.nanstd(a)
X=np.column_stack([z(rfsG),z(V["D"].values),z(V["N"].values),z(-V["E"].values),z(Hraw)])
gm={k:i for i,k in enumerate(uq)}
gi=np.array([gm[k] for k in (obs.wt_aa+obs.residue.astype(str)+obs.alt_aa)]); pi=np.array([pidx[p] for p in obs.patient])
def nll(b):
    Z=G[pi]*np.exp(X@b)[None,:]
    return -(np.log(Z[np.arange(len(gi)),gi])-np.log(Z.sum(1))).sum()
best=None
for s0 in ([1,.8,-.1,.4,.1],[0,0,0,0,0],[1.0031,0.7109,-0.0912,0.4530,0.0876]):
    r=minimize(nll,s0,method="L-BFGS-B",options={"ftol":1e-15,"gtol":1e-12,"maxiter":10000})
    if best is None or r.fun<best.fun: best=r
b=best.x
Z=G[pi]*np.exp(X@b)[None,:]; Pb=Z/Z.sum(1,keepdims=True); m1=Pb@X
I=np.zeros((5,5))
for a_ in range(5):
    for c_ in range(5):
        I[a_,c_]=((Pb*(X[:,a_]*X[:,c_])[None,:]).sum(1)-m1[:,a_]*m1[:,c_]).sum()
se=np.sqrt(np.diag(np.linalg.inv(I)))
nm=["R","D","N","E","H"]
out=dict(m4_coefficients_full={k:float(v) for k,v in zip(nm,b)},
         m4_exact_se_full={k:float(v) for k,v in zip(nm,se)},
         m4_loglik_full=float(-best.fun), max_gradient=float(np.abs(best.jac).max()))
print(json.dumps(out,indent=2))
prev=dict(R=1.0031,D=0.7109,N=-0.0912,E=0.4530,H=0.0876)
print("\ndelta vs frozen-as-used:")
for k,v in zip(nm,b): print("  %s  %+.8f  vs %+.4f   diff %.2e" % (k,v,prev[k],v-prev[k]))
json.dump(out,open("/home/raza/p53_thermal_xai/config/m4_full_precision.json","w"),indent=2)
