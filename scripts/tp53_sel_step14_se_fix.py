#!/usr/bin/env python3
# Step 14: reconcile Wald SE with the LRT for beta_HSF1. Exact observed information + profile likelihood.
import io, zipfile, numpy as np, pandas as pd, urllib.request
from scipy.optimize import minimize, brentq
from scipy.stats import chi2, norm

P="/home/raza/p53_thermal_xai/data/external/pcawg_signatures"
D="/home/raza/p53_thermal_xai/data/external/tp53_selection"
GIA="/home/raza/p53_thermal_xai/data/external/giacomelli2018/41588_2018_204_MOESM6_ESM.xlsx"
HSF="/home/raza/p53_thermal_xai/data/external/halim2026_hsf1/NIHMS2136175-supplement-4.xlsx"
CHR="17"; ORDER=["LUSC","OV","HNSC","BRCA","COADREAD","LUAD","STAD","BLCA","UCEC","GBM","LIHC","PAAD"]
SIG={"LUSC":"Lung-SCC","OV":"Ovary-AdenoCa","HNSC":"Head-SCC","BRCA":"Breast-cancer",
     "COADREAD":"ColoRect-AdenoCa","LUAD":"Lung-AdenoCa","STAD":"Stomach-AdenoCa",
     "BLCA":"Transitional-cell-carcinoma","UCEC":"Uterus-AdenoCa","GBM":"CNS-GBM",
     "LIHC":"Liver-HCC","PAAD":"Panc-AdenoCa"}
COMP={"A":"T","C":"G","G":"C","T":"A"}

uni=pd.read_csv(D+"/tp53_universe_126_292.tsv",sep="\t"); uni=uni[uni.CRISPR_RFS.notna()].reset_index(drop=True)
obs=pd.read_csv(D+"/tcga_observed_events_126_292.tsv",sep="\t"); obs=obs[obs.cancer.isin(ORDER)].reset_index(drop=True)
lo,hi=uni.gpos.min()-3,uni.gpos.max()+3
gseq=urllib.request.urlopen("https://grch37.rest.ensembl.org/sequence/region/human/%s:%d..%d:1?content-type=text/plain"%(CHR,lo,hi),timeout=60).read().decode().strip().upper()
def pent(gp,rp,ap):
    w=gseq[gp-2-lo:gp+3-lo]
    return (w,"%s>%s"%(rp,ap)) if rp in "CT" else ("".join(COMP[b] for b in reversed(w)),"%s>%s"%(COMP[rp],COMP[ap]))
uni["k1536"]=[b+"|"+a for a,b in (pent(r.gpos,r.ref_plus,r.alt_plus) for r in uni.itertuples())]
z_=zipfile.ZipFile(P+"/WES_TCGA_2018_03_09.zip")
c15=pd.read_csv(io.BytesIO(z_.read("WES_TCGA.1536.csv"))); c15["k1536"]=c15["Mutation type"]+"|"+c15["Pentanucleotide"]
samp=[c for c in c15.columns if "::" in c]; lab=pd.Series([s.split("::")[0] for s in samp],index=samp)
c15=c15.set_index("k1536")[samp]
O15=pd.read_csv(P+"/context_counts_GRCh37_1536_exome.csv",index_col=0).sum(1)
mt=c15.index.str.split("|").str[0]; pn=c15.index.str.split("|").str[1]
Ok=pd.Series(O15.reindex(pn).values,index=c15.index).astype(float)
par=pd.Series((pn.str[1]+"["+mt+"]"+pn.str[3]).values,index=c15.index)
obs_k=obs.merge(uni[["residue","wt_aa","alt_aa","k1536"]].drop_duplicates(["residue","wt_aa","alt_aa"]),
                on=["residue","wt_aa","alt_aa"],how="left")
F={}
for ct in ORDER:
    N=c15[lab[lab==SIG[ct]].index].sum(1).astype(float)
    N=N.subtract(obs_k[obs_k.cancer==ct].k1536.value_counts().reindex(N.index).fillna(0),fill_value=0).clip(lower=0)
    d=pd.DataFrame({"N":N,"O":Ok,"par":par})
    F[ct]=((d.N+0.5)/d.groupby("par").N.transform("sum").add(8.0))/(d.O/d.groupby("par").O.transform("sum"))
Fm=pd.DataFrame(F)
chan=pd.read_csv(P+"/COSMIC_v3_SBS_GRCh37_WESrenorm.txt",sep="\t",index_col=0); sigs=list(chan.columns)
Ex=pd.read_csv(P+"/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
Ex["patient"]=Ex["Sample Names"].str.slice(0,12); Ex=Ex.drop_duplicates("patient").set_index("patient")
Oex=pd.read_csv(P+"/context_counts_GRCh37_96_exome.csv",index_col=0).sum(1)
tri=chan.index.str[0]+chan.index.str[2]+chan.index.str[6]
Ovec=Oex.reindex(tri.values).values.astype(float)
ci=np.array([list(chan.index).index(c) for c in uni.sbs96])
obs=obs[obs.patient.isin(Ex.index)].reset_index(drop=True); pl=obs.patient.unique()
r_ic=(Ex.loc[pl,sigs].astype(float).values@chan.values.T/Ovec[None,:])[:,ci]
pidx={p:k for k,p in enumerate(pl)}
pct=obs.drop_duplicates("patient").set_index("patient").cancer.to_dict(); cidx={c:i for i,c in enumerate(ORDER)}
Fk=np.column_stack([Fm[c].reindex(uni.k1536).values for c in ORDER])
Wm=r_ic*np.array([Fk[:,cidx[pct[p]]] for p in pl]); Pn=Wm/Wm.sum(1,keepdims=True)
codes,uq=pd.factorize((uni.wt_aa+uni.residue.astype(str)+uni.alt_aa).values)
G=np.zeros((len(pl),len(uq))); np.add.at(G.T,codes,Pn.T); G=G/G.sum(1,keepdims=True)
rfsG=pd.Series(uni.CRISPR_RFS.values).groupby(codes).first().values
g=pd.read_excel(GIA,header=1); g["Mutation"]=g.AA_wt.astype(str)+g.Position.astype(int).astype(str)+g.AA_variant.astype(str)
c3={"D":"A549_p53WT_Nutlin-3_Z-score","N":"A549_p53NULL_Nutlin-3_Z-score","E":"A549_p53NULL_Etoposide_Z-score"}
V=g.assign(**{k:pd.to_numeric(g[c],errors="coerce") for k,c in c3.items()}).groupby("Mutation")[list(c3)].mean().reindex(uq)
h=pd.read_excel(HSF,sheet_name="mut_log2FC_summary",header=1); h=h[h.Mut_type=="Missense"]
HC="mean_nutlin3-dox-37_vs_nutlin3-37"
Hraw=pd.Series(uq).map(h.groupby("allele")[HC].mean()).values
zz=lambda a:(a-np.nanmean(a))/np.nanstd(a)
X3=np.column_stack([zz(rfsG),zz(V["D"].values),zz(V["N"].values),zz(-V["E"].values)])
X4=np.column_stack([X3,zz(Hraw)])
gmap={k:i for i,k in enumerate(uq)}
gi=np.array([gmap[k] for k in (obs.wt_aa+obs.residue.astype(str)+obs.alt_aa)])
pi=np.array([pidx[p] for p in obs.patient])

def negll(b,X):
    Z=G[pi]*np.exp(X@np.atleast_1d(b))[None,:]
    return -(np.log(Z[np.arange(len(gi)),gi])-np.log(Z.sum(1))).sum()
def info(b,X):
    """Exact observed information: sum_i [E_i(xx') - E_i(x)E_i(x)'] under the fitted multinomial."""
    Z=G[pi]*np.exp(X@b)[None,:]; Pb=Z/Z.sum(1,keepdims=True)
    m1=Pb@X                                   # events x k
    I=np.zeros((X.shape[1],X.shape[1]))
    for j in range(X.shape[1]):
        for k in range(X.shape[1]):
            I[j,k]=((Pb*(X[:,j]*X[:,k])[None,:]).sum(1)-m1[:,j]*m1[:,k]).sum()
    return I
r3=minimize(negll,[1,.8,-.1,.4],args=(X3,),method="L-BFGS-B")
r4=minimize(negll,[1,.8,-.1,.4,.1],args=(X4,),method="L-BFGS-B")
b4=r4.x; ll3,ll4=-r3.fun,-r4.fun
I=info(b4,X4); C=np.linalg.inv(I); seH=np.sqrt(C[4,4])
Happ=r4.hess_inv.todense() if hasattr(r4.hess_inv,"todense") else np.asarray(r4.hess_inv)
print("HSF1 variable entering optimizer:")
print("raw or z-standardized = z-standardized over the 990 alleles")
print("HSF1 raw SD           = %.4f" % np.nanstd(Hraw, ddof=1))
print("\nOptimizer beta_H      = %.4f" % b4[4])
print("Hessian SE_H (exact observed information) = %.4f" % seH)
print("  [previous L-BFGS-B hess_inv approximation gave %.4f - that value was wrong]" % np.sqrt(Happ[4,4]))
wz=b4[4]/seH
print("Wald z  = %.4f" % wz); print("Wald p  = %.4g" % (2*norm.sf(abs(wz))))
print("95%% CI  = %.4f - %.4f" % (b4[4]-1.96*seH, b4[4]+1.96*seH))
lr=2*(ll4-ll3)
print("\nM3 logLik = %.3f\nM4 logLik = %.3f\nLR = %.3f\nLR p = %.4g" % (ll3,ll4,lr,chi2.sf(lr,1)))
ev=np.linalg.eigvalsh(I)
print("\nCondition number of Hessian  = %.4g" % (ev.max()/ev.min()))
print("Smallest Hessian eigenvalue  = %.4g" % ev.min())
def prof(bh):
    f=lambda a: negll(np.concatenate([a,[bh]]),X4)
    return -minimize(f,b4[:4],method="L-BFGS-B").fun
crit=ll4-chi2.ppf(0.95,1)/2
lo_=brentq(lambda x: prof(x)-crit, b4[4]-1.0, b4[4]); hi_=brentq(lambda x: prof(x)-crit, b4[4], b4[4]+1.0)
print("Profile-likelihood 95%% CI for beta_H = %.4f - %.4f" % (lo_,hi_))
print("\nAll exact SEs for M4:")
for n_,v_,s_ in zip(["beta_R","beta_D","beta_N","beta_E","beta_HSF1"],b4,np.sqrt(np.diag(C))):
    print("  %-10s = %+.4f  SE %.4f  z %+.2f  95%% CI %+.4f - %+.4f" % (n_,v_,s_,v_/s_,v_-1.96*s_,v_+1.96*s_))
