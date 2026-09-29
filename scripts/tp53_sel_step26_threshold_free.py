#!/usr/bin/env python3
# Step 26: threshold-free external validation. TCGA residual rank-percentile scored in frozen PCAWG null.
import io, zipfile, numpy as np, pandas as pd, urllib.request
from pathlib import Path
from scipy.optimize import minimize

P=Path("/home/raza/p53_thermal_xai/data/external/pcawg_signatures")
B=Path("/home/raza/p53_thermal_xai/data/external/tp53_selection")
GIA="/home/raza/p53_thermal_xai/data/external/giacomelli2018/41588_2018_204_MOESM6_ESM.xlsx"
HSF="/home/raza/p53_thermal_xai/data/external/halim2026_hsf1/NIHMS2136175-supplement-4.xlsx"
BETA=dict(R=1.0031,D=0.7109,N=-0.0912,E=0.4530,H=0.0876)
NB, SEED, CHR = 1000000, 20260929, "17"
COMP={"A":"T","C":"G","G":"C","T":"A"}
TCGA_ORD=["LUSC","OV","HNSC","BRCA","COADREAD","LUAD","STAD","BLCA","UCEC","GBM","LIHC","PAAD"]
TSIG={"LUSC":"Lung-SCC","OV":"Ovary-AdenoCa","HNSC":"Head-SCC","BRCA":"Breast-cancer",
      "COADREAD":"ColoRect-AdenoCa","LUAD":"Lung-AdenoCa","STAD":"Stomach-AdenoCa",
      "BLCA":"Transitional-cell-carcinoma","UCEC":"Uterus-AdenoCa","GBM":"CNS-GBM",
      "LIHC":"Liver-HCC","PAAD":"Panc-AdenoCa"}
rng=np.random.default_rng(SEED)

uni=pd.read_csv(B/"tp53_universe_126_292.tsv",sep="\t"); uni=uni[uni.CRISPR_RFS.notna()].reset_index(drop=True)
lo,hi=uni.gpos.min()-3,uni.gpos.max()+3
gseq=urllib.request.urlopen("https://grch37.rest.ensembl.org/sequence/region/human/%s:%d..%d:1?content-type=text/plain"%(CHR,lo,hi),timeout=60).read().decode().strip().upper()
def pent(gp,rp,ap):
    w=gseq[gp-2-lo:gp+3-lo]
    return (w,"%s>%s"%(rp,ap)) if rp in "CT" else ("".join(COMP[b] for b in reversed(w)),"%s>%s"%(COMP[rp],COMP[ap]))
uni["k1536"]=[b+"|"+a for a,b in (pent(r.gpos,r.ref_plus,r.alt_plus) for r in uni.itertuples())]
uni["Mutation"]=uni.wt_aa+uni.residue.astype(str)+uni.alt_aa
codes,uq=pd.factorize(uni.Mutation.values)
rfsG=pd.Series(uni.CRISPR_RFS.values).groupby(codes).first().values
g=pd.read_excel(GIA,header=1); g["Mutation"]=g.AA_wt.astype(str)+g.Position.astype(int).astype(str)+g.AA_variant.astype(str)
c3={"D":"A549_p53WT_Nutlin-3_Z-score","N":"A549_p53NULL_Nutlin-3_Z-score","E":"A549_p53NULL_Etoposide_Z-score"}
V=g.assign(**{k:pd.to_numeric(g[c],errors="coerce") for k,c in c3.items()}).groupby("Mutation")[list(c3)].mean().reindex(uq)
h=pd.read_excel(HSF,sheet_name="mut_log2FC_summary",header=1); h=h[h.Mut_type=="Missense"]
Hraw=pd.Series(uq).map(h.groupby("allele")["mean_nutlin3-dox-37_vs_nutlin3-37"].mean()).values
z=lambda a:(a-np.nanmean(a))/np.nanstd(a)
X=np.column_stack([z(rfsG),z(V["D"].values),z(V["N"].values),z(-V["E"].values),z(Hraw)])
gvec_allele=np.exp(X@np.array([BETA[k] for k in "RDNEH"]))

# ---------- TCGA side: cross-fitted probabilities -> E_m, V_m, Z_m ----------
chanW=pd.read_csv(P/"COSMIC_v3_SBS_GRCh37_WESrenorm.txt",sep="\t",index_col=0); sigs=list(chanW.columns)
ExT=pd.read_csv(P/"TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
ExT["patient"]=ExT["Sample Names"].str.slice(0,12); ExT=ExT.drop_duplicates("patient").set_index("patient")
OexT=pd.read_csv(P/"context_counts_GRCh37_96_exome.csv",index_col=0).sum(1)
triW=chanW.index.str[0]+chanW.index.str[2]+chanW.index.str[6]
OvecT=OexT.reindex(triW.values).values.astype(float)
ciW=np.array([list(chanW.index).index(c) for c in uni.sbs96])
zT=zipfile.ZipFile(P/"WES_TCGA_2018_03_09.zip")
c15T=pd.read_csv(io.BytesIO(zT.read("WES_TCGA.1536.csv")))
c15T["k"]=c15T["Mutation type"]+"|"+c15T["Pentanucleotide"]
sT=[c for c in c15T.columns if "::" in c]; labT=pd.Series([s.split("::")[0] for s in sT],index=sT)
c15T=c15T.set_index("k")[sT]
O15T=pd.read_csv(P/"context_counts_GRCh37_1536_exome.csv",index_col=0).sum(1)
pnT=c15T.index.str.split("|").str[1]; mtT=c15T.index.str.split("|").str[0]
OkT=pd.Series(O15T.reindex(pnT).values,index=c15T.index).astype(float)
parT=pd.Series((pnT.str[1]+"["+mtT+"]"+pnT.str[3]).values,index=c15T.index)
obsT=pd.read_csv(B/"tcga_observed_events_126_292.tsv",sep="\t"); obsT=obsT[obsT.cancer.isin(TCGA_ORD)]
obsT=obsT[obsT.patient.isin(ExT.index)].reset_index(drop=True)
obsTk=obsT.merge(uni[["residue","wt_aa","alt_aa","k1536"]].drop_duplicates(["residue","wt_aa","alt_aa"]),
                 on=["residue","wt_aa","alt_aa"],how="left")
FT={}
for ct in TCGA_ORD:
    N=c15T[labT[labT==TSIG[ct]].index].sum(1).astype(float)
    N=N.subtract(obsTk[obsTk.cancer==ct].k1536.value_counts().reindex(N.index).fillna(0),fill_value=0).clip(lower=0)
    d=pd.DataFrame({"N":N,"O":OkT,"par":parT})
    FT[ct]=((d.N+0.5)/d.groupby("par").N.transform("sum").add(8.0))/(d.O/d.groupby("par").O.transform("sum"))
FTm=pd.DataFrame(FT)
plT=obsT.patient.unique(); pidxT={p:i for i,p in enumerate(plT)}
rT=(ExT.loc[plT,sigs].astype(float).values@chanW.values.T/OvecT[None,:])[:,ciW]
pctT=obsT.drop_duplicates("patient").set_index("patient").cancer.to_dict()
FkT=np.column_stack([FTm[c].reindex(uni.k1536).values for c in TCGA_ORD])
cidxT={c:i for i,c in enumerate(TCGA_ORD)}
WT_=rT*np.array([FkT[:,cidxT[pctT[p]]] for p in plT]); PnT=WT_/WT_.sum(1,keepdims=True)
GT=np.zeros((len(plT),len(uq))); np.add.at(GT.T,codes,PnT.T); GT=GT/GT.sum(1,keepdims=True)
giT=np.array([list(uq).index(k) for k in (obsT.wt_aa+obsT.residue.astype(str)+obsT.alt_aa)])
piT=np.array([pidxT[p] for p in obsT.patient])
upT=obsT.patient.unique().copy(); rng.shuffle(upT)
foT=np.array([{p:i%5 for i,p in enumerate(upT)}[p] for p in obsT.patient])
def nllT(mask):
    def f(b):
        Z=GT[piT[mask]]*np.exp(X@b)[None,:]
        return -(np.log(Z[np.arange(mask.sum()),giT[mask]])-np.log(Z.sum(1))).sum()
    return f
PcfT=np.zeros((len(obsT),len(uq)))
for f_ in range(5):
    tr,te=foT!=f_,foT==f_
    a=minimize(nllT(tr),[1,.8,-.1,.4,.1],method="L-BFGS-B").x
    Z=GT[piT[te]]*np.exp(X@a)[None,:]; PcfT[te]=Z/Z.sum(1,keepdims=True)
Em=PcfT.sum(0); Vm=(PcfT*(1-PcfT)).sum(0); Om=np.bincount(giT,minlength=len(uq))
Zm=(Om-Em)/np.sqrt(Vm)
Rm=(pd.Series(Zm).rank(method="average").values-0.5)/len(uq)

# ---------- PCAWG side: frozen probabilities ----------
chanG=pd.read_csv(P/"COSMIC_v3_SBS_GRCh37.txt",sep="\t").set_index("Type")
exP=pd.read_csv(B/"PCAWG_sigProfiler_SBS_signatures_in_samples.csv"); exP["Sample Names"]=exP["Sample Names"].astype(str)
OgG=pd.read_csv(P/"context_counts_GRCh37_96.csv",index_col=0).sum(1)
triG=chanG.index.str[0]+chanG.index.str[2]+chanG.index.str[6]
OvecG=OgG.reindex(triG.values).values.astype(float)
ciG=np.array([list(chanG.index).index(c) for c in uni.sbs96])
Fp=pd.read_csv(B/"pcawg_F1536_by_cancer.tsv",sep="\t",index_col=0)
ev=pd.read_csv(B/"pcawg_validation_events_bridged.tsv",sep="\t"); ev["spec"]=ev.icgc_specimen_id.astype(str)
ev["ct"]=ev.spec.map(exP.set_index("Sample Names")["Cancer Types"])
don=ev.groupby("Donor_ID").agg(spec=("spec","first"),ct=("ct","first"),k=("Mutation","nunique")).reset_index()
Ep=exP.set_index("Sample Names").loc[don.spec,sigs].astype(float).values
rp=(Ep@chanG[sigs].values.T)/OvecG[None,:]
Fkp=np.column_stack([Fp[c].reindex(uni.k1536).values for c in don.ct])
wp=rp[:,ciG]*Fkp.T*gvec_allele[codes][None,:]
pp=wp/wp.sum(1,keepdims=True)
Pall=np.zeros((len(don),len(uq))); 
for j in range(len(uq)): Pall[:,j]=pp[:,codes==j].sum(1)
ev_allele=np.array([list(uq).index(k) for k in ev.Mutation])
Tobs=Rm[ev_allele].sum()

print("THRESHOLD-FREE PCAWG VALIDATION\n")
print("Alleles ranked from TCGA = %d" % len(uq))
print("PCAWG events = %d" % len(ev)); print("PCAWG donors = %d" % len(don))
print("\nTCGA residual Z:\nmin = %.3f\nmedian = %.3f\nmax = %.3f" % (Zm.min(), np.median(Zm), Zm.max()))

ks=don.k.values; cum=np.cumsum(Pall,axis=1)
def draw():
    tot=0.0
    for i in range(len(don)):
        if ks[i]==1:
            tot+=Rm[np.searchsorted(cum[i],rng.random())]
        else:
            a=np.searchsorted(cum[i],rng.random()); q=Pall[i].copy(); q[a]=0; q/=q.sum()
            b=np.searchsorted(np.cumsum(q),rng.random()); tot+=Rm[a]+Rm[b]
    return tot
Tb=np.array([draw() for _ in range(NB//100)])   # scaled; see note
print("\nObserved PCAWG rank-score T = %.3f" % Tobs)
print("Null mean T = %.3f" % Tb.mean()); print("Null SD = %.3f" % Tb.std(ddof=1))
print("Null 2.5%%-97.5%% interval = %.3f - %.3f" % (np.percentile(Tb,2.5),np.percentile(Tb,97.5)))
Zv=(Tobs-Tb.mean())/Tb.std(ddof=1)
print("Validation Z = %.3f" % Zv)
print("Empirical one-sided P = %.4g  (n_null = %d)" % ((1+(Tb>=Tobs).sum())/(len(Tb)+1), len(Tb)))
print("\nMean TCGA residual percentile among:")
print("all observed PCAWG events = %.4f" % Rm[ev_allele].mean())
print("null expectation          = %.4f" % (Tb.mean()/len(ev)))
EX={"R175H","R248Q","Y220C"}
keep=~ev.Mutation.isin(EX).values
exc=np.array([list(uq).index(m) for m in EX])
Pk=Pall.copy(); Pk[:,exc]=0; Pk=Pk/Pk.sum(1,keepdims=True); cum2=np.cumsum(Pk,axis=1)
Tobs2=Rm[ev_allele[keep]].sum()
def draw2():
    tot=0.0
    for i in range(len(don)):
        for _ in range(ks[i]): tot+=Rm[np.searchsorted(cum2[i],rng.random())]
    return tot
Tb2=np.array([draw2() for _ in range(NB//100)])
print("\nSensitivity excluding R175H, R248Q, Y220C:")
print("Observed T = %.3f" % Tobs2); print("Null mean = %.3f" % Tb2.mean())
print("Validation Z = %.3f" % ((Tobs2-Tb2.mean())/Tb2.std(ddof=1)))
print("P = %.4g" % ((1+(Tb2>=Tobs2).sum())/(len(Tb2)+1)))
pd.DataFrame({"Mutation":uq,"O":Om,"E":Em,"V":Vm,"Z":Zm,"R":Rm}).round(6).to_csv(B/"tcga_residual_ranks_990.tsv",sep="\t",index=False)
print("\nSaved: tcga_residual_ranks_990.tsv")
