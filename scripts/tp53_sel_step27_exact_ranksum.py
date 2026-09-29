#!/usr/bin/env python3
# Step 27: exact null distribution of the summed rank score by donor-PMF convolution.
import io, zipfile, numpy as np, pandas as pd, urllib.request
from pathlib import Path
from scipy.stats import norm
import sys as _s; _s.path.insert(0,'/home/raza/p53_thermal_xai/lib')

P=Path("/home/raza/p53_thermal_xai/data/external/pcawg_signatures")
B=Path("/home/raza/p53_thermal_xai/data/external/tp53_selection")
GIA="/home/raza/p53_thermal_xai/data/external/giacomelli2018/41588_2018_204_MOESM6_ESM.xlsx"
HSF="/home/raza/p53_thermal_xai/data/external/halim2026_hsf1/NIHMS2136175-supplement-4.xlsx"
BETA=dict(R=1.0031,D=0.7109,N=-0.0912,E=0.4530,H=0.0876)
CHR="17"; COMP={"A":"T","C":"G","G":"C","T":"A"}

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
gvec=np.exp(X@np.array([BETA[k] for k in "RDNEH"]))

rk=pd.read_csv(B/"tcga_residual_ranks_990.tsv",sep="\t").set_index("Mutation")
Rm=rk.R.reindex(uq).values
IR=np.rint(Rm*990+0.5).astype(int)          # integer ranks 1..990

chanG=pd.read_csv(P/"COSMIC_v3_SBS_GRCh37.txt",sep="\t").set_index("Type")
exP=pd.read_csv(B/"PCAWG_sigProfiler_SBS_signatures_in_samples.csv"); exP["Sample Names"]=exP["Sample Names"].astype(str)
sigs=[c for c in exP.columns if c.startswith("SBS")]
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
wp=rp[:,ciG]*Fkp.T*gvec[codes][None,:]
pp=wp/wp.sum(1,keepdims=True)
Pall=np.zeros((len(don),len(uq)))
for j in range(len(uq)): Pall[:,j]=pp[:,codes==j].sum(1)
Pall=Pall/Pall.sum(1,keepdims=True)
amap={m:i for i,m in enumerate(uq)}
ev_a=ev.Mutation.map(amap).values
ks=don.k.values

def donor_pmf(pv,k,maxr):
    if k==1:
        out=np.zeros(maxr+1); np.add.at(out,IR,pv); return out
    idx=np.argsort(-pv)[:400]; q=pv[idx]; q=q/pv.sum()
    out=np.zeros(2*maxr+1)
    denom=1-q
    for a in range(len(idx)):
        w=q[a]*q/np.clip(1-q[a],1e-300,None)
        w[a]=0
        np.add.at(out, IR[idx[a]]+IR[idx], w)
    s=out.sum()
    return out/s if s>0 else out
def convolve(pmfs):
    while len(pmfs)>1:
        nxt=[]
        for i in range(0,len(pmfs),2):
            nxt.append(np.convolve(pmfs[i],pmfs[i+1]) if i+1<len(pmfs) else pmfs[i])
        pmfs=nxt
    return pmfs[0]
def exact(mask_alleles=None):
    Pm=Pall.copy()
    if mask_alleles is not None:
        Pm[:,mask_alleles]=0; Pm=Pm/Pm.sum(1,keepdims=True)
    pmfs=[donor_pmf(Pm[i],ks[i],990) for i in range(len(don))]
    pmf=convolve([np.asarray(q,dtype=np.longdouble) for q in pmfs])
    pmf=np.clip(pmf,0,None); pmf/=pmf.sum()
    x=np.arange(len(pmf)); mu=float((x*pmf).sum()); sd=float(np.sqrt(((x-mu)**2*pmf).sum()))
    return pmf,x,mu,sd

print("EXACT THRESHOLD-FREE PCAWG VALIDATION\n")
pmf,x,mu,sd=exact()
Tobs=int(IR[ev_a].sum())
cdf=np.cumsum(pmf)
print("Observed rank-sum  = %d" % Tobs)
print("Exact null mean    = %.3f" % mu)
print("Exact null SD      = %.3f" % sd)
print("Exact 95%% interval = %d - %d" % (np.searchsorted(cdf,0.025), np.searchsorted(cdf,0.975)))
p1=pmf[Tobs:].sum()
print("Exact upper-tail P = %.4g" % p1)
print("Equivalent validation Z = %.3f" % ((Tobs-mu)/sd))
print("  (normal-tail check on Z: P = %.4g)" % norm.sf(float((Tobs-mu)/sd)))

EX=[amap[m] for m in ["R175H","R248Q","Y220C"]]
keep=~ev.Mutation.isin(["R175H","R248Q","Y220C"]).values
pmf2,x2,mu2,sd2=exact(np.array(EX))
Tobs2=int(IR[ev_a[keep]].sum())
print("\nTop-3-removed:")
print("Observed rank-sum  = %d" % Tobs2)
print("Exact null mean    = %.3f" % mu2)
print("Exact null SD      = %.3f" % sd2)
print("Exact upper-tail P = %.4g" % pmf2[Tobs2:].sum())
print("Equivalent Z       = %.3f" % ((Tobs2-mu2)/sd2))
print("  (normal-tail check on Z: P = %.4g)" % norm.sf(float((Tobs2-mu2)/sd2)))
