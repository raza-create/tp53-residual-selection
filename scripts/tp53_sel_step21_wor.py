#!/usr/bin/env python3
# Step 21: without-replacement sensitivity for the 8 two-event donors. Nothing refitted.
import numpy as np, pandas as pd, urllib.request
from pathlib import Path

P = Path("/home/raza/p53_thermal_xai/data/external/pcawg_signatures")
B = Path("/home/raza/p53_thermal_xai/data/external/tp53_selection")
GIA = "/home/raza/p53_thermal_xai/data/external/giacomelli2018/41588_2018_204_MOESM6_ESM.xlsx"
HSF = "/home/raza/p53_thermal_xai/data/external/halim2026_hsf1/NIHMS2136175-supplement-4.xlsx"
BETA = dict(R=1.0031, D=0.7109, N=-0.0912, E=0.4530, H=0.0876)
A23 = {"R158L","V157F","H179R","R273L","C176F","G266V","Y220C","I195T","R175H","C242F",
       "R249S","R248Q","D281Y","R175C","E171K","D281E","D281V","G245D","H193R","Y163C",
       "S241F","V173L","Y234C"}
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
gm=pd.Series(np.exp(BETA["R"]*z(rfsG)+BETA["D"]*z(V["D"].values)+BETA["N"]*z(V["N"].values)
                    +BETA["E"]*z(-V["E"].values)+BETA["H"]*z(Hraw)),index=uq)
gvec=gm.reindex(uni.Mutation).values
chan=pd.read_csv(P/"COSMIC_v3_SBS_GRCh37.txt",sep="\t").set_index("Type")
ex=pd.read_csv(B/"PCAWG_sigProfiler_SBS_signatures_in_samples.csv"); ex["Sample Names"]=ex["Sample Names"].astype(str)
sigs=[c for c in ex.columns if c.startswith("SBS")]
Og=pd.read_csv(P/"context_counts_GRCh37_96.csv",index_col=0).sum(1)
tri=chan.index.str[0]+chan.index.str[2]+chan.index.str[6]
Ovec=Og.reindex(tri.values).values.astype(float)
ci=np.array([list(chan.index).index(c) for c in uni.sbs96])
Fm=pd.read_csv(B/"pcawg_F1536_by_cancer.tsv",sep="\t",index_col=0)
ev=pd.read_csv(B/"pcawg_validation_events_bridged.tsv",sep="\t"); ev["spec"]=ev.icgc_specimen_id.astype(str)
ev["ct"]=ev.spec.map(ex.set_index("Sample Names")["Cancer Types"])
don=ev.groupby("Donor_ID").agg(spec=("spec","first"),ct=("ct","first"),k=("Mutation","nunique"),
                               carrier=("Mutation",lambda s:int(len(set(s)&A23)>0))).reset_index()
E=ex.set_index("Sample Names").loc[don.spec,sigs].astype(float).values
r=(E@chan[sigs].values.T)/Ovec[None,:]
Fk=np.column_stack([Fm[c].reindex(uni.k1536).values for c in don.ct])
w=r[:,ci]*Fk.T*gvec[None,:]; p=w/w.sum(1,keepdims=True)
inA=uni.Mutation.isin(A23).values
s=p[:,inA].sum(1)
pi_wr=1-(1-s)**don.k.values
pi_wor=s.copy()
two=np.where(don.k.values==2)[0]
for i in two:
    pm=p[i]; si=s[i]; nm=pm[~inA]
    P0=np.sum(nm*(1-si-nm)/(1-nm))
    pi_wor[i]=1-P0
def pb(pv):
    N=len(pv)+1; wr=np.exp(2j*np.pi*np.arange(N)/N)
    pmf=np.real(np.fft.fft(np.prod(1-pv[None,:]+pv[None,:]*wr[:,None],axis=1))/N)
    return np.clip(pmf,0,None)/np.clip(pmf,0,None).sum()
pmf=pb(pi_wor); cdf=np.cumsum(pmf); x=np.arange(len(pmf)); obs=int(don.carrier.sum())
mu=(x*pmf).sum(); sd=np.sqrt(((x-mu)**2*pmf).sum())
print("WITHOUT-REPLACEMENT SENSITIVITY\n")
print("Observed carriers = %d\n" % obs)
print("Expected carriers:")
print("with replacement    = %.3f" % pi_wr.sum())
print("without replacement = %.3f" % pi_wor.sum())
print("\nO/E without replacement = %.3f" % (obs/pi_wor.sum()))
print("\nPoisson-binomial:")
print("mean = %.3f" % mu); print("SD   = %.3f" % sd)
print("95%% interval = %d - %d" % (np.searchsorted(cdf,0.025), np.searchsorted(cdf,0.975)))
print("P(X >= %d) = %.4g" % (obs, pmf[obs:].sum()))
d=np.abs(pi_wr-pi_wor)
print("\nMaximum absolute change in pi among %d two-event donors = %.3g" % (len(two), d[two].max()))
print("  (their pi_wr -> pi_wor: %s)" % ", ".join("%.4f->%.4f"%(pi_wr[i],pi_wor[i]) for i in two))
