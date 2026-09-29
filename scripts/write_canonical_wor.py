#!/usr/bin/env python3
import numpy as np, pandas as pd, urllib.request, yaml
from pathlib import Path
P=Path("/home/raza/p53_thermal_xai/data/external/pcawg_signatures")
B=Path("/home/raza/p53_thermal_xai/data/external/tp53_selection")
C=Path("/home/raza/p53_thermal_xai/config")
GIA="/home/raza/p53_thermal_xai/data/external/giacomelli2018/41588_2018_204_MOESM6_ESM.xlsx"
HSF="/home/raza/p53_thermal_xai/data/external/halim2026_hsf1/NIHMS2136175-supplement-4.xlsx"
cfg=yaml.safe_load(open(C/"frozen_model.yaml"))
BETA=cfg["m4_coefficients_validation_frozen"]; A23=set(cfg["residual23"])
CHR="17"; COMP={"A":"T","C":"G","G":"C","T":"A"}
uni=pd.read_csv(B/"tp53_universe_126_292.tsv",sep="\t"); uni=uni[uni.CRISPR_RFS.notna()].reset_index(drop=True)
lo,hi=uni.gpos.min()-3,uni.gpos.max()+3
gs=urllib.request.urlopen("https://grch37.rest.ensembl.org/sequence/region/human/%s:%d..%d:1?content-type=text/plain"%(CHR,lo,hi),timeout=60).read().decode().strip().upper()
def pent(gp,rp,ap):
    w=gs[gp-2-lo:gp+3-lo]
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
gvec=np.exp(BETA["R"]*z(rfsG)+BETA["D"]*z(V["D"].values)+BETA["N"]*z(V["N"].values)
            +BETA["E"]*z(-V["E"].values)+BETA["H"]*z(Hraw))
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
w=r[:,ci]*Fk.T*gvec[codes][None,:]; p=w/w.sum(1,keepdims=True)
inA=uni.Mutation.isin(A23).values
s=p[:,inA].sum(1); ks=don.k.values
pi_wr=1-(1-s)**ks
pi=s.copy()
for i in np.where(ks==2)[0]:
    nm=p[i][~inA]; pi[i]=1-np.sum(nm*(1-s[i]-nm)/(1-nm))
import sys as _s; _s.path.insert(0,"/home/raza/p53_thermal_xai/lib")
from pbtail import poisson_binomial_upper_tail
obs=int(don.carrier.sum())
don.assign(s=s,pi_wr=pi_wr,pi=pi).to_csv(B/"pcawg_validation_donor_probs_WOR.tsv",sep="\t",index=False,float_format="%.17g")
pi=pd.read_csv(B/"pcawg_validation_donor_probs_WOR.tsv",sep="\t").pi.to_numpy(dtype=float)
out=dict(observed_carriers=obs, expected_carriers_wor=float(pi.sum()),
         expected_carriers_wr=float(pi_wr.sum()), OE_wor=float(obs/pi.sum()),
         exact_P_wor=poisson_binomial_upper_tail(pi,obs))
don.assign(s=s,pi_wr=pi_wr,pi=pi).to_csv(B/"pcawg_validation_donor_probs_WOR.tsv",sep="\t",index=False,float_format="%.17g")
cfg["reported_results"].update(pcawg_expected_carriers=out["expected_carriers_wor"],
    pcawg_OE=out["OE_wor"], pcawg_exact_carrier_P=out["exact_P_wor"],
    pcawg_expected_carriers_with_replacement=out["expected_carriers_wr"])
cfg["canonical_validation_file"]="pcawg_validation_donor_probs_WOR.tsv (without-replacement; primary)"
yaml.safe_dump(cfg,open(C/"frozen_model.yaml","w"),sort_keys=False)
for k,v in out.items(): print("%-32s %.10g"%(k,v))
