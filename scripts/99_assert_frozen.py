#!/usr/bin/env python3
import sys, yaml, numpy as np, pandas as pd
from pathlib import Path
B=Path("/home/raza/p53_thermal_xai/data/external/tp53_selection")
c=yaml.safe_load(open("/home/raza/p53_thermal_xai/config/frozen_model.yaml"))
R=c["reported_results"]; A=set(c["residual23"]); ok=[]
def chk(n,got,exp):
    r=got==exp; print("%-38s %12s  expected %-12s %s"%(n,got,exp,"OK" if r else "FAIL")); ok.append(r)
def close(n,got,exp,rtol):
    r=abs(got-exp)<=rtol*abs(exp)
    print("%-38s %12.6g  expected %-12.6g %s"%(n,got,exp,"OK" if r else "FAIL")); ok.append(r)
uni=pd.read_csv(B/"tp53_universe_126_292.tsv",sep="\t"); uni=uni[uni.CRISPR_RFS.notna()]
m=pd.read_csv(B/"tp53_master_table_990.tsv",sep="\t")
ev=pd.read_csv(B/"pcawg_validation_events_bridged.tsv",sep="\t")
ob=pd.read_csv(B/"tcga_observed_events_126_292.tsv",sep="\t")
T12=['LUSC','OV','HNSC','BRCA','COADREAD','LUAD','STAD','BLCA','UCEC','GBM','LIHC','PAAD']
chk("nucleotide candidates",len(uni),c["n_nucleotide_candidates"])
chk("protein alleles",len(m),c["n_protein_alleles"])
chk("TCGA events",len(ob[ob.cancer.isin(T12)]),c["tcga_events"])
chk("Residual23 size",int(m.Residual23.sum()),23)
chk("PCAWG donors",ev.Donor_ID.nunique(),c["pcawg_donors"])
chk("PCAWG events",len(ev),c["pcawg_events"])
chk("PCAWG residual23 carriers",ev[ev.Mutation.isin(A)].Donor_ID.nunique(),c["pcawg_residual23_carriers"])
d=pd.read_csv(B/"pcawg_validation_donor_probs_WOR.tsv",sep="\t")
close("PCAWG expected carriers (WOR)",d.pi.sum(),R["pcawg_expected_carriers"],1e-6)
close("PCAWG O/E (WOR)",d.carrier.sum()/d.pi.sum(),R["pcawg_OE"],1e-6)
sys.path.insert(0,"/home/raza/p53_thermal_xai/lib")
from pbtail import poisson_binomial_upper_tail
close("PCAWG exact carrier P (WOR)",poisson_binomial_upper_tail(d.pi.to_numpy(dtype=float),int(d.carrier.sum())),R["pcawg_exact_carrier_P"],1e-10)
rk=pd.read_csv(B/"tcga_residual_ranks_990.tsv",sep="\t").set_index("Mutation")
IR=np.rint(rk.R.reindex(ev.Mutation).values*990+0.5).astype(int)
chk("threshold-free observed rank-sum",int(IR.sum()),R["threshold_free_observed_ranksum"])
close("threshold-free exact P",R["threshold_free_exact_P"],2.743e-24,1e-3)
close("top-3-removed exact P",R["top3_removed_exact_P"],1.631e-09,1e-3)
print("\nALL PASS" if all(ok) else "\nFAILURES PRESENT")
sys.exit(0 if all(ok) else 1)
