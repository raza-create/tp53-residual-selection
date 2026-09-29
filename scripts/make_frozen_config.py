#!/usr/bin/env python3
import numpy as np, pandas as pd, yaml
from pathlib import Path
B=Path("/home/raza/p53_thermal_xai/data/external/tp53_selection")
GIA="/home/raza/p53_thermal_xai/data/external/giacomelli2018/41588_2018_204_MOESM6_ESM.xlsx"
HSF="/home/raza/p53_thermal_xai/data/external/halim2026_hsf1/NIHMS2136175-supplement-4.xlsx"
m=pd.read_csv(B/"tp53_master_table_990.tsv",sep="\t")
cols={"R":"CRISPR_RFS","D":"WT_Nutlin_DNE","N":"Null_Nutlin","E":"reversed_Null_Etoposide","H":"HSF1_potentiation"}
std={k:{"mean":float(np.nanmean(m[c])),"sd_ddof0":float(np.nanstd(m[c]))} for k,c in cols.items()}
cfg=dict(
  genome_build="GRCh37/hg19", transcript="ENST00000269305 / NM_000546",
  residue_interval=[126,292], n_nucleotide_candidates=1102, n_protein_alleles=990,
  tcga_events=1819, tcga_cancer_types=12, pcawg_donors=359, pcawg_events=367,
  pcawg_residual23_carriers=88, seed=20260929, n_folds=5,
  m4_coefficients=dict(R=1.0031,D=0.7109,N=-0.0912,E=0.4530,H=0.0876),
  m4_exact_se=dict(R=0.0752,D=0.0437,N=0.0535,E=0.0519,H=0.0262),
  m4_loglik=-9793.690, m3_loglik=-9799.317,
  standardization_tcga_990=std,
  predictor_orientation={
    "R":"CRISPR rfs_crispr; higher = greater loss of function",
    "D":"Giacomelli A549_p53WT_Nutlin-3_Z-score; higher = stronger dominant-negative",
    "N":"Giacomelli A549_p53NULL_Nutlin-3_Z-score; as published",
    "E":"NEGATED Giacomelli A549_p53NULL_Etoposide_Z-score; higher = greater disruption",
    "H":"Halim mean_nutlin3-dox-37_vs_nutlin3-37; higher = greater HSF1 potentiation"},
  signature_basis=dict(tcga="COSMIC_v3_SBS_GRCh37 WES-renormalised; exome 96/1536 opportunities",
                       pcawg="COSMIC_v3_SBS_GRCh37 genome; genome 96/1536 opportunities"),
  residual23=sorted(m.loc[m.Residual23==1,"Mutation"].tolist()))
Path("/home/raza/p53_thermal_xai/config").mkdir(exist_ok=True)
yaml.safe_dump(cfg,open("/home/raza/p53_thermal_xai/config/frozen_model.yaml","w"),sort_keys=False)
print(yaml.safe_dump(cfg,sort_keys=False))
assert len(cfg["residual23"])==23
