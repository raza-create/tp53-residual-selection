#!/usr/bin/env python3
# Step 24: frozen 990-allele master table. Joins only; no analysis of ddG or RSA.
import numpy as np, pandas as pd
from pathlib import Path

B   = Path("/home/raza/p53_thermal_xai/data/external/tp53_selection")
GIA = "/home/raza/p53_thermal_xai/data/external/giacomelli2018/41588_2018_204_MOESM6_ESM.xlsx"
H   = Path("/home/raza/p53_thermal_xai/data/external/halim2026_hsf1")
A23 = {"R158L","V157F","H179R","R273L","C176F","G266V","Y220C","I195T","R175H","C242F",
       "R249S","R248Q","D281Y","R175C","E171K","D281E","D281V","G245D","H193R","Y163C",
       "S241F","V173L","Y234C"}

uni = pd.read_csv(B/"tp53_universe_126_292.tsv", sep="\t"); uni = uni[uni.CRISPR_RFS.notna()]
uni["Mutation"] = uni.wt_aa + uni.residue.astype(str) + uni.alt_aa
m = (uni.groupby("Mutation")
        .agg(residue=("residue","first"), WT_AA=("wt_aa","first"), ALT_AA=("alt_aa","first"),
             CRISPR_RFS=("CRISPR_RFS","first")).reset_index())

d = pd.read_csv(B/"residual_M4.tsv", sep="\t")
disc = pd.read_csv(B/"residual_sbs1536.tsv", sep="\t")[["Mutation","Observed","Expected_1536","OE_1536","FDR_1536"]]
m = m.merge(disc.rename(columns={"Observed":"TCGA_observed","Expected_1536":"TCGA_expected_M4_SBS1536",
                                 "OE_1536":"TCGA_OE","FDR_1536":"TCGA_FDR"}), on="Mutation", how="left")
m = m.merge(d[["Mutation","P_M4"]].rename(columns={"P_M4":"TCGA_P"}), on="Mutation", how="left")
m["Residual23"] = m.Mutation.isin(A23).astype(int)

ev = pd.read_csv(B/"pcawg_validation_events_bridged.tsv", sep="\t")
pc = ev.groupby("Mutation").agg(PCAWG_observed_donors=("Donor_ID","nunique"),
                                PCAWG_projects=("Project_Code","nunique")).reset_index()
m = m.merge(pc, on="Mutation", how="left").fillna({"PCAWG_observed_donors":0,"PCAWG_projects":0})
m[["PCAWG_observed_donors","PCAWG_projects"]] = m[["PCAWG_observed_donors","PCAWG_projects"]].astype(int)
m["PCAWG_observed_at_least_once"] = (m.PCAWG_observed_donors > 0).astype(int)

g = pd.read_excel(GIA, header=1)
g["Mutation"] = g.AA_wt.astype(str) + g.Position.astype(int).astype(str) + g.AA_variant.astype(str)
c3 = {"WT_Nutlin_DNE":"A549_p53WT_Nutlin-3_Z-score","Null_Nutlin":"A549_p53NULL_Nutlin-3_Z-score",
      "Null_Etoposide":"A549_p53NULL_Etoposide_Z-score"}
gg = g.assign(**{k: pd.to_numeric(g[v], errors="coerce") for k,v in c3.items()}).groupby("Mutation")[list(c3)].mean()
m = m.merge(gg, on="Mutation", how="left")
m["reversed_Null_Etoposide"] = -m.Null_Etoposide
m = m.drop(columns=["Null_Etoposide"])

hs = pd.read_excel(H/"NIHMS2136175-supplement-4.xlsx", sheet_name="mut_log2FC_summary", header=1)
hs = hs[hs.Mut_type == "Missense"]
m["HSF1_potentiation"] = m.Mutation.map(hs.groupby("allele")["mean_nutlin3-dox-37_vs_nutlin3-37"].mean())

dg = pd.read_excel(H/"NIHMS2136175-supplement-8.xlsx", header=1)
dg.columns = ["allele","ddG"]
m["Rosetta_ddG"] = m.Mutation.map(pd.to_numeric(dg.set_index("allele").ddG, errors="coerce"))
rs = pd.read_excel(H/"NIHMS2136175-supplement-7.xlsx", header=1)
rs.columns = ["site","RSA"]
m["RSA"] = m.residue.map(pd.to_numeric(rs.set_index("site").RSA, errors="coerce"))

cols = ["Mutation","residue","WT_AA","ALT_AA","TCGA_observed","TCGA_expected_M4_SBS1536","TCGA_OE",
        "TCGA_P","TCGA_FDR","Residual23","PCAWG_observed_donors","PCAWG_projects",
        "PCAWG_observed_at_least_once","CRISPR_RFS","WT_Nutlin_DNE","Null_Nutlin",
        "reversed_Null_Etoposide","HSF1_potentiation","Rosetta_ddG","RSA"]
m = m[cols].sort_values(["residue","ALT_AA"]).reset_index(drop=True)
m.round(6).to_csv(B/"tp53_master_table_990.tsv", sep="\t", index=False)

print("MASTER TABLE QC\n")
print("Rows = %d" % len(m))
print("Unique mutations = %d" % m.Mutation.nunique())
print("\nResidual23 = %d" % m.Residual23.sum())
print("\nRosetta ddG available = %d / 990" % m.Rosetta_ddG.notna().sum())
print("RSA available         = %d / 990" % m.RSA.notna().sum())
print("\nComplete data for:")
for lbl, c in [("CRISPR","CRISPR_RFS"),("DNE","WT_Nutlin_DNE"),("NullNutlin","Null_Nutlin"),
               ("NullEtoposide","reversed_Null_Etoposide"),("HSF1","HSF1_potentiation")]:
    print("%-14s = %d" % (lbl, m[c].notna().sum()))
print("\nPCAWG observed >=1 = %d / 990" % m.PCAWG_observed_at_least_once.sum())
print("Residual23 observed in PCAWG = %d / 23" % m[(m.Residual23==1)&(m.PCAWG_observed_at_least_once==1)].shape[0])
print("\nAny duplicate Mutation rows = %s" % ("YES" if m.Mutation.duplicated().any() else "NO"))
missing = A23 - set(m.loc[m.Residual23==1, "Mutation"])
print("Any missing Residual23 annotations = %s%s" % ("YES " + str(sorted(missing)) if missing else "NO", ""))
print("\nSaved: %s/tp53_master_table_990.tsv" % B)
