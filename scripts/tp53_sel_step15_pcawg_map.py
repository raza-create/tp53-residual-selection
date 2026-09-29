import pandas as pd
from pathlib import Path

base = Path("/home/raza/p53_thermal_xai/data/external/tp53_selection")
maf_path = Path("/home/raza/p53_thermal_xai/data/external/pcawg_icgc/tp53_icgc_complete.maf")

maf = pd.read_csv(maf_path, sep="\t")
uni = pd.read_csv(base / "tp53_universe_126_292.tsv", sep="\t")

x = maf[(maf["Variant_Classification"] == "Missense_Mutation") &
        (maf["Variant_Type"] == "SNP")].copy()

def get_alt(row):
    ref = str(row["Reference_Allele"]).upper()
    a1 = str(row["Tumor_Seq_Allele1"]).upper()
    a2 = str(row["Tumor_Seq_Allele2"]).upper()
    alts = []
    for a in [a1, a2]:
        if a != ref and a not in alts:
            alts.append(a)
    return alts[0] if len(alts) == 1 else None

x["alt_derived"] = x.apply(get_alt, axis=1)
print("PCAWG missense SNP rows =", len(x))
print("Rows with unambiguous ALT =", x["alt_derived"].notna().sum())
print("Rows with ambiguous ALT =", x["alt_derived"].isna().sum())

u = uni.copy()
u["gpos"] = u["gpos"].astype(int)
u["ref_plus"] = u["ref_plus"].str.upper()
u["alt_plus"] = u["alt_plus"].str.upper()
u["Mutation"] = u["wt_aa"].astype(str) + u["residue"].astype(str) + u["alt_aa"].astype(str)
print("\nUniverse rows =", len(u))
print("Unique genomic SNVs =", u[["gpos","ref_plus","alt_plus"]].drop_duplicates().shape[0])
print("Unique protein alleles =", u["Mutation"].nunique())

m = x.merge(u, left_on=["Start_position", "Reference_Allele", "alt_derived"],
            right_on=["gpos", "ref_plus", "alt_plus"], how="left", indicator=True)
matched = m[m["_merge"] == "both"].copy()
print("\n=== EXACT 126-292 MAPPING ===")
print("Mapped PCAWG rows =", len(matched))
print("Unique tumour samples =", matched["Tumor_Sample_Barcode"].nunique())
print("Unique donors =", matched["Donor_ID"].nunique())
print("Unique protein alleles =", matched["Mutation"].nunique())
multimap = matched.groupby(["Donor_ID","Tumor_Sample_Barcode","Start_position",
                           "Reference_Allele","alt_derived"]).size()
print("Events mapping >1 universe row =", (multimap > 1).sum())

un = m[m["_merge"] == "left_only"]
inwin = un[(un.Start_position >= u.gpos.min()) & (un.Start_position <= u.gpos.max())]
print("Unmapped rows =", len(un), "| unmapped INSIDE universe span =", len(inwin), "<-- investigate if >0")
if len(inwin):
    print(inwin[["Start_position","Reference_Allele","alt_derived","Genome_Change"]].head(20).to_string(index=False))

donor_events = matched.drop_duplicates(
    subset=["Donor_ID","Start_position","Reference_Allele","alt_derived"]).copy()
print("\n=== DONOR-LEVEL DATASET ===")
print("Unique donor-variant events =", len(donor_events))
print("Unique donors =", donor_events["Donor_ID"].nunique())
nalleles = donor_events.groupby("Donor_ID")["Mutation"].nunique()
print("\nDonors with >1 distinct TP53 126-292 missense allele:")
print(nalleles[nalleles > 1].sort_values(ascending=False).to_string() if (nalleles > 1).any() else "  none")
print("\nProjects after donor-event deduplication:")
print(donor_events["Project_Code"].value_counts().to_string())

residual23 = {"R158L","V157F","H179R","R273L","C176F","G266V","Y220C","I195T","R175H","C242F",
              "R249S","R248Q","D281Y","R175C","E171K","D281E","D281V","G245D","H193R","Y163C",
              "S241F","V173L","Y234C"}
hit = donor_events[donor_events["Mutation"].isin(residual23)].copy()
print("\n=== FROZEN 23-ALLELE OVERLAP ===")
print("Residual alleles observed =", hit["Mutation"].nunique(), "/ 23")
print("Donor-level residual events =", len(hit))
print("Unique donors carrying >=1 residual allele =", hit["Donor_ID"].nunique())

counts = (donor_events.groupby("Mutation")
          .agg(PCAWG_donor_count=("Donor_ID","nunique"), N_projects=("Project_Code","nunique"))
          .reset_index())
table23 = (pd.DataFrame({"Mutation": sorted(residual23)})
           .merge(counts, on="Mutation", how="left")
           .fillna({"PCAWG_donor_count":0, "N_projects":0}))
table23["PCAWG_donor_count"] = table23["PCAWG_donor_count"].astype(int)
table23["N_projects"] = table23["N_projects"].astype(int)
print()
print(table23.sort_values(["PCAWG_donor_count","Mutation"], ascending=[False,True]).to_string(index=False))
print("\nTop 15 alleles overall in PCAWG 126-292:")
print(counts.sort_values("PCAWG_donor_count", ascending=False).head(15).to_string(index=False))

donor_events.to_csv(base / "pcawg_icgc_tp53_126_292_donor_events.tsv", sep="\t", index=False)
table23.to_csv(base / "pcawg_icgc_residual23_overlap.tsv", sep="\t", index=False)
print("\nSaved to", base)
