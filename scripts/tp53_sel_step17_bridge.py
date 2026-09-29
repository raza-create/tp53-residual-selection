#!/usr/bin/env python3
# Step 17: bridge Tumor_Sample_Barcode -> aliquot_id -> icgc_specimen_id -> signature Sample Names. Mapping QC only.
import pandas as pd
from pathlib import Path

base = Path("/home/raza/p53_thermal_xai/data/external/tp53_selection")
ev  = pd.read_csv(base / "pcawg_icgc_tp53_126_292_donor_events.tsv", sep="\t")
ss  = pd.read_csv(base / "pcawg_sample_sheet.tsv", sep="\t")
sig = pd.read_csv(base / "PCAWG_sigProfiler_SBS_signatures_in_samples.csv")

print("validation events         =", len(ev))
print("validation donors         =", ev.Donor_ID.nunique())
print("validation tumour aliquots=", ev.Tumor_Sample_Barcode.nunique())

al = set(ev.Tumor_Sample_Barcode.astype(str))
ss["aliquot_id"] = ss.aliquot_id.astype(str)
ss_sub = ss[ss.aliquot_id.isin(al)].copy()
print("\n=== ALIQUOT -> SAMPLE-SHEET MATCH ===")
print("sample-sheet rows matched          =", len(ss_sub))
print("unique validation aliquots matched =", ss_sub.aliquot_id.nunique())
print("unique ICGC specimens              =", ss_sub.icgc_specimen_id.nunique())
print("unique ICGC donors                 =", ss_sub.icgc_donor_id.nunique())
print("validation aliquots missing from sample sheet =", len(al - set(ss_sub.aliquot_id)))
print("\nSpecimen types:")
print(ss_sub.dcc_specimen_type.value_counts().to_string())
per = ss_sub.groupby("aliquot_id").icgc_specimen_id.nunique()
print("\nAliquots mapping to >1 specimen =", int((per > 1).sum()))

bridge = ss_sub[["aliquot_id","icgc_specimen_id","icgc_donor_id","dcc_project_code"]].drop_duplicates("aliquot_id")
ev2 = ev.merge(bridge, left_on="Tumor_Sample_Barcode", right_on="aliquot_id", how="left")
print("\nevents with specimen ID =", ev2.icgc_specimen_id.notna().sum(), "/", len(ev2))
print("donors with specimen ID =", ev2[ev2.icgc_specimen_id.notna()].Donor_ID.nunique(), "/", ev.Donor_ID.nunique())

sig["Sample Names"] = sig["Sample Names"].astype(str)
sigset = set(sig["Sample Names"])
ev2["has_sbs"] = ev2.icgc_specimen_id.astype(str).isin(sigset)
print("\nevents with SBS exposure    =", int(ev2.has_sbs.sum()))
print("events without SBS exposure =", int((~ev2.has_sbs).sum()))
print("donors with SBS exposure    =", ev2[ev2.has_sbs].Donor_ID.nunique(), "/", ev.Donor_ID.nunique())
print("tumour specimens with SBS exposure =", ev2[ev2.has_sbs].icgc_specimen_id.nunique())

mm = ev2.dropna(subset=["icgc_donor_id"])
print("\nDonor-ID mismatches =", int((mm.Donor_ID.astype(str) != mm.icgc_donor_id.astype(str)).sum()))

print("\nMAF Project_Code -> signature Cancer Types")
j = ev2[ev2.has_sbs].merge(sig[["Sample Names","Cancer Types"]],
                           left_on=ev2[ev2.has_sbs].icgc_specimen_id.astype(str),
                           right_on="Sample Names", how="left")
print(pd.crosstab(j.Project_Code, j["Cancer Types"]).replace(0, "").to_string())

if (~ev2.has_sbs).any():
    print("\nEvents lacking SBS, by project:")
    print(ev2[~ev2.has_sbs].Project_Code.value_counts().to_string())
ev2.to_csv(base / "pcawg_validation_events_bridged.tsv", sep="\t", index=False)
print("\nSaved: pcawg_validation_events_bridged.tsv")
