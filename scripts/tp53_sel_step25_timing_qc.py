#!/usr/bin/env python3
# Step 25: map frozen PCAWG TP53 events to the panorama driver table. QC only, no statistics.
import gzip, numpy as np, pandas as pd
from pathlib import Path

B = Path("/home/raza/p53_thermal_xai/data/external/tp53_selection")
Q = Path("/home/raza/p53_thermal_xai/data/external/pcawg_icgc")
A23 = {"R158L","V157F","H179R","R273L","C176F","G266V","Y220C","I195T","R175H","C242F",
       "R249S","R248Q","D281Y","R175C","E171K","D281E","D281V","G245D","H193R","Y163C",
       "S241F","V173L","Y234C"}

t = pd.read_csv(Q/"TableS3_panorama_driver_mutations_ICGC_samples.public.tsv.gz", sep="\t", compression="gzip")
print("TIMING MATCH QC\n")
print("Timing-table rows        =", len(t))
tp = t[t.gene == "TP53"].copy()
print("TP53 rows in timing table =", len(tp))
tp["key"] = (tp.sample_id.astype(str) + "|" + tp.chr.astype(str) + "|" +
             tp.pos.astype(str) + "|" + tp.ref.astype(str) + "|" + tp.alt.astype(str))
print("\nUnique timing mutation keys    =", tp.key.nunique())
print("Duplicate timing mutation keys =", int(tp.key.duplicated().sum()))

ev = pd.read_csv(B/"pcawg_validation_events_bridged.tsv", sep="\t")
ev["spec"] = ev.icgc_specimen_id.astype(str)
alt = "alt_derived" if "alt_derived" in ev.columns else "ALT"
ev["key"] = (ev.spec + "|" + ev.Chromosome.astype(str) + "|" + ev.Start_position.astype(str) + "|" +
             ev.Reference_Allele.astype(str) + "|" + ev[alt].astype(str))
ev["is23"] = ev.Mutation.isin(A23)
print("\nsample_id examples in driver table:", tp.sample_id.astype(str).head(3).tolist())
print("specimen examples in validation set:", ev.spec.head(3).tolist())
print("overlap of sample_id with icgc_specimen_id: %d / %d" %
      (len(set(tp.sample_id.astype(str)) & set(ev.spec)), ev.spec.nunique()))
print("overlap of sample_id with Tumor_Sample_Barcode: %d" %
      len(set(tp.sample_id.astype(str)) & set(ev.Tumor_Sample_Barcode.astype(str))))

m = ev.merge(tp.drop_duplicates("key"), on="key", how="left", indicator="_match")
ok = m._match == "both"
print("\nValidation TP53 events =", len(ev))
print("Matched exactly        =", int(ok.sum()))
print("Unmatched              =", int((~ok).sum()))
print("Matched %%              = %.1f" % (100*ok.mean()))
print("\nResidual23 events      =", int(ev.is23.sum()))
print("Residual23 matched     =", int((ok & m.is23).sum()))
print("Non-residual events    =", int((~ev.is23).sum()))
print("Non-residual matched   =", int((ok & ~m.is23).sum()))
print("\nUnique validation donors represented =", m.loc[ok, "Donor_ID"].nunique())

for c in ["category", "top_category", "driver", "driver_statement", "biallelic"]:
    if c in m.columns:
        print("\n%s values:" % c)
        print(m.loc[ok, c].value_counts(dropna=False).to_string())

cols = [c for c in ["Donor_ID","sample_id","Mutation","chr","pos","ref","alt",
                    "category","top_category","driver","driver_statement","biallelic"] if c in m.columns]
print("\nFirst 10 matched rows:")
print(m.loc[ok, cols].head(10).to_string(index=False))
m.loc[ok].to_csv(B/"pcawg_tp53_driver_annotated.tsv", sep="\t", index=False)
print("\nSaved: pcawg_tp53_driver_annotated.tsv")
