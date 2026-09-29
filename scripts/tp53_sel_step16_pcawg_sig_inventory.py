import pandas as pd
from pathlib import Path

base = Path("/home/raza/p53_thermal_xai/data/external/tp53_selection")
sig = pd.read_csv(base / "PCAWG_sigProfiler_SBS_signatures_in_samples.csv")
samples = pd.read_csv(base / "pcawg_sample_sheet.tsv", sep="\t")
ev = pd.read_csv(base / "pcawg_icgc_tp53_126_292_donor_events.tsv", sep="\t")

print("SIGNATURE FILE")
print("shape =", sig.shape)
print("columns =", sig.columns.tolist()[:8], "...", sig.columns.tolist()[-3:])
print(sig.iloc[:3, :8].to_string(index=False))

print("\nSAMPLE SHEET")
print("shape =", samples.shape)
print("columns =", samples.columns.tolist())

print("\nVALIDATION DATA")
print("donors =", ev["Donor_ID"].nunique())
print("tumour barcodes =", ev["Tumor_Sample_Barcode"].nunique())

idcols = [c for c in samples.columns
          if any(k in c.lower() for k in ["donor", "sample", "aliquot", "specimen", "project"])]
print("\nRelevant sample-sheet columns:")
print(idcols)
print(samples[idcols].head(10).to_string(index=False))

# how the signature file's own identifier columns relate to what we have
idc = [c for c in sig.columns if not c.startswith("SBS")]
print("\nSignature file non-SBS columns:", idc)
for c in idc:
    v = sig[c].astype(str)
    print("  %-22s n_unique=%-6d example=%s" % (c, v.nunique(), v.iloc[0]))
    print("      overlap with Tumor_Sample_Barcode: %d / %d" %
          (len(set(v) & set(ev.Tumor_Sample_Barcode.astype(str))), ev.Tumor_Sample_Barcode.nunique()))
    print("      overlap with Donor_ID:             %d / %d" %
          (len(set(v) & set(ev.Donor_ID.astype(str))), ev.Donor_ID.nunique()))
