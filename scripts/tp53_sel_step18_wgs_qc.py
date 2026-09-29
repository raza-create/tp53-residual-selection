#!/usr/bin/env python3
# Step 18: PCAWG WGS reconstruction QC using the GENOME-normalised COSMIC v3 matrix.
import io, zipfile, numpy as np, pandas as pd
from scipy.stats import spearmanr
from pathlib import Path

P = Path("/home/raza/p53_thermal_xai/data/external/pcawg_signatures")
B = Path("/home/raza/p53_thermal_xai/data/external/tp53_selection")

z = zipfile.ZipFile(P / "WGS_PCAWG_2018_02_09.zip")
name = [n for n in z.namelist() if not n.startswith("__") and "96" in n][0]
obs = pd.read_csv(io.BytesIO(z.read(name)))
print("PCAWG WGS RECONSTRUCTION QC\n")
print("catalogue file used             =", name)
meta = [c for c in obs.columns if "::" not in c and obs[c].dtype == object][:2]
mt  = [c for c in obs.columns if "utation" in c][0]
tri = [c for c in obs.columns if "rinucleotide" in c][0]
lab = obs[tri].str[0] + "[" + obs[mt] + "]" + obs[tri].str[2]
samp = [c for c in obs.columns if c not in (mt, tri)]
obs = obs.set_index(lab)[samp]; obs.index.name = "Type"
obs.columns = [c.split("::", 1)[-1] for c in obs.columns]
print("Observed SBS96 catalogue shape  =", obs.shape)
print("Observed SBS96 samples          =", obs.shape[1])

ex = pd.read_csv(B / "PCAWG_sigProfiler_SBS_signatures_in_samples.csv")
ex["Sample Names"] = ex["Sample Names"].astype(str)
print("Exposure samples                =", len(ex))
shared = sorted(set(ex["Sample Names"]) & set(obs.columns))
print("Exact sample-name matches       =", len(shared))

S = pd.read_csv(P / "COSMIC_v3_SBS_GRCh37.txt", sep="\t").set_index("Type")
sigs = [c for c in ex.columns if c.startswith("SBS")]
print("\nSignatures in exposure matrix          =", len(sigs))
print("Signatures matched to COSMIC genome mx =", len([c for c in sigs if c in S.columns]))
print("Unmatched exposure signatures          =", [c for c in sigs if c not in S.columns] or "none")
print("COSMIC signatures absent from exposures=", [c for c in S.columns if c not in sigs] or "none")
print("full signature list:", sigs)

S = S.loc[obs.index, sigs]
E = ex.set_index("Sample Names").loc[shared, sigs].astype(float).values
X = obs[shared].astype(float).values
R = (E @ S.values.T).T
acc = ex.set_index("Sample Names").loc[shared, "Accuracy"].astype(float).values
cos = (X*R).sum(0)/(np.linalg.norm(X, axis=0)*np.linalg.norm(R, axis=0))

def block(title, idx):
    c = cos[idx]; a = acc[idx]
    print("\n%s" % title)
    print("n = %d" % len(c))
    print("Median observed mutation count      = %.0f" % np.median(X[:, idx].sum(0)))
    print("Median reconstructed mutation count = %.0f" % np.median(R[:, idx].sum(0)))
    print("Cosine similarity:")
    print("  median = %.4f" % np.median(c))
    print("  IQR    = %.4f - %.4f" % (np.percentile(c,25), np.percentile(c,75)))
    print("  5th pct= %.4f" % np.percentile(c,5))
    print("  minimum= %.4f" % c.min())
    print("Reported Accuracy: median = %.4f" % np.median(a))
    print("Median |own cosine - reported Accuracy| = %.4f" % np.median(np.abs(c-a)))
    print("Spearman(own cosine, reported Accuracy) = %.4f" % float(spearmanr(c, a).statistic))

block("ALL %d PCAWG WGS" % len(shared), np.arange(len(shared)))

ev = pd.read_csv(B / "pcawg_validation_events_bridged.tsv", sep="\t")
vs = set(ev.icgc_specimen_id.astype(str))
idx = np.array([i for i, s in enumerate(shared) if s in vs])
print("\nvalidation specimens found in QC set = %d / %d" % (len(idx), len(vs)))
block("%d VALIDATION SPECIMENS" % len(idx), idx)

pd.DataFrame({"sample": shared, "cosine": cos, "reported_accuracy": acc}).round(5).to_csv(
    B / "pcawg_wgs_reconstruction_qc.tsv", sep="\t", index=False)
print("\nSaved: pcawg_wgs_reconstruction_qc.tsv")
