#!/usr/bin/env python3
# Step 3 QC: does the supplied exposure vector reconstruct the WES catalogue it was fitted to?
# No opportunity correction here. No TP53 projection. No simulation.
import io, zipfile, numpy as np, pandas as pd
from scipy.stats import spearmanr

P = "/home/raza/p53_thermal_xai/data/external/pcawg_signatures"
EXPO = P + "/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv"
SIGM = P + "/COSMIC_v3_SBS_GRCh37.txt"
CAT  = P + "/WES_TCGA_2018_03_09.zip"

z = zipfile.ZipFile(CAT)
obs = pd.read_csv(io.BytesIO(z.read("WES_TCGA.96.csv")))
meta = [c for c in obs.columns if "::" not in c]
samp = [c for c in obs.columns if "::" in c]
print("WES_TCGA.96.csv shape:          %s" % (obs.shape,))
print("Metadata columns:               %s" % meta)
print("Number of sample columns:       %d" % len(samp))

mt = [c for c in meta if "utation" in c][0]
tri = [c for c in meta if "rinucleotide" in c][0]
lab = obs[tri].str[0] + "[" + obs[mt] + "]" + obs[tri].str[2]
print("Number of unique SBS96 channels: %d" % lab.nunique())
print("First 10 channel labels:        %s" % list(lab[:10]))
dup = lab[lab.duplicated()].tolist()
print("Duplicate channel labels:       %s" % (dup if dup else "none"))

obs = obs.set_index(lab)[samp]
obs.index.name = "Type"
obs.columns = [c.split("::", 1)[1] for c in obs.columns]

ex = pd.read_csv(EXPO)
ex_ids = ex["Sample Names"].tolist()
print("\nExposure samples:               %d" % len(ex_ids))
print("Observed SBS96 samples:         %d" % obs.shape[1])
shared = sorted(set(ex_ids) & set(obs.columns))
print("Exact matches:                  %d" % len(shared))
print("Exposure-only:                  %d" % len(set(ex_ids) - set(obs.columns)))
print("SBS96-only:                     %d" % len(set(obs.columns) - set(ex_ids)))

S = pd.read_csv(SIGM, sep="\t").set_index("Type")
sigs = [c for c in ex.columns if c.startswith("SBS")]
S = S.loc[obs.index, sigs]                                  # align channel order
E = ex.set_index("Sample Names").loc[shared, sigs].astype(float)
X = obs[shared].astype(float)
R = pd.DataFrame((E.values @ S.values.T).T, index=obs.index, columns=shared)

tot = X.sum(0)
keep = tot[tot > 0].index.tolist()
print("\nSamples shared:                 %d" % len(shared))
print("Samples with >=1 observed SBS mutation: %d" % len(keep))
print("Median observed mutation count:      %.0f" % tot[keep].median())
print("Median reconstructed mutation count: %.0f" % R[keep].sum(0).median())

a = X[keep].values; b = R[keep].values
cos = (a*b).sum(0) / (np.linalg.norm(a, axis=0) * np.linalg.norm(b, axis=0))
cos = pd.Series(cos, index=keep)
print("\nCosine similarity:")
print("  median = %.4f" % cos.median())
print("  IQR    = %.4f - %.4f" % (cos.quantile(.25), cos.quantile(.75)))
print("  5th pct= %.4f" % cos.quantile(.05))
print("  minimum= %.4f" % cos.min())

acc = ex.set_index("Sample Names").loc[keep, "Accuracy"].astype(float)
print("\nReported reconstruction accuracy:")
print("  median = %.4f" % acc.median())
r = spearmanr(cos.values, acc.values)
print("\nSpearman(your cosine, reported accuracy) = %.4f (p = %.3g, n = %d)"
      % (float(r.statistic), float(r.pvalue), len(keep)))

print("\nSignature matrix used: COSMIC_v3_SBS_GRCh37.txt")
print("  source: https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3_SBS_GRCh37.txt")
print("  96 channels x 67 signatures, genome (GRCh37) normalised, all columns sum to 1.0")
print("  %d of %d exposure signatures matched into the matrix" % (len(sigs), len(sigs)))
pd.DataFrame({"sample": keep, "cosine": cos.values, "reported_accuracy": acc.values,
              "obs_muts": tot[keep].values}).round(5).to_csv(
    "/home/raza/p53_thermal_xai/data/external/tp53_selection/reconstruction_qc.tsv", sep="\t", index=False)
print("\nSaved: /home/raza/p53_thermal_xai/data/external/tp53_selection/reconstruction_qc.tsv")
