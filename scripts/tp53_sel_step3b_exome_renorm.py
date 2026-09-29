#!/usr/bin/env python3
# Step 3b: genome vs WES-renormalised COSMIC v3 basis, same exposures, same catalogues.
import io, zipfile, numpy as np, pandas as pd
from scipy.stats import spearmanr

P = "/home/raza/p53_thermal_xai/data/external/pcawg_signatures"
z = zipfile.ZipFile(P + "/WES_TCGA_2018_03_09.zip")
obs = pd.read_csv(io.BytesIO(z.read("WES_TCGA.96.csv")))
lab = obs["Trinucleotide"].str[0] + "[" + obs["Mutation type"] + "]" + obs["Trinucleotide"].str[2]
samp = [c for c in obs.columns if "::" in c]
obs = obs.set_index(lab)[samp]; obs.index.name = "Type"
obs.columns = [c.split("::", 1)[1] for c in obs.columns]

ex = pd.read_csv(P + "/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
S  = pd.read_csv(P + "/COSMIC_v3_SBS_GRCh37.txt", sep="\t").set_index("Type")
sigs = [c for c in ex.columns if c.startswith("SBS")]
absent = [c for c in S.columns if c not in sigs]
print("COSMIC signatures in matrix but absent from exposures: %s" % absent)

G = pd.read_csv(P + "/context_counts_GRCh37_96.csv", index_col=0).sum(1)
E = pd.read_csv(P + "/context_counts_GRCh37_96_exome.csv", index_col=0).sum(1)
tri = obs.index.str[0] + obs.index.str[2] + obs.index.str[6]     # e.g. A[C>A]A -> ACA
ratio = (E / G).reindex(tri.values).values
print("exome/genome trinucleotide ratio: min %.4f median %.4f max %.4f"
      % (np.nanmin(ratio), np.nanmedian(ratio), np.nanmax(ratio)))

Sg = S.loc[obs.index, sigs]
Sw = Sg.mul(ratio, axis=0); Sw = Sw / Sw.sum(0)
Sw.to_csv(P + "/COSMIC_v3_SBS_GRCh37_WESrenorm.txt", sep="\t")

shared = sorted(set(ex["Sample Names"]) & set(obs.columns))
Em = ex.set_index("Sample Names").loc[shared, sigs].astype(float).values
X  = obs[shared].astype(float).values
acc = ex.set_index("Sample Names").loc[shared, "Accuracy"].astype(float).values

def cosines(Smat):
    R = Em @ Smat.values.T
    return ((X.T*R).sum(1) / (np.linalg.norm(X, axis=0)*np.linalg.norm(R, axis=1)))

cg, cw = cosines(Sg), cosines(Sw)
def block(c):
    return dict(median=np.median(c), q1=np.percentile(c,25), q3=np.percentile(c,75),
                p5=np.percentile(c,5), mn=c.min(),
                gap=np.median(np.abs(c-acc)), rho=float(spearmanr(c,acc).statistic))
g, w = block(cg), block(cw)
print("\n%-34s %16s %24s" % ("", "Genome matrix", "WES-renormalized matrix"))
for k, lbl in [("median","Median cosine"),("q1","IQR low"),("q3","IQR high"),
               ("p5","5th percentile"),("mn","Minimum"),
               ("gap","Median |cosine - reported|"),("rho","Spearman(cosine, reported)")]:
    print("%-34s %16.4f %24.4f" % (lbl, g[k], w[k]))
print("\nn samples = %d" % len(shared))
pd.DataFrame({"sample":shared,"cos_genome":cg,"cos_wes":cw,"reported":acc}).round(5).to_csv(
    "/home/raza/p53_thermal_xai/data/external/tp53_selection/reconstruction_qc_genome_vs_wes.tsv",
    sep="\t", index=False)
print("Saved: .../reconstruction_qc_genome_vs_wes.tsv and COSMIC_v3_SBS_GRCh37_WESrenorm.txt")
