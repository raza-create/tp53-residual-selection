#!/usr/bin/env python3
# Step 12 QC: validate HSF1 potentiation score before adding it as a predictor. No model fitting.
import numpy as np, pandas as pd
from scipy.stats import spearmanr

D   = "/home/raza/p53_thermal_xai/data/external/tp53_selection"
GIA = "/home/raza/p53_thermal_xai/data/external/giacomelli2018/41588_2018_204_MOESM6_ESM.xlsx"
HSF = "/home/raza/p53_thermal_xai/data/external/halim2026_hsf1/NIHMS2136175-supplement-4.xlsx"
HCOL = "mean_nutlin3-dox-37_vs_nutlin3-37"
REPS = ["nutlin3-dox-37_vs_nutlin3-37_rep%d" % i for i in (1,2,3)]
FOCUS = ["Y220C","V157F","R158L","I195T","R248Q","R273L","H179R","G266V",
         "R175H","R249S","G245S","R282W"]

h = pd.read_excel(HSF, sheet_name="mut_log2FC_summary", header=1)
h = h[h.Mut_type == "Missense"].copy()
for c in [HCOL] + REPS:
    h[c] = pd.to_numeric(h[c], errors="coerce")
print("HSF1 score source:\n%s\n" % HCOL)
print("DMS missense rows: %d | unique alleles: %d" % (len(h), h.allele.nunique()))

uni = pd.read_csv(D + "/tp53_universe_126_292.tsv", sep="\t")
uni = uni[uni.CRISPR_RFS.notna()]
alle = (uni.wt_aa + uni.residue.astype(str) + uni.alt_aa)
codes, uq = pd.factorize(alle.values)
rfsG = pd.Series(uni.CRISPR_RFS.values).groupby(codes).first().values

hm = h.dropna(subset=[HCOL]).groupby("allele")[HCOL].mean()
H = pd.Series(uq).map(hm).values
matched = ~np.isnan(H)
print("\nCandidate alleles matched = %d / %d" % (matched.sum(), len(uq)))
obs = pd.read_csv(D + "/tcga_observed_events_126_292.tsv", sep="\t")
ok = (obs.wt_aa + obs.residue.astype(str) + obs.alt_aa)
keepset = set(np.array(uq)[matched])
print("Observed TCGA events retained = %d / %d" % (ok.isin(keepset).sum(), len(obs)))
miss = sorted(set(uq) - set(hm.index))
print("Missing candidate alleles = %d%s" % (len(miss), (": " + ", ".join(miss[:25])) if miss else ""))

r = h.dropna(subset=REPS)
print("\nSpearman rep1 vs rep2 = %.4f" % spearmanr(r[REPS[0]], r[REPS[1]]).statistic)
print("Spearman rep1 vs rep3 = %.4f" % spearmanr(r[REPS[0]], r[REPS[2]]).statistic)
print("Spearman rep2 vs rep3 = %.4f" % spearmanr(r[REPS[1]], r[REPS[2]]).statistic)
print("  (n = %d missense alleles with all three replicates)" % len(r))

v = pd.Series(H[matched])
print("\nHSF1_potentiation (on the %d matched candidate alleles):" % matched.sum())
print("mean   = %.4f" % v.mean()); print("SD     = %.4f" % v.std(ddof=1))
print("median = %.4f" % v.median())
print("IQR    = %.4f - %.4f" % (v.quantile(.25), v.quantile(.75)))
print("range  = %.4f - %.4f" % (v.min(), v.max()))

g = pd.read_excel(GIA, header=1)
g["Mutation"] = g.AA_wt.astype(str) + g.Position.astype(int).astype(str) + g.AA_variant.astype(str)
cols3 = {"D":"A549_p53WT_Nutlin-3_Z-score","N":"A549_p53NULL_Nutlin-3_Z-score","E":"A549_p53NULL_Etoposide_Z-score"}
V = g.assign(**{k: pd.to_numeric(g[c], errors="coerce") for k,c in cols3.items()}) \
     .groupby("Mutation")[list(cols3)].mean().reindex(uq)
z = lambda a: (a-np.nanmean(a))/np.nanstd(a)
Hz = z(H)
m = matched & V.notna().all(1).values
print("\n(n = %d alleles with HSF1 and all Giacomelli scores)" % m.sum())
print("Spearman HSF1 vs CRISPR_RFS             = %+.4f" % spearmanr(Hz[m], rfsG[m]).statistic)
print("Spearman HSF1 vs WT_Nutlin_DNE          = %+.4f" % spearmanr(Hz[m], V["D"].values[m]).statistic)
print("Spearman HSF1 vs Null_Nutlin            = %+.4f" % spearmanr(Hz[m], V["N"].values[m]).statistic)
print("Spearman HSF1 vs reversed_Null_Etoposide= %+.4f" % spearmanr(Hz[m], -V["E"].values[m]).statistic)

T = pd.DataFrame({"Mutation": uq, "HSF1_potentiation": H, "HSF1_z": Hz, "CRISPR_RFS": rfsG})
for c in REPS:
    T[c] = pd.Series(uq).map(h.dropna(subset=[c]).groupby("allele")[c].mean()).values
T.round(5).to_csv(D + "/hsf1_scores_990.tsv", sep="\t", index=False)
print("\nSanity check on prespecified alleles")
print(T.set_index("Mutation").reindex(FOCUS)[["HSF1_potentiation","HSF1_z"] + REPS]
        .to_string(float_format=lambda x: "%.4f" % x))
print("\nSaved: %s/hsf1_scores_990.tsv" % D)
