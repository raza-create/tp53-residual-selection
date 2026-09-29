#!/usr/bin/env python3
# Replica fragility of mutant-vs-WT conclusions. 20 mutants x 8 metrics = 160 conclusions.
# ROBUST only if the two mutant replicas lie entirely on one side of both WT replicas.
import sys
from pathlib import Path
import pandas as pd

src = Path("/home/raza/p53_thermal_xai/analysis/replica_icc/trajectory_summaries_20_100ns.tsv")
outdir = Path("/home/raza/p53_thermal_xai/analysis/replica_icc")
metrics = ["RMSD", "Q", "beta_core", "L1", "L2_Zn", "L3", "H2", "other_DBD"]

df = pd.read_csv(src, sep="\t")
w = df.pivot(index="Mutation", columns="Replica", values=metrics)
if "WT" not in w.index:
    sys.exit("WT row missing")
mutants = [m for m in df.Mutation.unique() if m != "WT"]
print("mutants = %d | metrics = %d | expected conclusions = %d"
      % (len(mutants), len(metrics), len(mutants) * len(metrics)))

results = []
for mut in mutants:
    for met in metrics:
        wt1, wt2 = w.loc["WT", (met, "r1")], w.loc["WT", (met, "r2")]
        m1, m2 = w.loc[mut, (met, "r1")], w.loc[mut, (met, "r2")]
        if pd.isna([wt1, wt2, m1, m2]).any():
            continue
        d11, d12 = m1 - wt1, m1 - wt2
        d21, d22 = m2 - wt1, m2 - wt2
        if min(m1, m2) > max(wt1, wt2):
            cls = "ROBUST_INCREASE"
        elif max(m1, m2) < min(wt1, wt2):
            cls = "ROBUST_DECREASE"
        else:
            cls = "FRAGILE"
        results.append({
            "Mutation": mut, "Metric": met,
            "WT_r1": wt1, "WT_r2": wt2, "Mut_r1": m1, "Mut_r2": m2,
            "d_Mr1_WTr1": d11, "d_Mr1_WTr2": d12,
            "d_Mr2_WTr1": d21, "d_Mr2_WTr2": d22,
            "Classification": cls,
        })

res = pd.DataFrame(results)
res.round(4).to_csv(outdir / "replica_fragility_160_conclusions.tsv", sep="\t", index=False)

summary = (res.groupby("Metric")["Classification"].value_counts()
             .unstack(fill_value=0).reindex(metrics))
for c in ["ROBUST_INCREASE", "ROBUST_DECREASE", "FRAGILE"]:
    if c not in summary.columns:
        summary[c] = 0
summary = summary[["ROBUST_INCREASE", "ROBUST_DECREASE", "FRAGILE"]].fillna(0).astype(int)
summary["Total"] = summary.sum(axis=1)
summary["Fragility_percent"] = (100 * summary["FRAGILE"] / summary["Total"]).round(1)
summary.to_csv(outdir / "replica_fragility_summary.tsv", sep="\t")

robust = int((res.Classification != "FRAGILE").sum())
fragile = int((res.Classification == "FRAGILE").sum())
total = len(res)
print("\nSUMMARY BY METRIC")
print(summary.to_string())
print("\nOVERALL")
print("Robust conclusions  = %d / %d" % (robust, total))
print("Fragile conclusions = %d / %d" % (fragile, total))
print("Fragility rate      = %.1f%%" % (100 * fragile / total))
print("\nSaved:")
print(outdir / "replica_fragility_160_conclusions.tsv")
print(outdir / "replica_fragility_summary.tsv")
