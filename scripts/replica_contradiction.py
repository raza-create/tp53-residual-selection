#!/usr/bin/env python3
# Split fragile conclusions into true directional contradictions vs mere overlap.
# 20 mutants x 8 metrics = 160 mutant-vs-WT conclusions.
import sys
from pathlib import Path
import pandas as pd

src = Path("/home/raza/p53_thermal_xai/analysis/replica_icc/trajectory_summaries_20_100ns.tsv")
outdir = Path("/home/raza/p53_thermal_xai/analysis/replica_icc")
metrics = ["RMSD", "Q", "beta_core", "L1", "L2_Zn", "L3", "H2", "other_DBD"]
CATS = ["ROBUST_INCREASE", "ROBUST_DECREASE", "CONTRADICTORY", "INDETERMINATE"]

df = pd.read_csv(src, sep="\t")
w = df.pivot(index="Mutation", columns="Replica", values=metrics)
if "WT" not in w.index:
    sys.exit("WT row missing")
mutants = [m for m in df.Mutation.unique() if m != "WT"]

rows = []
for mut in mutants:
    for met in metrics:
        wt_r1, wt_r2 = w.loc["WT", (met, "r1")], w.loc["WT", (met, "r2")]
        m_r1, m_r2 = w.loc[mut, (met, "r1")], w.loc[mut, (met, "r2")]
        if pd.isna([wt_r1, wt_r2, m_r1, m_r2]).any():
            continue
        wt_min, wt_max = min(wt_r1, wt_r2), max(wt_r1, wt_r2)
        mut_min, mut_max = min(m_r1, m_r2), max(m_r1, m_r2)
        if mut_min > wt_max:
            cls = "ROBUST_INCREASE"
        elif mut_max < wt_min:
            cls = "ROBUST_DECREASE"
        elif mut_min < wt_min and mut_max > wt_max:
            cls = "CONTRADICTORY"
        else:
            cls = "INDETERMINATE"
        rows.append({"Mutation": mut, "Metric": met,
                     "WT_r1": wt_r1, "WT_r2": wt_r2, "WT_min": wt_min, "WT_max": wt_max,
                     "Mut_r1": m_r1, "Mut_r2": m_r2, "Mut_min": mut_min, "Mut_max": mut_max,
                     "below_WT_by": wt_min - mut_min, "above_WT_by": mut_max - wt_max,
                     "Classification": cls})

res = pd.DataFrame(rows)
res.round(4).to_csv(outdir / "replica_contradiction_160_conclusions.tsv", sep="\t", index=False)

summary = (res.groupby("Metric")["Classification"].value_counts()
             .unstack(fill_value=0).reindex(metrics))
for c in CATS:
    if c not in summary.columns:
        summary[c] = 0
summary = summary[CATS].fillna(0).astype(int)
summary["Total"] = summary.sum(axis=1)
summary["Contradictory_percent"] = (100 * summary["CONTRADICTORY"] / summary["Total"]).round(1)
summary.to_csv(outdir / "replica_contradiction_summary.tsv", sep="\t")

n = len(res)
print("SUMMARY BY METRIC")
print(summary.to_string())
print("\nOVERALL (n = %d)" % n)
for c in CATS:
    k = int((res.Classification == c).sum())
    print("%-16s = %3d /%d  (%.1f%%)" % (c, k, n, 100 * k / n))

con = res[res.Classification == "CONTRADICTORY"].sort_values(["Metric", "Mutation"])
print("\nCONTRADICTORY cases (n = %d)" % len(con))
if len(con):
    print(con[["Mutation", "Metric", "WT_r1", "WT_r2", "Mut_r1", "Mut_r2",
               "below_WT_by", "above_WT_by"]].round(4).to_string(index=False))
    print("\nby metric: %s" % con.Metric.value_counts().to_dict())
    print("by mutant: %s" % con.Mutation.value_counts().to_dict())
else:
    print("none")

print("\nSaved:")
print(outdir / "replica_contradiction_160_conclusions.tsv")
print(outdir / "replica_contradiction_summary.tsv")
