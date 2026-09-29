#!/usr/bin/env python3
# Step 19: PCAWG WGS SBS1536 background modifier. QC only; no model, no allele test.
import io, zipfile, numpy as np, pandas as pd, urllib.request
from pathlib import Path

P = Path("/home/raza/p53_thermal_xai/data/external/pcawg_signatures")
B = Path("/home/raza/p53_thermal_xai/data/external/tp53_selection")
CHR = "17"; COMP = {"A":"T","C":"G","G":"C","T":"A"}

print("PCAWG SBS1536 BACKGROUND QC\n")
z = zipfile.ZipFile(P / "WGS_PCAWG_2018_02_09.zip")
nm = [n for n in z.namelist() if not n.startswith("__") and "1536" in n][0]
c15 = pd.read_csv(io.BytesIO(z.read(nm)))
mt  = [c for c in c15.columns if "utation" in c][0]
pn  = [c for c in c15.columns if "entanucleotide" in c or "enta" in c][0]
c15["k1536"] = c15[mt] + "|" + c15[pn]
samp = [c for c in c15.columns if c not in (mt, pn, "k1536")]
c15 = c15.set_index("k1536")[samp]
c15.columns = [c.split("::", 1)[-1] for c in c15.columns]
print("catalogue file                  =", nm)
print("Observed WGS1536 shape          =", c15.shape)
print("Observed WGS1536 samples        =", c15.shape[1])

ex = pd.read_csv(B / "PCAWG_sigProfiler_SBS_signatures_in_samples.csv")
ex["Sample Names"] = ex["Sample Names"].astype(str)
ct_of = ex.set_index("Sample Names")["Cancer Types"]
matched = [s for s in c15.columns if s in ct_of.index]
print("Samples matched to Cancer Types =", len(matched))
print("Cancer Types represented        =", ct_of.reindex(matched).nunique())

O = pd.read_csv(P / "context_counts_GRCh37_1536_genome.csv", index_col=0).sum(1)
mt_of = pd.Series(c15.index.str.split("|").str[0], index=c15.index)
pn_of = pd.Series(c15.index.str.split("|").str[1], index=c15.index)
Ok = pd.Series(O.reindex(pn_of.values).values, index=c15.index).astype(float)
par = pd.Series((pn_of.str[1] + "[" + mt_of + "]" + pn_of.str[3]).values, index=c15.index)
print("\nGenome 1536 opportunities total =", f"{int(O.sum()):,}")
print("Number of 1536 mutation channels=", len(c15))
fail = int((pn_of.str[2] != mt_of.str[0]).sum()) + int(Ok.isna().sum())
print("Centre-base/context QC failures =", fail)

uni = pd.read_csv(B / "tp53_universe_126_292.tsv", sep="\t")
uni = uni[uni.CRISPR_RFS.notna()]
lo, hi = uni.gpos.min()-3, uni.gpos.max()+3
gseq = urllib.request.urlopen("https://grch37.rest.ensembl.org/sequence/region/human/%s:%d..%d:1?content-type=text/plain"
                              % (CHR, lo, hi), timeout=60).read().decode().strip().upper()
def pent(gp, rp, ap):
    w = gseq[gp-2-lo:gp+3-lo]
    return (w, "%s>%s" % (rp, ap)) if rp in "CT" else ("".join(COMP[b] for b in reversed(w)), "%s>%s" % (COMP[rp], COMP[ap]))
uni["k1536"] = [b + "|" + a for a, b in (pent(r.gpos, r.ref_plus, r.alt_plus) for r in uni.itertuples())]
ev = pd.read_csv(B / "pcawg_validation_events_bridged.tsv", sep="\t")
ev = ev.merge(uni[["residue","wt_aa","alt_aa","k1536"]].drop_duplicates(["residue","wt_aa","alt_aa"]),
              on=["residue","wt_aa","alt_aa"], how="left")
ev["ct"] = ev.icgc_specimen_id.astype(str).map(ct_of)
print("TP53 SNPs subtracted from background =", int(ev.k1536.notna().sum()))

VAL = sorted(ev.ct.dropna().unique())
rows, F = [], {}
for t in VAL:
    cols = [s for s in matched if ct_of[s] == t]
    N = c15[cols].sum(1).astype(float)
    sub = ev[ev.ct == t].k1536.value_counts()
    N = N.subtract(sub.reindex(N.index).fillna(0), fill_value=0).clip(lower=0)
    rows.append((t, len(cols), int(N.sum())))
    d = pd.DataFrame({"N": N, "O": Ok, "par": par})
    Phat = (d.N + 0.5)/d.groupby("par").N.transform("sum").add(8.0)
    Q = d.O/d.groupby("par").O.transform("sum")
    F[t] = Phat/Q
print("\nCancer_Type           PCAWG_WGS_samples  Background_SBS_mutations")
for t, n, m in rows: print("%-22s %10d %25s" % (t, n, f"{m:,}"))

Fm = pd.DataFrame(F); v = Fm.values.ravel()
print("\nF_tk distribution across the %d validation cancer types:" % len(VAL))
print("median          = %.4f" % np.median(v))
print("5th percentile  = %.4f" % np.percentile(v, 5))
print("95th percentile = %.4f" % np.percentile(v, 95))
print("minimum         = %.4f" % v.min())
print("maximum         = %.4f" % v.max())
print("Number F_tk < 0.25 =", int((v < 0.25).sum()))
print("Number F_tk > 4    =", int((v > 4).sum()))

Qs = (Ok/pd.DataFrame({"O":Ok,"par":par}).groupby("par").O.transform("sum"))
err = 0.0
for t in VAL:
    s = (Qs*Fm[t]).groupby(par).sum()
    err = max(err, float((s-1).abs().max()))
print("\nMaximum absolute normalization error = %.3g" % err)
Fm.to_csv(B / "pcawg_F1536_by_cancer.tsv", sep="\t")
print("Saved: pcawg_F1536_by_cancer.tsv")
