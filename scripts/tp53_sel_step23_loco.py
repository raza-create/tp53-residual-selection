#!/usr/bin/env python3
# Step 23: leave-one-cancer-type-out on the FROZEN donor probabilities. Nothing recomputed.
import numpy as np, pandas as pd
from pathlib import Path

B = Path("/home/raza/p53_thermal_xai/data/external/tp53_selection")
A23 = {"R158L","V157F","H179R","R273L","C176F","G266V","Y220C","I195T","R175H","C242F",
       "R249S","R248Q","D281Y","R175C","E171K","D281E","D281V","G245D","H193R","Y163C",
       "S241F","V173L","Y234C"}

d = pd.read_csv(B / "pcawg_validation_donor_probs.tsv", sep="\t")
ev = pd.read_csv(B / "pcawg_validation_events_bridged.tsv", sep="\t")
ev["spec"] = ev.icgc_specimen_id.astype(str)
ex = pd.read_csv(B / "PCAWG_sigProfiler_SBS_signatures_in_samples.csv")
ex["Sample Names"] = ex["Sample Names"].astype(str)
ct = ex.set_index("Sample Names")["Cancer Types"]
if "ct" not in d.columns:
    d["ct"] = d.spec.astype(str).map(ct)
d["carrier"] = d.carrier.astype(int)

def pb(pv):
    N = len(pv) + 1
    w = np.exp(2j*np.pi*np.arange(N)/N)
    pmf = np.real(np.fft.fft(np.prod(1-pv[None,:]+pv[None,:]*w[:,None], axis=1))/N)
    pmf = np.clip(pmf, 0, None); return pmf/pmf.sum()
def test(sub):
    obs = int(sub.carrier.sum()); pi = sub.pi.values
    pmf = pb(pi); return len(sub), obs, pi.sum(), obs/pi.sum(), pmf[obs:].sum()

print("LEAVE-ONE-CANCER-TYPE-OUT PCAWG VALIDATION\n")
print("Removed cancer         Remaining_n   Observed   Expected      O/E            P")
rows = []
for c in sorted(d.ct.dropna().unique()):
    n, o, e, oe, p = test(d[d.ct != c])
    rows.append((c, n, o, e, oe, p))
    print("%-22s %10d %10d %10.3f %8.3f %12.4g" % (c, n, o, e, oe, p))

oes = [r[4] for r in rows]; ps = [r[5] for r in rows]
worst = rows[int(np.argmax(ps))]
print("\nMinimum leave-one-cancer O/E = %.3f" % min(oes))
print("Maximum leave-one-cancer P   = %.4g" % max(ps))
print("Cancer whose removal weakens replication most = %s (O/E %.3f, P %.4g)" % (worst[0], worst[4], worst[5]))
print("All 18 remain P < 0.05:  %s" % ("YES" if max(ps) < 0.05 else "NO"))
print("All 18 remain P < 0.001: %s" % ("YES" if max(ps) < 0.001 else "NO"))
print("\nFor reference, per-cancer contribution (not a test):")
g = d.groupby("ct").agg(n=("carrier","size"), obs=("carrier","sum"), exp=("pi","sum"))
g["O_over_E"] = g.obs/g.exp
print(g.round(3).sort_values("n", ascending=False).to_string())
