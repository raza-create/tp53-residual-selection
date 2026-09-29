#!/usr/bin/env python3
# Prespecified endpoints: PAb240 epitope (213-217) and aggregation segment S9 (251-257) exposure.
# Validation panel: WT, R175H, R248Q, V143A x r1,r2. Every 1 ns. No thresholds, no states, no stats.
# Self-contained Shrake-Rupley (numpy/scipy); mdtraj deliberately unused (ABI mismatch with numpy 1.26).
import os, sys, numpy as np, pandas as pd, MDAnalysis as mda
from scipy.spatial import cKDTree
from MDAnalysis.coordinates.XTC import XTCReader
from MDAnalysis.lib.distances import calc_bonds

B   = "/home/raza/p53_thermal_xai/data/structures"
TOP = "/home/raza/p53_thermal_xai/analysis/qc_native_contacts/top"
OUT = "/home/raza/p53_thermal_xai/analysis/epitope_exposure"
os.makedirs(OUT, exist_ok=True)
SYS, REPS   = ["WT", "R175H", "R248Q", "V143A"], ["r1", "r2"]
PAB, S9     = list(range(213, 218)), list(range(251, 258))
STRIDE      = 100
PROBE, NPTS = 1.4, 256
CUT, SEQSEP = 4.5, 3
BINS = [(0,10),(10,20),(20,40),(40,60),(60,80),(80,100)]
VDW = {"C":1.70,"N":1.55,"O":1.52,"S":1.80,"H":1.20,"ZN":1.39,"P":1.80}

def elements(ag):
    out = []
    for a in ag:
        n = a.name.upper().strip()
        if n.startswith("ZN"):
            out.append("ZN"); continue
        c = n.lstrip("0123456789")[:1]
        out.append(c if c in VDW else "C")
    return np.array(out)

_i  = np.arange(NPTS) + 0.5
_ph = np.arccos(1 - 2*_i/NPTS)
_th = np.pi*(1+5**0.5)*_i
SPHERE = np.column_stack([np.cos(_th)*np.sin(_ph), np.sin(_th)*np.sin(_ph), np.cos(_ph)])

def sasa_atoms(xyz, radii, want):
    R = radii + PROBE
    tree = cKDTree(xyz)
    rmax = R.max()
    out = np.zeros(len(want))
    for k, i in enumerate(want):
        nb = np.array([j for j in tree.query_ball_point(xyz[i], R[i] + rmax) if j != i])
        pts = xyz[i] + R[i]*SPHERE
        if len(nb):
            d = np.linalg.norm(pts[:, None, :] - xyz[nb][None, :, :], axis=2)
            acc = (d >= R[nb][None, :]).all(1)
        else:
            acc = np.ones(NPTS, bool)
        out[k] = 4*np.pi*R[i]**2 * acc.mean()
    return out

print("# Shrake-Rupley: probe %.1f A, %d sphere points" % (PROBE, NPTS))
rows, ref = [], {}
for s in SYS:
    pdb = TOP + "/" + s + "_t0.pdb"
    if not os.path.exists(pdb):
        sys.exit("missing " + pdb)
    u = mda.Universe(pdb)
    ag = u.atoms
    el = elements(ag)
    radii = np.array([VDW[e] for e in el])
    resid = np.array([a.resid for a in ag])
    i_pab = np.where(np.isin(resid, PAB))[0]
    i_s9  = np.where(np.isin(resid, S9))[0]
    i_np  = np.array([k for k in i_s9 if el[k] in ("C", "S")])
    heavy_other = np.array([k for k in range(len(ag)) if el[k] != "H" and resid[k] not in S9
                            and min(abs(resid[k]-r) for r in S9) > SEQSEP])
    heavy_s9 = np.array([k for k in i_s9 if el[k] != "H"])
    pairs = np.array([(i, j) for i in heavy_s9 for j in heavy_other])
    d0 = calc_bonds(ag.positions[pairs[:,0]], ag.positions[pairs[:,1]])
    ref[s] = int((d0 < CUT).sum())
    a0 = sasa_atoms(ag.positions, radii, np.concatenate([i_pab, i_s9]))
    if s == "WT":
        uniq, cnt = np.unique(el, return_counts=True)
        print("# element assignment (WT):", dict(zip(uniq, cnt.tolist())))
        tot = sasa_atoms(ag.positions, radii, np.arange(len(ag))).sum()
        print("# total protein SASA at t0 = %.0f A^2 (%d atoms)" % (tot, len(ag)))
    print("# %s: PAb240 atoms %d (SASA t0 %.1f) | S9 atoms %d, nonpolar %d (SASA t0 %.1f) | S9 contacts t0 = %d"
          % (s, len(i_pab), a0[:len(i_pab)].sum(), len(i_s9), len(i_np), a0[len(i_pab):].sum(), ref[s]))
    for r in REPS:
        xtc = B + "/c36_pbcfix/" + s + "_" + r + ".xtc"
        if not os.path.exists(xtc):
            print("MISSING " + xtc); continue
        rd = XTCReader(xtc)
        n = 0
        for ts in rd[::STRIDE]:
            x = ts.positions
            a = sasa_atoms(x, radii, np.concatenate([i_pab, i_s9, i_np]))
            p  = a[:len(i_pab)].sum()
            t9 = a[len(i_pab):len(i_pab)+len(i_s9)].sum()
            n9 = a[-len(i_np):].sum()
            d = calc_bonds(x[pairs[:,0]], x[pairs[:,1]], box=ts.dimensions)
            rows.append(dict(System=s, Replica=r, Time_ns=round(ts.time/1000.0, 2),
                             PAb240_SASA_213_217=p, S9_total_SASA_251_257=t9,
                             S9_nonpolar_SASA_251_257=n9,
                             S9_shielding_contact_fraction=(d < CUT).sum()/ref[s]))
            n += 1
        rd.close()
        print("  %s %s: %d sampled frames" % (s, r, n), file=sys.stderr)

df = pd.DataFrame(rows)
df.round(4).to_csv(OUT + "/epitope_s9_timeseries_1ns_350K.tsv", sep="\t", index=False)
cols = ["PAb240_SASA_213_217","S9_total_SASA_251_257","S9_nonpolar_SASA_251_257","S9_shielding_contact_fraction"]
out = []
for (s, r), g in df.groupby(["System","Replica"], sort=False):
    for lo, hi in BINS:
        w = g[(g.Time_ns >= lo) & ((g.Time_ns < hi) if hi < 100 else (g.Time_ns <= hi))]
        rec = dict(System=s, Replica=r, Window="%d-%d" % (lo, hi), n=len(w))
        for c in cols:
            rec[c] = w[c].median()
        out.append(rec)
med = pd.DataFrame(out)
med.round(3).to_csv(OUT + "/epitope_s9_binned_medians_350K.tsv", sep="\t", index=False)
print("\nMedians per window (SASA in A^2; shielding fraction relative to each system's own t0)")
print(med.round(3).to_string(index=False))
print("\nS9 t0 reference contact counts: %s" % ref)
print("\nSaved: " + OUT + "/epitope_s9_timeseries_1ns_350K.tsv")
print("       " + OUT + "/epitope_s9_binned_medians_350K.tsv")
