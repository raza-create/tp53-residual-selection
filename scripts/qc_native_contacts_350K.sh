#!/usr/bin/env bash
set -u
GMX=${GMX:-gmx}
B=/home/raza/p53_thermal_xai/data/structures
OUT=/home/raza/p53_thermal_xai/analysis/qc_native_contacts
LOG=$OUT/gmx_trjconv.log; : > "$LOG"
for s in $(ls "$B"/c36_pbcfix/*_r1.xtc | xargs -n1 basename | sed 's/_r1\.xtc$//' | sort); do
  tpr=""; for d in c36_prod c36_build; do [ -f "$B/$d/${s}_r1.tpr" ] && { tpr="$B/$d/${s}_r1.tpr"; break; }; done
  [ -z "$tpr" ] && { echo "MISSING TPR $s"; continue; }
  echo Protein | $GMX trjconv -s "$tpr" -f "$B/c36_pbcfix/${s}_r1.xtc" -dump 0 -o "$OUT/top/${s}_t0.pdb" >> "$LOG" 2>&1 \
    || echo "FAILED trjconv $s (see $LOG)"
done
python3 - << 'PY'
import glob, os, sys, numpy as np, pandas as pd, MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCReader
from MDAnalysis.lib.distances import distance_array, calc_bonds
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
B="/home/raza/p53_thermal_xai/data/structures"; OUT="/home/raza/p53_thermal_xai/analysis/qc_native_contacts"
CUT, SEQSEP, BETA, LAM = 4.5, 3, 5.0, 1.8   # Best-Hummer-Eaton heavy-atom Q; |i-j|>3
BB = {"N","CA","C","O"}
systems = sorted(os.path.basename(f)[:-7] for f in glob.glob(f"{OUT}/top/*_t0.pdb"))
if "WT" not in systems: sys.exit("WT reference PDB missing")
U = {s: mda.Universe(f"{OUT}/top/{s}_t0.pdb") for s in systems}
keys = {}
for s,u in U.items():
    ag = u.select_atoms("protein and not name H*")
    keys[s] = {(int(a.resid), a.name): (a.resname, int(a.index)) for a in ag}
    rr = sorted({k[0] for k in keys[s]})
    print(f"# {s}: heavy atoms={len(ag)} residues={len(rr)} resid {rr[0]}-{rr[-1]} total_atoms={len(u.atoms)}")
wt = keys["WT"]
mutpos = sorted({k[0] for s in systems for k,(rn,_) in keys[s].items() if k in wt and wt[k][0] != rn})
common = []
for k,(rn,_) in wt.items():
    ok = True
    for s in systems:
        v = keys[s].get(k)
        if v is None or (v[0] != rn and k[1] not in BB): ok = False; break
        if k[0] in mutpos and k[1] not in BB: ok = False; break
    if ok: common.append(k)
common.sort()
ref = U["WT"].atoms[[wt[k][1] for k in common]].positions
D = distance_array(ref, ref)
res = np.array([k[0] for k in common])
iu, ju = np.where(np.triu((D < CUT) & (np.abs(res[:,None]-res[None,:]) > SEQSEP), 1))
r0 = D[iu, ju]
pairs = pd.DataFrame({"res_i":res[iu], "atom_i":[common[i][1] for i in iu],
                      "res_j":res[ju], "atom_j":[common[j][1] for j in ju], "r0_A":np.round(r0,3)})
pairs.to_csv(f"{OUT}/WT_native_contact_set.tsv", sep="\t", index=False)
print(f"# Mutated positions (side chains excluded, backbone N/CA/C/O kept): {mutpos}")
print(f"# Common heavy atoms: {len(common)} | WT native heavy-atom contacts: {len(r0)} "
      f"(cutoff {CUT} Å, |i-j|>{SEQSEP}, beta={BETA}/Å, lambda={LAM}) | residue pairs: {len(set(zip(pairs.res_i,pairs.res_j)))}")
rows, data = [], {}
for s in systems:
    ii = np.array([keys[s][common[i]][1] for i in iu]); jj = np.array([keys[s][common[j]][1] for j in ju])
    for r in ("r1","r2"):
        xtc = f"{B}/c36_pbcfix/{s}_{r}.xtc"
        if not os.path.exists(xtc): print(f"MISSING {xtc}"); continue
        rd = XTCReader(xtc)
        if rd.n_atoms < len(U[s].atoms): print(f"ATOM MISMATCH {s} {r}: xtc {rd.n_atoms} < pdb {len(U[s].atoms)}"); continue
        t = np.empty(rd.n_frames); Q = np.empty(rd.n_frames)
        for n,ts in enumerate(rd):
            d = calc_bonds(ts.positions[ii], ts.positions[jj], box=ts.dimensions)
            Q[n] = np.mean(1.0/(1.0+np.exp(BETA*(d-LAM*r0)))); t[n] = ts.time/1000.0
        rd.close(); data[(s,r)] = (t,Q)
        def ms(lo,hi):
            m=(t>=lo)&(t<=hi); return f"{Q[m].mean():.3f} ± {Q[m].std(ddof=1):.3f}"
        rows.append({"System":s,"Rep":r,"Q_t0":round(Q[0],3),"Q_0-20ns":ms(0,20),"Q_20-100ns":ms(20,100),
                     "Q_last20ns":ms(t[-1]-20,t[-1]),"Min_Q":round(Q.min(),3),"End_ns":round(t[-1],2)})
        print(f"done {s} {r}", file=sys.stderr)
np.savez_compressed(f"{OUT}/Q_timeseries_350K.npz", **{f"{s}_{r}_t":v[0] for (s,r),v in data.items()},
                    **{f"{s}_{r}_Q":v[1] for (s,r),v in data.items()})
df = pd.DataFrame(rows); order = ["WT"]+[x for x in systems if x!="WT"]
df["System"] = pd.Categorical(df.System, order); df = df.sort_values(["System","Rep"])
df.to_csv(f"{OUT}/Q_qc_table_350K.tsv", sep="\t", index=False)
print("\nNative-contact fraction Q (common WT heavy-atom contact set), 350 K, mean ± SD per interval")
print(df.to_string(index=False))
n=len(order); cols=4; rws=int(np.ceil(n/cols))
fig,axes=plt.subplots(rws,cols,figsize=(4*cols,2.6*rws),sharex=True,sharey=True)
for ax,s in zip(axes.flat,order):
    for r,c in (("r1","tab:blue"),("r2","tab:orange")):
        if (s,r) in data: ax.plot(*data[(s,r)],lw=0.5,color=c,label=r)
    ax.set_title(s,fontsize=10); ax.axvline(20,ls="--",lw=0.6,color="grey"); ax.axvline(80,ls=":",lw=0.6,color="grey")
for ax in list(axes.flat)[n:]: ax.axis("off")
axes.flat[0].legend(fontsize=8)
fig.supxlabel("Time (ns)"); fig.supylabel("Native-contact fraction Q")
fig.suptitle("p53 DBD WT and mutants, 350 K: Q relative to common WT contact set")
fig.tight_layout(); fig.savefig(f"{OUT}/Q_timeseries_all_350K.png",dpi=300)
os.makedirs(f"{OUT}/per_system",exist_ok=True)
for s in order:
    plt.figure(figsize=(6,3))
    for r,c in (("r1","tab:blue"),("r2","tab:orange")):
        if (s,r) in data: plt.plot(*data[(s,r)],lw=0.6,color=c,label=r)
    plt.xlabel("Time (ns)"); plt.ylabel("Q"); plt.title(f"{s} — 350 K"); plt.legend(); plt.tight_layout()
    plt.savefig(f"{OUT}/per_system/{s}_Q_350K.png",dpi=200); plt.close()
print(f"\nSaved: {OUT}/Q_qc_table_350K.tsv | {OUT}/WT_native_contact_set.tsv | {OUT}/Q_timeseries_350K.npz\n       {OUT}/Q_timeseries_all_350K.png | {OUT}/per_system/*.png")
PY
