#!/usr/bin/env bash
set -u
GMX=${GMX:-gmx}
B=/home/raza/p53_thermal_xai/data/structures
OUT=/home/raza/p53_thermal_xai/analysis/qc_rmsd
LOG=$OUT/gmx_rms.log; : > "$LOG"
systems=$(ls "$B"/c36_pbcfix/*_r1.xtc | xargs -n1 basename | sed 's/_r1\.xtc$//' | sort)
for s in $systems; do
  for r in r1 r2; do
    tpr=""; for d in c36_prod c36_build; do [ -f "$B/$d/${s}_${r}.tpr" ] && { tpr="$B/$d/${s}_${r}.tpr"; break; }; done
    xtc="$B/c36_pbcfix/${s}_${r}.xtc"
    [ -z "$tpr" ] || [ ! -f "$xtc" ] && { echo "MISSING $s $r" | tee -a "$LOG"; continue; }
    echo "== $s $r" | tee -a "$LOG"
    printf "C-alpha\nC-alpha\n" | $GMX rms -s "$tpr" -f "$xtc" -o "$OUT/xvg/${s}_${r}_rmsd.xvg" -tu ns -xvg none >> "$LOG" 2>&1 \
      || echo "FAILED $s $r (see $LOG)"
  done
done
python3 - << 'PY'
import glob, os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT="/home/raza/p53_thermal_xai/analysis/qc_rmsd"
rows=[]; data={}
for f in sorted(glob.glob(f"{OUT}/xvg/*_rmsd.xvg")):
    s,r=os.path.basename(f).replace("_rmsd.xvg","").rsplit("_",1)
    a=np.loadtxt(f); t=a[:,0]; d=a[:,1]*10.0  # nm -> Å
    data[(s,r)]=(t,d)
    def ms(lo,hi):
        m=(t>=lo)&(t<=hi); return f"{d[m].mean():.2f} ± {d[m].std(ddof=1):.2f}" if m.sum()>1 else "NA"
    rows.append(dict(System=s,Rep=r,N_frames=len(t),End_ns=round(t[-1],2),
        **{"RMSD_0-20ns":ms(0,20),"RMSD_20-100ns":ms(20,100),"RMSD_last20ns":ms(t[-1]-20,t[-1])},
        Max_RMSD=round(d.max(),2)))
df=pd.DataFrame(rows)
order=["WT"]+sorted(x for x in df.System.unique() if x!="WT")
df["System"]=pd.Categorical(df.System,order); df=df.sort_values(["System","Rep"])
df.to_csv(f"{OUT}/rmsd_qc_table_350K.tsv",sep="\t",index=False)
print("Cα RMSD (Å) vs production start structure, 350 K, mean ± SD per interval")
print(df.to_string(index=False))
n=len(order); cols=4; rws=int(np.ceil(n/cols))
fig,axes=plt.subplots(rws,cols,figsize=(4*cols,2.6*rws),sharex=True,sharey=True)
for ax,s in zip(axes.flat,order):
    for r,c in (("r1","tab:blue"),("r2","tab:orange")):
        if (s,r) in data: t,d=data[(s,r)]; ax.plot(t,d,lw=0.5,color=c,label=r)
    ax.set_title(s,fontsize=10); ax.axvline(20,ls="--",lw=0.6,color="grey"); ax.axvline(80,ls=":",lw=0.6,color="grey")
for ax in list(axes.flat)[n:]: ax.axis("off")
axes.flat[0].legend(fontsize=8)
fig.supxlabel("Time (ns)"); fig.supylabel("Cα RMSD (Å)")
fig.suptitle("p53 DBD WT and mutants, CHARMM36/TIP3P, 350 K: Cα RMSD per replica")
fig.tight_layout(); fig.savefig(f"{OUT}/rmsd_timeseries_all_350K.png",dpi=300)
os.makedirs(f"{OUT}/per_system",exist_ok=True)
for s in order:
    plt.figure(figsize=(6,3))
    for r,c in (("r1","tab:blue"),("r2","tab:orange")):
        if (s,r) in data: t,d=data[(s,r)]; plt.plot(t,d,lw=0.6,color=c,label=r)
    plt.xlabel("Time (ns)"); plt.ylabel("Cα RMSD (Å)"); plt.title(f"{s} — 350 K"); plt.legend(); plt.tight_layout()
    plt.savefig(f"{OUT}/per_system/{s}_rmsd_350K.png",dpi=200); plt.close()
print(f"\nSaved: {OUT}/rmsd_qc_table_350K.tsv\n       {OUT}/rmsd_timeseries_all_350K.png\n       {OUT}/per_system/*.png")
PY
