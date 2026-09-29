import sys; sys.path.insert(0,'/home/raza/p53_thermal_xai/analysis')
import config as c
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
import numpy as np, pandas as pd

skip = int(c.EQUIL_DISCARD_NS*1000/c.FRAME_DT_PS)
rows = []
for s in c.SYSTEMS:
    for rep in c.REPLICATES:
        gro, xtc = c.traj(s, rep)
        u = mda.Universe(gro, xtc)
        ca = u.select_atoms('protein and name CA')
        ref = u.select_atoms('protein and name CA').positions.copy()
        aligner = align.AlignTraj(u, u, select='protein and name CA',
                                   in_memory=True).run(start=skip)
        R = rms.RMSF(ca).run(start=skip)
        for resid, val in zip(ca.resids, R.results.rmsf):
            rows.append(dict(system=s, replicate=rep, resid=int(resid), rmsf=float(val)))
        print(f"{s}_{rep} done", flush=True)

df = pd.DataFrame(rows)
df.to_csv("/home/raza/p53_thermal_xai/results/features/rmsf_per_residue.csv", index=False)
print("saved")
