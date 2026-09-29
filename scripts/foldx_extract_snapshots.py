import sys; sys.path.insert(0,'/home/raza/p53_thermal_xai/analysis')
import config as c
import MDAnalysis as mda
import numpy as np, os

OUT = "/home/raza/p53_thermal_xai/analysis/foldx/snapshots_all"
os.makedirs(OUT, exist_ok=True)
N_SNAP = 20

for s in c.SYSTEMS:
    for rep in c.REPLICATES:
        gro, xtc = c.traj(s, rep)
        u = mda.Universe(gro, xtc)
        prot = u.select_atoms('protein')
        f0 = int(c.EQUIL_DISCARD_NS*1000/u.trajectory.dt)
        frames = np.linspace(f0, len(u.trajectory)-1, N_SNAP, dtype=int)
        for k, fi in enumerate(frames):
            u.trajectory[fi]
            prot.write(f'{OUT}/{s}_{rep}_snap{k:02d}.pdb')
        print(f'{s}_{rep}: wrote {N_SNAP} snapshots', flush=True)

print('total files:', len(os.listdir(OUT)))
