import MDAnalysis as mda
import numpy as np, pandas as pd, os

PBC = "/home/raza/p53_thermal_xai/data/structures/c36_pbcfix"
NC  = "/home/raza/p53_thermal_xai/results/features/native_contacts.npz"
OUT = "/home/raza/p53_thermal_xai/results/features"
EQUIL_NS, FRAME_PS = 20.0, 10.0
BETA, LAMB = 5.0, 1.8

SYSTEMS = ["WT","R175H","Y220C","G245S","R248Q","R248W","R273H","R282W",
           "R175A","F134L","R249S","M237I","C242S","V143A","L145Q",
           "P151S","V157F","I195T","I232T","I255F","F270C"]

nc = np.load(NC)
resid_i, name_i = nc['resid_i'], nc['name_i']
resid_j, name_j = nc['resid_j'], nc['name_j']
r0 = nc['r0']
all_resids = nc['all_resids']
skip = int(EQUIL_NS*1000/FRAME_PS)
res_index = {r:k for k,r in enumerate(all_resids)}
nres = len(all_resids)
npair = len(r0)

rows, perres = [], {}
missing_log = []

for s in SYSTEMS:
    for rep in ['r1','r2']:
        g, x = f"{PBC}/{s}_{rep}.gro", f"{PBC}/{s}_{rep}.xtc"
        if not os.path.exists(x):
            print(f"MISSING {s}_{rep}", flush=True); continue
        u = mda.Universe(g, x)
        p = u.select_atoms("protein and not name H*")

        # build (resid,name)->index lookup for THIS system
        lut = {(rid,nm):k for k,(rid,nm) in enumerate(zip(p.resids, p.names))}

        idx_i = np.full(npair, -1, dtype=int)
        idx_j = np.full(npair, -1, dtype=int)
        for k in range(npair):
            idx_i[k] = lut.get((resid_i[k], name_i[k]), -1)
            idx_j[k] = lut.get((resid_j[k], name_j[k]), -1)

        valid = (idx_i >= 0) & (idx_j >= 0)
        n_broken_atom = (~valid).sum()
        if n_broken_atom:
            missing_log.append((s, rep, int(n_broken_atom)))

        vi, vj, vr0 = idx_i[valid], idx_j[valid], r0[valid]
        v_resid_i = np.array([res_index[r] for r in resid_i[valid]])
        v_resid_j = np.array([res_index[r] for r in resid_j[valid]])

        Qs, acc, cnt = [], np.zeros(nres), 0
        for ts in u.trajectory[skip:]:
            pos = p.positions
            d = np.linalg.norm(pos[vi]-pos[vj], axis=1)
            q = 1.0/(1.0+np.exp(BETA*(d - LAMB*vr0)))
            # pairs whose atom vanished at the mutation site count as broken (q=0)
            q_full = np.zeros(npair)
            q_full[valid] = q
            Qs.append(q_full.mean())
            acc += np.bincount(v_resid_i, weights=q, minlength=nres)
            acc += np.bincount(v_resid_j, weights=q, minlength=nres)
            cnt += 1

        Q = float(np.mean(Qs))
        rows.append(dict(system=s, replicate=rep, n_frames=cnt,
                         Q_mean=Q, Q_sd=float(np.std(Qs)),
                         n_pairs_atom_missing=int(n_broken_atom)))
        perres[f"{s}_{rep}"] = dict(zip(all_resids, acc/max(cnt,1)))
        print(f"{s}_{rep}: frames={cnt}  Q={Q:.4f}  missing_atom_pairs={n_broken_atom}", flush=True)

pd.DataFrame(rows).to_csv(f"{OUT}/global_Q.csv", index=False)
pd.DataFrame(perres).T.to_csv(f"{OUT}/per_residue_Q.csv")
print("Saved global_Q.csv and per_residue_Q.csv")
print("Systems with missing-atom contacts:", missing_log)
