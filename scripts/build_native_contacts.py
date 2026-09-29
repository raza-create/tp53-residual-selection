import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.distances import distance_array

REF = "/home/raza/p53_thermal_xai/data/structures/c36_pbcfix/WT_r1.gro"
OUT = "/home/raza/p53_thermal_xai/results/features/native_contacts.npz"
CUTOFF = 8.0
MIN_SEQSEP = 3

u = mda.Universe(REF)
prot = u.select_atoms("protein and not name H*")
pos = prot.positions
resid = prot.resids
name = prot.names          # atom name, e.g. CA, CB, N, O
d = distance_array(pos, pos)

n = len(pos)
i_idx, j_idx = np.triu_indices(n, k=1)
mask = (np.abs(resid[i_idx]-resid[j_idx]) > MIN_SEQSEP) & (d[i_idx,j_idx] < CUTOFF)

pair_i, pair_j = i_idx[mask], j_idx[mask]
r0 = d[pair_i, pair_j]

# store as (resid, atomname) keys instead of raw positional index
resid_i, name_i = resid[pair_i], name[pair_i]
resid_j, name_j = resid[pair_j], name[pair_j]
all_resids = np.unique(resid)

print(f"Reference: {REF}")
print(f"Heavy atoms: {prot.n_atoms}, residues: {prot.n_residues}")
print(f"Native contacts: {len(pair_i)}")

np.savez(OUT, resid_i=resid_i, name_i=name_i.astype(str),
         resid_j=resid_j, name_j=name_j.astype(str), r0=r0,
         all_resids=all_resids)
print(f"Saved: {OUT}")
