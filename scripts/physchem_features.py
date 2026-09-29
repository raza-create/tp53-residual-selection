import MDAnalysis as mda
import numpy as np, pandas as pd
from scipy.stats import spearmanr

VOL = {'A':88.6,'R':173.4,'N':114.1,'D':111.1,'C':108.5,'Q':143.8,'E':138.4,
       'G':60.1,'H':153.2,'I':166.7,'L':166.7,'K':168.6,'M':162.9,'F':189.9,
       'P':112.7,'S':89.0,'T':116.1,'W':227.8,'Y':193.6,'V':140.0}
KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
      'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,
      'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}

MUTATIONS = {"R175H":('R','H',175),"Y220C":('Y','C',220),"G245S":('G','S',245),
             "R248Q":('R','Q',248),"R248W":('R','W',248),"R273H":('R','H',273),
             "R282W":('R','W',282),"R175A":('R','A',175),"F134L":('F','L',134),
             "R249S":('R','S',249),"M237I":('M','I',237),"C242S":('C','S',242),
             "V143A":('V','A',143),"L145Q":('L','Q',145),"P151S":('P','S',151),
             "V157F":('V','F',157),"I195T":('I','T',195),"I232T":('I','T',232),
             "I255F":('I','F',255),"F270C":('F','C',270)}

REF = "/home/raza/p53_thermal_xai/data/structures/c36_pbcfix/WT_r1.gro"
u = mda.Universe(REF)
prot = u.select_atoms("protein")

import subprocess, os
os.chdir("/tmp")
u.select_atoms("protein").write("/tmp/wt_ref.pdb")
try:
    import freesasa
    struct = freesasa.Structure("/tmp/wt_ref.pdb")
    result = freesasa.calc(struct)
    have_sasa = True
except ImportError:
    have_sasa = False
    print("freesasa not installed - RSA feature will be skipped")

rsa_by_resid = {}
if have_sasa:
    residue_areas = result.residueAreas()
    for chain in residue_areas:
        for resnum, area in residue_areas[chain].items():
            rsa_by_resid[int(resnum)] = area.relativeTotal

rows = []
for mut, (wt_aa, mut_aa, resid) in MUTATIONS.items():
    dvol = VOL[mut_aa] - VOL[wt_aa]
    dhyd = KD[mut_aa] - KD[wt_aa]
    rsa = rsa_by_resid.get(resid, np.nan)
    rows.append(dict(mutant=mut, dVolume=dvol, dHydrophobicity=dhyd,
                     RSA=rsa, dVol_x_burial=dvol*(1-rsa) if not np.isnan(rsa) else np.nan))

feat = pd.DataFrame(rows)
feat.to_csv("/home/raza/p53_thermal_xai/results/features/physchem_features.csv", index=False)
print(feat.to_string(index=False))

exp = pd.read_csv("/home/raza/p53_thermal_xai/experimental/exp_stability.csv")
merged = exp.merge(feat, on='mutant', how='left')
test = merged[merged['verified']=='yes'].dropna(subset=['ddG_kcal_mol'])

def perm_test(x, y, n=10000, seed=0):
    mask = ~(np.isnan(x)|np.isnan(y))
    x,y = x[mask], y[mask]
    if len(x)<5: return np.nan,np.nan,len(x)
    rho,_ = spearmanr(x,y)
    rng = np.random.default_rng(seed)
    null = np.array([spearmanr(x, rng.permutation(y))[0] for _ in range(n)])
    p = (np.sum(np.abs(null)>=np.abs(rho))+1)/(n+1)
    return rho,p,len(x)

print(f"\n{'feature':16s} {'rho':>8s} {'perm_p':>8s} {'n':>4s}")
for col in ['dVolume','dHydrophobicity','RSA','dVol_x_burial']:
    rho,p,n = perm_test(test[col].values, test['ddG_kcal_mol'].values)
    print(f"{col:16s} {rho:+8.3f} {p:8.4f} {n:4d}")
