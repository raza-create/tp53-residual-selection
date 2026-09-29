import sys; sys.path.insert(0,'/home/raza/p53_thermal_xai/analysis')
import config as c
import glob, re, numpy as np, pandas as pd
from scipy.stats import spearmanr

RD = "/home/raza/p53_thermal_xai/analysis/foldx/repaired_all"
EXP = "/home/raza/p53_thermal_xai/experimental/exp_stability.csv"

rows = []
for f in sorted(glob.glob(f"{RD}/*_ST.fxout")):
    m = re.search(r'/([A-Za-z0-9]+)_r([12])_snap(\d+)_Repair_0_ST\.fxout', f)
    if not m:
        print("UNPARSED:", f); continue
    system, rep, snap = m.group(1), f"r{m.group(2)}", int(m.group(3))
    with open(f) as fh:
        line = fh.readline().strip().split('\t')
    total_energy = float(line[1])
    rows.append(dict(system=system, replicate=rep, snap=snap, total_energy=total_energy))

df = pd.DataFrame(rows)
print(f"parsed {len(df)} / 840 fxout files")
df.to_csv("/home/raza/p53_thermal_xai/results/features/foldx_raw.csv", index=False)

sys_mean = df.groupby('system')['total_energy'].agg(['mean','std','count'])
print(sys_mean)

wt_mean = sys_mean.loc['WT','mean']
sys_mean['ddG_foldx'] = sys_mean['mean'] - wt_mean   # positive = less stable than WT

exp = pd.read_csv(EXP)
merged = exp.merge(sys_mean[['ddG_foldx']], left_on='mutant', right_index=True, how='left')
merged.to_csv("/home/raza/p53_thermal_xai/results/features/foldx_merged.csv", index=False)
print(merged[['mutant','ddG_kcal_mol','ddG_foldx','verified']].to_string(index=False))

test = merged[merged['verified']=='yes'].dropna(subset=['ddG_foldx','ddG_kcal_mol'])
print(f"\nn = {len(test)}")

def perm_test(x, y, n=10000, seed=0):
    rho, _ = spearmanr(x, y)
    rng = np.random.default_rng(seed)
    null = np.array([spearmanr(x, rng.permutation(y))[0] for _ in range(n)])
    p = (np.sum(np.abs(null) >= np.abs(rho)) + 1) / (n + 1)
    return rho, p

rho, p = perm_test(test['ddG_foldx'].values, test['ddG_kcal_mol'].values)
print(f"\nFoldX ddG vs experimental ddG: rho={rho:+.3f}  perm_p={p:.4f}")
