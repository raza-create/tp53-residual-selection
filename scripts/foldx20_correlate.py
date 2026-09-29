import glob, re, numpy as np, pandas as pd
from scipy.stats import spearmanr

RD = "/home/raza/p53_thermal_xai/analysis/foldx/buildmodel_clean"
EXP = "/home/raza/p53_thermal_xai/experimental/exp_stability.csv"

CODE2SYS = {"RX175H":"R175H","YX220C":"Y220C","GX245S":"G245S","RX248Q":"R248Q",
            "RX248W":"R248W","RX273H":"R273H","RX282W":"R282W","RX175A":"R175A",
            "FX134L":"F134L","RX249S":"R249S","MX237I":"M237I","CX242S":"C242S",
            "VX143A":"V143A","LX145Q":"L145Q","PX151S":"P151S","VX157F":"V157F",
            "IX195T":"I195T","IX232T":"I232T","IX255F":"I255F","FX270C":"F270C"}
MUT_ORDER = list(CODE2SYS.keys())

rows = []
for f in sorted(glob.glob(f"{RD}/Dif_WT_r*_snap*_Repair.fxout")):
    m = re.search(r'Dif_WT_(r[12])_snap(\d+)_Repair\.fxout', f)
    rep, snap = m.group(1), int(m.group(2))
    with open(f) as fh:
        lines = [l for l in fh if l.startswith('WT_')]
    for i, line in enumerate(lines):
        parts = line.strip().split('\t')
        total_energy = float(parts[1])
        mut_idx = i % 20
        rows.append(dict(system=CODE2SYS[MUT_ORDER[mut_idx]], replicate=rep,
                         snap=snap, ddG_foldx=total_energy))

df = pd.DataFrame(rows)
print(f"parsed rows: {len(df)} (expect 40*20=800)")
df.to_csv("/home/raza/p53_thermal_xai/results/features/foldx20_raw.csv", index=False)

sys_mean = df.groupby('system')['ddG_foldx'].agg(['mean','std','count'])
print(sys_mean)

exp = pd.read_csv(EXP)
merged = exp.merge(sys_mean[['mean']], left_on='mutant', right_index=True, how='left')
merged = merged.rename(columns={'mean':'ddG_foldx'})
merged.to_csv("/home/raza/p53_thermal_xai/results/features/foldx20_merged.csv", index=False)
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
print(f"\nFoldX BuildModel ddG vs experimental ddG: rho={rho:+.3f}  perm_p={p:.4f}")
