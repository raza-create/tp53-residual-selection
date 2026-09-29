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

f0 = sorted(glob.glob(f"{RD}/Dif_WT_r*_snap*_Repair.fxout"))[0]
with open(f0) as fh:
    for line in fh:
        if line.startswith('Pdb\t'):
            cols = line.strip().split('\t')[1:]
            break

rows = []
for f in sorted(glob.glob(f"{RD}/Dif_WT_r*_snap*_Repair.fxout")):
    with open(f) as fh:
        lines = [l for l in fh if l.startswith('WT_')]
    for i, line in enumerate(lines):
        parts = line.strip().split('\t')
        vals = [float(x) for x in parts[1:]]
        row = dict(zip(cols, vals))
        row['system'] = CODE2SYS[MUT_ORDER[i % 20]]
        rows.append(row)

df = pd.DataFrame(rows)
sys_mean = df.groupby('system')[cols].mean()

exp = pd.read_csv(EXP)
merged = exp.merge(sys_mean, left_on='mutant', right_index=True, how='left')
base = merged[(merged['verified']=='yes') & merged['mutant'].isin(sys_mean.index)]
print(f"n={len(base)} mutants with both FoldX and experimental ddG (WT excluded)")

def perm_test(x, y, n=10000, seed=0):
    rho, _ = spearmanr(x, y)
    if np.isnan(rho): return rho, np.nan
    rng = np.random.default_rng(seed)
    null = np.array([spearmanr(x, rng.permutation(y))[0] for _ in range(n)])
    p = (np.sum(np.abs(null) >= np.abs(rho)) + 1) / (n + 1)
    return rho, p

results = []
for col in cols:
    sub = base.dropna(subset=[col, 'ddG_kcal_mol'])
    x = sub[col].values
    if len(x) < 5 or np.std(x) < 1e-10: continue
    rho, p = perm_test(x, sub['ddG_kcal_mol'].values)
    results.append((col, rho, p, len(sub)))

results.sort(key=lambda r: r[2])
m = len(results)
print(f"\n{'term':28s} {'rho':>8s} {'perm_p':>8s} {'p_holm':>8s} {'n':>4s}")
for k,(col,rho,p,n) in enumerate(results):
    p_holm = min(1.0, p*(m-k))
    print(f"{col:28s} {rho:+8.3f} {p:8.4f} {p_holm:8.4f} {n:4d}")
