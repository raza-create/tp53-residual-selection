import pandas as pd, numpy as np
from scipy.stats import spearmanr

exp = pd.read_csv("/home/raza/p53_thermal_xai/experimental/exp_stability.csv")
merged_feat = pd.read_csv("/home/raza/p53_thermal_xai/results/features/merged_features.csv")
foldx = pd.read_csv("/home/raza/p53_thermal_xai/results/features/foldx20_merged.csv")
physchem = pd.read_csv("/home/raza/p53_thermal_xai/results/features/physchem_features.csv")

df = exp[['mutant','pct_wt_dna_binding','verified']].copy()
df = df.merge(merged_feat[['system','dQ_local','dQ_site']], left_on='mutant', right_on='system', how='left')
df = df.merge(foldx[['mutant','ddG_foldx']], on='mutant', how='left')
df = df.merge(physchem[['mutant','RSA']], on='mutant', how='left')

test = df[df['verified']=='yes'].dropna(subset=['pct_wt_dna_binding'])
print(test[['mutant','pct_wt_dna_binding','dQ_local','dQ_site','ddG_foldx','RSA']].to_string(index=False))
print(f"\nn = {len(test)}")

def perm_test(x, y, n=10000, seed=0):
    mask = ~(np.isnan(x)|np.isnan(y))
    x,y = x[mask], y[mask]
    rho,_ = spearmanr(x,y)
    rng = np.random.default_rng(seed)
    null = np.array([spearmanr(x, rng.permutation(y))[0] for _ in range(n)])
    p = (np.sum(np.abs(null)>=np.abs(rho))+1)/(n+1)
    return rho,p,len(x)

print(f"\n{'feature':12s} {'rho':>8s} {'perm_p':>8s} {'n':>4s}")
for col in ['dQ_local','dQ_site','ddG_foldx','RSA']:
    rho,p,n = perm_test(test[col].values, test['pct_wt_dna_binding'].values, )
    print(f"{col:12s} {rho:+8.3f} {p:8.4f} {n:4d}")
