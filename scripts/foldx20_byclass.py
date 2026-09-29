import pandas as pd, numpy as np
from scipy.stats import spearmanr

merged = pd.read_csv("/home/raza/p53_thermal_xai/results/features/foldx20_merged.csv")
test = merged[merged['verified']=='yes'].dropna(subset=['ddG_foldx','ddG_kcal_mol'])

def perm_test(x, y, n=10000, seed=0):
    if len(x) < 4:
        return np.nan, np.nan
    rho, _ = spearmanr(x, y)
    rng = np.random.default_rng(seed)
    null = np.array([spearmanr(x, rng.permutation(y))[0] for _ in range(n)])
    p = (np.sum(np.abs(null) >= np.abs(rho)) + 1) / (n + 1)
    return rho, p

print(test[['mutant','class','ddG_kcal_mol','ddG_foldx']].to_string(index=False))

for cls in ['contact', None]:
    if cls == 'contact':
        sub = test[test['class']=='contact']
        label = "DNA-contact class only"
    else:
        sub = test[test['class']!='contact']
        label = "non-contact (structural/DNA_region/Zn/beta_sandwich)"
    rho, p = perm_test(sub['ddG_foldx'].values, sub['ddG_kcal_mol'].values)
    print(f"\n{label}: n={len(sub)}  rho={rho:+.3f}  perm_p={p}")
