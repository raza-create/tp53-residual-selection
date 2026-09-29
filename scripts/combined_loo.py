import pandas as pd, numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

MUT_SITE = {"R175H":175,"Y220C":220,"G245S":245,"R248Q":248,"R248W":248,
            "R273H":273,"R282W":282,"R175A":175,"F134L":134,"R249S":249,
            "M237I":237,"C242S":242,"V143A":143,"L145Q":145,"P151S":151,
            "V157F":157,"I195T":195,"I232T":232,"I255F":255,"F270C":270}

exp = pd.read_csv("/home/raza/p53_thermal_xai/experimental/exp_stability.csv")
merged_feat = pd.read_csv("/home/raza/p53_thermal_xai/results/features/merged_features.csv")
foldx = pd.read_csv("/home/raza/p53_thermal_xai/results/features/foldx20_merged.csv")
physchem = pd.read_csv("/home/raza/p53_thermal_xai/results/features/physchem_features.csv")
rmsf = pd.read_csv("/home/raza/p53_thermal_xai/results/features/rmsf_per_residue.csv")

rmsf_wt = rmsf[rmsf['system']=='WT'].groupby('resid')['rmsf'].mean()
rmsf_mut = rmsf[rmsf['system']!='WT'].groupby(['system','resid'])['rmsf'].mean().reset_index()

drmsf_rows = []
for mut, site in MUT_SITE.items():
    sub = rmsf_mut[(rmsf_mut['system']==mut) & (rmsf_mut['resid']==site)]
    if len(sub) and site in rmsf_wt.index:
        drmsf_rows.append(dict(mutant=mut, dRMSF_site=float(sub['rmsf'].iloc[0]-rmsf_wt[site])))
drmsf = pd.DataFrame(drmsf_rows)

df = exp[['mutant','ddG_kcal_mol','verified']].copy()
df = df.merge(merged_feat[['system','dQ_site']], left_on='mutant', right_on='system', how='left')
df = df.merge(foldx[['mutant','ddG_foldx']], on='mutant', how='left')
df = df.merge(physchem[['mutant','RSA']], on='mutant', how='left')
df = df.merge(drmsf, on='mutant', how='left')

test = df[df['verified']=='yes'].dropna(subset=['ddG_kcal_mol','dQ_site','ddG_foldx','RSA','dRMSF_site'])
print(f"n = {len(test)} (need all 4 features non-null)")
print(test[['mutant','ddG_kcal_mol','dQ_site','ddG_foldx','RSA','dRMSF_site']].to_string(index=False))

X = test[['dQ_site','ddG_foldx','RSA','dRMSF_site']].values
y = test['ddG_kcal_mol'].values

loo = LeaveOneOut()
preds = np.zeros(len(y))
for tr, te in loo.split(X):
    sc = StandardScaler().fit(X[tr])
    Xtr, Xte = sc.transform(X[tr]), sc.transform(X[te])
    model = Ridge(alpha=1.0).fit(Xtr, y[tr])
    preds[te] = model.predict(Xte)

rho_loo, _ = spearmanr(preds, y)
print(f"\nLOO-CV predicted vs actual: rho={rho_loo:+.3f}")

def perm_test_loo(X, y, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    real_rho, _ = spearmanr(_loo_predict(X,y), y)
    null = []
    for i in range(n):
        yp = rng.permutation(y)
        pr = _loo_predict(X, yp)
        null.append(spearmanr(pr, yp)[0])
    null = np.array(null)
    p = (np.sum(np.abs(null) >= np.abs(real_rho)) + 1) / (n + 1)
    return real_rho, p

def _loo_predict(X, y):
    loo = LeaveOneOut()
    preds = np.zeros(len(y))
    for tr, te in loo.split(X):
        sc = StandardScaler().fit(X[tr])
        Xtr, Xte = sc.transform(X[tr]), sc.transform(X[te])
        m = Ridge(alpha=1.0).fit(Xtr, y[tr])
        preds[te] = m.predict(Xte)
    return preds

rho, p = perm_test_loo(X, y, n=2000)
print(f"Permutation test on LOO-CV: rho={rho:+.3f}  p={p:.4f}")
