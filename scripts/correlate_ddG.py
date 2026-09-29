import numpy as np, pandas as pd
from scipy.stats import spearmanr
from MDAnalysis.analysis.distances import distance_array
import MDAnalysis as mda

REF = "/home/raza/p53_thermal_xai/data/structures/c36_pbcfix/WT_r1.gro"
NC = "/home/raza/p53_thermal_xai/results/features/native_contacts.npz"
PERRES = "/home/raza/p53_thermal_xai/results/features/per_residue_Q.csv"
EXP = "/home/raza/p53_thermal_xai/experimental/exp_stability.csv"
SHELL_A = 8.0

MUT_SITE = {"WT":None,"R175H":175,"Y220C":220,"G245S":245,"R248Q":248,
            "R248W":248,"R273H":273,"R282W":282,"R175A":175,"F134L":134,
            "R249S":249,"M237I":237,"C242S":242,"V143A":143,"L145Q":145,
            "P151S":151,"V157F":157,"I195T":195,"I232T":232,"I255F":255,"F270C":270}

nc = np.load(NC)
resid_i, resid_j = nc['resid_i'], nc['resid_j']
all_resids = nc['all_resids']
res_index = {r:k for k,r in enumerate(all_resids)}
nres = len(all_resids)

degree = np.zeros(nres)
for r in resid_i: degree[res_index[r]] += 1
for r in resid_j: degree[res_index[r]] += 1
degree_s = pd.Series(degree, index=all_resids).replace(0, np.nan)  # 0-degree -> NaN, not divide-by-zero

u = mda.Universe(REF)
ca = u.select_atoms("protein and name CA")
ca_resid, ca_pos = ca.resids, ca.positions
d_ca = distance_array(ca_pos, ca_pos)
resid2idx = {r:k for k,r in enumerate(ca_resid)}

pr = pd.read_csv(PERRES, index_col=0)
pr.columns = pr.columns.astype(int)
pr_norm = pr.div(degree_s, axis=1)          # <-- the actual fix
pr_norm['system'] = [i.rsplit('_',1)[0] for i in pr_norm.index]
pr_mean = pr_norm.groupby('system').mean(numeric_only=True)
wt_q = pr_mean.loc['WT']

print("Sanity check - normalized Q should be ~0-1 range:")
print(f"  WT per-residue Q: min={wt_q.min():.3f} max={wt_q.max():.3f} mean={wt_q.mean():.3f}")

exp = pd.read_csv(EXP)
rows = []
for _, r in exp.iterrows():
    s = r['mutant']
    if s not in pr_mean.index:
        print(f"SKIP {s}: no MD data"); continue
    dq = wt_q - pr_mean.loc[s]
    site = MUT_SITE[s]
    if site is None:
        dq_local, dq_site = 0.0, 0.0
    else:
        ci = resid2idx[site]
        shell = [x for x in ca_resid[d_ca[ci] < SHELL_A] if x in dq.dropna().index]
        dq_local = float(dq[shell].mean())
        dq_site  = float(dq[site]) if site in dq.dropna().index else np.nan
    rows.append(dict(system=s, ddG=r['ddG_kcal_mol'], verified=r['verified'],
                     dQ_local=dq_local, dQ_site=dq_site,
                     dQ_max=float(dq.max()), dQ_global=float(dq.mean())))

df = pd.DataFrame(rows)
df.to_csv("/home/raza/p53_thermal_xai/results/features/merged_features.csv", index=False)
print(df.to_string(index=False))

test = df[df['verified']=='yes'].copy()
print(f"\nn = {len(test)} verified systems")

def perm_test(x, y, n=10000, seed=0):
    rho, _ = spearmanr(x, y)
    rng = np.random.default_rng(seed)
    null = np.array([spearmanr(x, rng.permutation(y))[0] for _ in range(n)])
    p = (np.sum(np.abs(null) >= np.abs(rho)) + 1) / (n + 1)
    return rho, p

print("\nFeature       rho      perm_p")
for feat in ['dQ_local','dQ_site','dQ_max','dQ_global']:
    rho, p = perm_test(test[feat].values, test['ddG'].values)
    print(f"{feat:12s}  {rho:+.3f}   {p:.4f}")
