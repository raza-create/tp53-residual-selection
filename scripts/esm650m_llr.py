import torch, esm, numpy as np, pandas as pd
from scipy.stats import spearmanr

MUTATIONS = {"R175H":('R','H',175),"Y220C":('Y','C',220),"G245S":('G','S',245),
             "R248Q":('R','Q',248),"R248W":('R','W',248),"R273H":('R','H',273),
             "R282W":('R','W',282),"R175A":('R','A',175),"F134L":('F','L',134),
             "R249S":('R','S',249),"M237I":('M','I',237),"C242S":('C','S',242),
             "V143A":('V','A',143),"L145Q":('L','Q',145),"P151S":('P','S',151),
             "V157F":('V','F',157),"I195T":('I','T',195),"I232T":('I','T',232),
             "I255F":('I','F',255),"F270C":('F','C',270)}

with open("/home/raza/p53_thermal_xai/data/p53_uniprot.fasta") as f:
    lines = f.readlines()
SEQ = "".join(l.strip() for l in lines if not l.startswith(">"))
print(f"SEQ length: {len(SEQ)}")

mismatches = 0
for mut, (wt_aa, mut_aa, resid) in MUTATIONS.items():
    if SEQ[resid-1] != wt_aa:
        print(f"MISMATCH at {mut}: expected {wt_aa}, sequence has {SEQ[resid-1]} at position {resid}")
        mismatches += 1
if mismatches == 0:
    print("All 20 mutation sites verified against sequence - proceeding")
else:
    raise SystemExit(f"{mismatches} mismatches found - STOPPING, do not trust results")

model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model.eval()
batch_converter = alphabet.get_batch_converter()

_, _, tokens = batch_converter([("wt", SEQ)])
with torch.no_grad():
    logits = model(tokens)['logits'][0]
logprobs = torch.log_softmax(logits, dim=-1)

rows = []
for mut, (wt_aa, mut_aa, resid) in MUTATIONS.items():
    tok_pos = resid  # +1 for BOS, but resid is 1-indexed so this is correct
    wt_tok = alphabet.get_idx(wt_aa)
    mut_tok = alphabet.get_idx(mut_aa)
    llr = (logprobs[tok_pos, mut_tok] - logprobs[tok_pos, wt_tok]).item()
    rows.append(dict(mutant=mut, LLR_esm650M=llr))

df = pd.DataFrame(rows)
df.to_csv("/home/raza/p53_thermal_xai/results/features/esm650m_llr.csv", index=False)
print(df.to_string(index=False))

exp = pd.read_csv("/home/raza/p53_thermal_xai/experimental/exp_stability.csv")
merged = exp.merge(df, on='mutant', how='left')
test = merged[merged['verified']=='yes'].dropna(subset=['LLR_esm650M','ddG_kcal_mol'])
print(f"\nn={len(test)}")

rho,_ = spearmanr(test['LLR_esm650M'], test['ddG_kcal_mol'])
rng = np.random.default_rng(0)
null = np.array([spearmanr(test['LLR_esm650M'], rng.permutation(test['ddG_kcal_mol']))[0] for _ in range(10000)])
p = (np.sum(np.abs(null)>=np.abs(rho))+1)/10001
print(f"ESM2-650M LLR vs experimental ddG: rho={rho:+.3f}  perm_p={p:.4f}")
