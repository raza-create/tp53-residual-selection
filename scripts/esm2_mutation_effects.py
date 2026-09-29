#!/usr/bin/env python3
"""ESM-2 Mutation Effect Predictions for p53 Cancer Mutations"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
from transformers import AutoTokenizer, EsmForMaskedLM
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

output_dir = Path("/home/raza/p53_thermal_xai/results/figures/xai")
output_dir.mkdir(parents=True, exist_ok=True)

print("="*70)
print("ESM-2 MUTATION EFFECT PREDICTIONS FOR p53")
print("="*70)

sequence = "SSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNT"

CANCER_MUTATIONS = [
    ('R', 81, 'H', 4.8, 'structural'),
    ('R', 154, 'Q', 3.6, 'contact'),
    ('R', 154, 'W', 2.8, 'contact'),
    ('R', 179, 'H', 3.0, 'contact'),
    ('R', 188, 'W', 1.8, 'structural'),
    ('G', 151, 'S', 1.5, 'structural'),
    ('Y', 126, 'C', 1.2, 'structural'),
]

print("Loading ESM-2...")
tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
model = EsmForMaskedLM.from_pretrained("facebook/esm2_t6_8M_UR50D")
model.eval()

def compute_llr(seq, pos, wt, mut):
    seq_list = list(seq)
    seq_list[pos] = tokenizer.mask_token
    inputs = tokenizer("".join(seq_list), return_tensors="pt")
    mask_idx = (inputs.input_ids == tokenizer.mask_token_id).nonzero(as_tuple=True)[1].item()
    with torch.no_grad():
        probs = torch.softmax(model(**inputs).logits[0, mask_idx], dim=-1)
    wt_p = probs[tokenizer.convert_tokens_to_ids(wt)].item()
    mut_p = probs[tokenizer.convert_tokens_to_ids(mut)].item()
    return np.log(mut_p/wt_p) if wt_p > 0 else -10

print("\nAnalyzing mutations...")
results = []
for wt, pos, mut, freq, mtype in CANCER_MUTATIONS:
    if pos < len(sequence) and sequence[pos] == wt:
        llr = compute_llr(sequence, pos, wt, mut)
        name = f"{wt}{pos+94}{mut}"
        results.append({"mutation": name, "LLR": llr, "frequency": freq, "type": mtype})
        tag = "[DELETERIOUS]" if llr < -1 else ""
        print(f"  {name}: LLR={llr:.3f} {tag}")

df = pd.DataFrame(results).sort_values("LLR")
fig, ax = plt.subplots(figsize=(10,6))
colors = [("#E74C3C" if t=="structural" else "#3498DB") for t in df["type"]]
ax.barh(df["mutation"], df["LLR"], color=colors)
ax.axvline(-1, color="red", linestyle="--")
ax.set_xlabel("Log-Likelihood Ratio")
ax.set_title("ESM-2 Mutation Effect Predictions", fontweight="bold")
plt.tight_layout()
plt.savefig(output_dir/"mutation_effects.png", dpi=200)
df.to_csv(output_dir/"mutation_predictions.csv", index=False)
print(f"\nSaved: mutation_effects.png, mutation_predictions.csv")
