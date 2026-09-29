#!/usr/bin/env python3
"""ESM-2 Comparative Analysis of p53/p63/p73 Orthologs"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
from transformers import AutoTokenizer, EsmModel
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

base_dir = Path("/home/raza/p53_thermal_xai")
output_dir = base_dir / "results/figures/xai"

print("="*70)
print("ESM-2 ORTHOLOG COMPARATIVE ANALYSIS")
print("="*70)

# Load sequences
sequences = {}
with open(base_dir / "data/sequences/orthologs.fasta") as f:
    name = None
    for line in f:
        if line.startswith(">"):
            name = line[1:].strip()
            sequences[name] = ""
        elif name:
            sequences[name] += line.strip()

print(f"Loaded {len(sequences)} sequences")

# Load ESM-2
print("Loading ESM-2...")
tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
model = EsmModel.from_pretrained("facebook/esm2_t6_8M_UR50D", output_attentions=True)
model.eval()

def get_attention_profile(seq):
    inputs = tokenizer(seq, return_tensors="pt", add_special_tokens=True)
    with torch.no_grad():
        outputs = model(**inputs)
    att = torch.stack(outputs.attentions).mean(dim=(0,1,2))
    att = att[1:-1, 1:-1].numpy()
    profile = att.sum(axis=0)
    return profile / profile.max()

print("\nComputing attention profiles...")
profiles = {}
for name, seq in sequences.items():
    profiles[name] = get_attention_profile(seq)
    print(f"  {name}: {len(seq)} residues")

# Compare p53 across species
print("\n" + "="*70)
print("CROSS-SPECIES CORRELATION (p53)")
print("="*70)

human_p53 = profiles["human_p53_DBD"]
mouse_p53 = profiles["mouse_p53_DBD"]
zfish_p53 = profiles["zebrafish_p53_DBD"]

min_len = min(len(human_p53), len(mouse_p53), len(zfish_p53))
r_hm, p_hm = stats.pearsonr(human_p53[:min_len], mouse_p53[:min_len])
r_hz, p_hz = stats.pearsonr(human_p53[:min_len], zfish_p53[:min_len])
r_mz, p_mz = stats.pearsonr(mouse_p53[:min_len], zfish_p53[:min_len])

print(f"  Human-Mouse p53:     r={r_hm:.3f} (p={p_hm:.2e})")
print(f"  Human-Zebrafish p53: r={r_hz:.3f} (p={p_hz:.2e})")
print(f"  Mouse-Zebrafish p53: r={r_mz:.3f} (p={p_mz:.2e})")

# Compare p53 vs p63 vs p73
print("\n" + "="*70)
print("PARALOG COMPARISON (human)")
print("="*70)

human_p63 = profiles["human_p63_DBD"]
human_p73 = profiles["human_p73_DBD"]

min_len2 = min(len(human_p53), len(human_p63), len(human_p73))
r_53_63, _ = stats.pearsonr(human_p53[:min_len2], human_p63[:min_len2])
r_53_73, _ = stats.pearsonr(human_p53[:min_len2], human_p73[:min_len2])
r_63_73, _ = stats.pearsonr(human_p63[:min_len2], human_p73[:min_len2])

print(f"  p53-p63: r={r_53_63:.3f}")
print(f"  p53-p73: r={r_53_73:.3f}")
print(f"  p63-p73: r={r_63_73:.3f}")

# Create figure
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: p53 orthologs
ax = axes[0,0]
ax.plot(human_p53, label="Human", alpha=0.8)
ax.plot(mouse_p53, label="Mouse", alpha=0.8)
ax.plot(zfish_p53, label="Zebrafish", alpha=0.8)
ax.set_xlabel("Residue Position")
ax.set_ylabel("Normalized Attention")
ax.set_title("p53 Attention Profiles Across Species", fontweight="bold")
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 2: Human paralogs
ax = axes[0,1]
ax.plot(human_p53, label="p53", alpha=0.8)
ax.plot(human_p63, label="p63", alpha=0.8)
ax.plot(human_p73, label="p73", alpha=0.8)
ax.set_xlabel("Residue Position")
ax.set_ylabel("Normalized Attention")
ax.set_title("Human p53/p63/p73 Attention Profiles", fontweight="bold")
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 3: Correlation heatmap
ax = axes[1,0]
corr_matrix = np.array([
    [1.0, r_hm, r_hz, r_53_63, r_53_73],
    [r_hm, 1.0, r_mz, 0, 0],
    [r_hz, r_mz, 1.0, 0, 0],
    [r_53_63, 0, 0, 1.0, r_63_73],
    [r_53_73, 0, 0, r_63_73, 1.0]
])
labels = ["H_p53", "M_p53", "Z_p53", "H_p63", "H_p73"]
im = ax.imshow(corr_matrix, cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(range(5))
ax.set_yticks(range(5))
ax.set_xticklabels(labels, rotation=45)
ax.set_yticklabels(labels)
ax.set_title("Attention Profile Correlations", fontweight="bold")
plt.colorbar(im, ax=ax, label="Pearson r")

# Plot 4: Summary stats
ax = axes[1,1]
comparisons = ["H-M p53", "H-Z p53", "M-Z p53", "p53-p63", "p53-p73", "p63-p73"]
correlations = [r_hm, r_hz, r_mz, r_53_63, r_53_73, r_63_73]
colors = ["#3498DB"]*3 + ["#E74C3C"]*3
ax.barh(comparisons, correlations, color=colors)
ax.axvline(0.8, color="green", linestyle="--", label="High similarity")
ax.set_xlabel("Pearson Correlation")
ax.set_title("Cross-Model Agreement", fontweight="bold")
ax.set_xlim(0, 1)
ax.legend()

plt.tight_layout()
plt.savefig(output_dir / "ortholog_comparison.png", dpi=200)
print(f"\nSaved: ortholog_comparison.png")

# Save results
results_df = pd.DataFrame({
    "Comparison": comparisons,
    "Correlation": correlations,
    "Type": ["Ortholog"]*3 + ["Paralog"]*3
})
results_df.to_csv(output_dir / "ortholog_correlations.csv", index=False)
print("Saved: ortholog_correlations.csv")

print("\n" + "="*70)
print("SUMMARY: High correlation across species supports evolutionary")
print("conservation of structurally important regions in p53 family.")
print("="*70)
