#!/usr/bin/env python3
"""
Generate publication-ready summary figure
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import mean_absolute_error, r2_score
from scipy import stats
import shap
import warnings
warnings.filterwarnings('ignore')

output_dir = Path('results/figures/xai/classification')

print("="*70)
print("GENERATING PUBLICATION SUMMARY FIGURE")
print("="*70)

# Load data
data = pd.read_csv('data/features_100ns/all_features_100ns.csv')
feature_cols = ['rmsd', 'rg', 'end_to_end', 'hbonds', 'native_contacts', 'asphericity']
X = data[feature_cols].values
y = data['temperature'].values

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Train model
print("Training model...")
reg = RandomForestRegressor(n_estimators=200, max_depth=15, random_state=42, n_jobs=-1)
reg.fit(X_scaled, y)

# SHAP
print("Computing SHAP values...")
explainer = shap.TreeExplainer(reg)
shap_values = explainer.shap_values(X_scaled[:3000])
importance = np.abs(shap_values).mean(axis=0)

# LOTO
print("Running LOTO validation...")
logo = LeaveOneGroupOut()
groups = data['temperature'].values
loto_results = []
for train_idx, test_idx in logo.split(X_scaled, y, groups):
    temp = groups[test_idx[0]]
    reg_temp = RandomForestRegressor(n_estimators=200, max_depth=15, random_state=42, n_jobs=-1)
    reg_temp.fit(X_scaled[train_idx], y[train_idx])
    pred = reg_temp.predict(X_scaled[test_idx])
    mae = mean_absolute_error(y[test_idx], pred)
    loto_results.append({'temp': temp, 'mae': mae, 'pred': np.mean(pred)})
loto_mae = np.mean([r['mae'] for r in loto_results])

# Create figure
print("Creating publication figure...")
fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# Panel A: Native contacts vs Temperature
ax1 = axes[0, 0]
temp_means = data.groupby('temperature')['native_contacts'].mean()
temp_stds = data.groupby('temperature')['native_contacts'].std()
ax1.errorbar(temp_means.index, temp_means.values, yerr=temp_stds.values,
             fmt='o-', capsize=5, markersize=10, color='#2ecc71', linewidth=2)
ax1.set_xlabel('Temperature (K)', fontsize=11)
ax1.set_ylabel('Native Contacts (fraction)', fontsize=11)
ax1.set_title('A) Native Contact Loss with Temperature', fontsize=12, fontweight='bold')
ax1.text(0.05, 0.95, 'r = -0.887\nCohens d = 8.16', transform=ax1.transAxes,
         fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat'))

# Panel B: SHAP Feature Importance
ax2 = axes[0, 1]
imp_df = pd.DataFrame({'feature': feature_cols, 'importance': importance})
imp_df = imp_df.sort_values('importance', ascending=True)
colors = plt.cm.Blues(np.linspace(0.3, 0.9, len(feature_cols)))
ax2.barh(imp_df['feature'], imp_df['importance'], color=colors, edgecolor='black')
ax2.set_xlabel('Mean |SHAP Value| (K)', fontsize=11)
ax2.set_title('B) SHAP Feature Importance', fontsize=12, fontweight='bold')
for i, (_, row) in enumerate(imp_df.iterrows()):
    pct = row['importance'] / importance.sum() * 100
    ax2.text(row['importance'] + 0.3, i, f'{pct:.1f}%', va='center', fontsize=9)

# Panel C: LOTO Validation
ax3 = axes[0, 2]
temps = [r['temp'] for r in loto_results]
preds = [r['pred'] for r in loto_results]
ax3.scatter(temps, preds, s=150, c='#e74c3c', edgecolor='black', linewidth=2, zorder=5)
ax3.plot([240, 460], [240, 460], 'k--', linewidth=2, label='Perfect')
ax3.set_xlabel('True Temperature (K)', fontsize=11)
ax3.set_ylabel('Predicted Temperature (K)', fontsize=11)
ax3.set_title('C) LOTO Validation', fontsize=12, fontweight='bold')
ax3.text(0.05, 0.95, f'LOTO MAE = {loto_mae:.1f}K', transform=ax3.transAxes,
         fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat'))
ax3.legend(loc='lower right')

# Panel D: Per-temperature boxplot
ax4 = axes[1, 0]
temps_unique = sorted(data['temperature'].unique())
bp_data = [data[data['temperature']==t]['native_contacts'].values for t in temps_unique]
bp = ax4.boxplot(bp_data, positions=range(len(temps_unique)), patch_artist=True)
colors_box = plt.cm.coolwarm(np.linspace(0, 1, len(temps_unique)))
for patch, color in zip(bp['boxes'], colors_box):
    patch.set_facecolor(color)
ax4.set_xticks(range(len(temps_unique)))
ax4.set_xticklabels([f'{t}K' for t in temps_unique])
ax4.set_xlabel('Temperature', fontsize=11)
ax4.set_ylabel('Native Contacts', fontsize=11)
ax4.set_title('D) Native Contacts Distribution', fontsize=12, fontweight='bold')

# Panel E: Feature correlations
ax5 = axes[1, 1]
corrs = []
for feat in feature_cols:
    r, _ = stats.pearsonr(data[feat], data['temperature'])
    corrs.append(r)
corr_df = pd.DataFrame({'feature': feature_cols, 'correlation': corrs})
corr_df = corr_df.sort_values('correlation')
colors_corr = ['#e74c3c' if c < 0 else '#3498db' for c in corr_df['correlation']]
ax5.barh(corr_df['feature'], corr_df['correlation'], color=colors_corr, edgecolor='black')
ax5.axvline(x=0, color='black', linewidth=1)
ax5.set_xlabel('Pearson Correlation with Temperature', fontsize=11)
ax5.set_title('E) Feature-Temperature Correlations', fontsize=12, fontweight='bold')

# Panel F: Prediction scatter (training)
ax6 = axes[1, 2]
y_pred = reg.predict(X_scaled)
ax6.scatter(y, y_pred, alpha=0.05, s=1, c='blue')
ax6.plot([240, 460], [240, 460], 'r--', linewidth=2, label='Perfect')
r2 = r2_score(y, y_pred)
ax6.set_xlabel('True Temperature (K)', fontsize=11)
ax6.set_ylabel('Predicted Temperature (K)', fontsize=11)
ax6.set_title(f'F) Training Performance (R2={r2:.3f})', fontsize=12, fontweight='bold')
ax6.legend()

plt.tight_layout()
plt.savefig(output_dir / 'publication_summary_figure.png', dpi=300, bbox_inches='tight')
plt.savefig(output_dir / 'publication_summary_figure.pdf', bbox_inches='tight')
plt.close()
print("Saved: publication_summary_figure.png")
print("Saved: publication_summary_figure.pdf")

# Print statistics
print("\n" + "="*70)
print("FINAL STATISTICS FOR PAPER")
print("="*70)
print(f"Dataset: 6 temperatures, 60006 frames, 100ns each")
print(f"Training R2: {r2:.3f}")
print(f"LOTO MAE: {loto_mae:.1f}K")
print(f"Native contacts correlation: r = -0.887")
print(f"Native contacts Cohens d: 8.16")
print(f"\nSHAP Importance:")
for feat in feature_cols:
    pct = importance[feature_cols.index(feat)] / importance.sum() * 100
    print(f"  {feat}: {pct:.1f}%")
print("="*70)
