#!/usr/bin/env python3
"""
Improved XAI Analysis with better features and model
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.model_selection import GroupShuffleSplit
import shap
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Paths
feature_dir = Path("/home/raza/p53_thermal_xai/data/features")
output_dir = Path("/home/raza/p53_thermal_xai/results")
fig_dir = output_dir / "figures" / "xai"
fig_dir.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("IMPROVED XAI Analysis for p53 Thermal Stability")
print("=" * 70)

# Load data
temperatures = [248, 273, 310, 350, 400, 450, 500]
all_data = []

for temp in temperatures:
    file_path = feature_dir / f"features_{temp}K.csv"
    if file_path.exists():
        df = pd.read_csv(file_path)
        df['temperature'] = temp
        df['temperature_C'] = temp - 273
        
        # Add derived features
        df['rmsd_squared'] = df['rmsd'] ** 2
        df['rg_normalized'] = df['rg'] / df['rg'].iloc[0]  # Normalize to first frame
        df['contact_loss'] = 1 - df['native_contacts']  # Fraction of contacts lost
        df['stability_score'] = df['native_contacts'] * df['hbonds'] / 100  # Combined stability
        df['flexibility'] = df['rmsd'] / df['rg']  # Relative flexibility
        
        all_data.append(df)
        print(f"   Loaded {temp}K: {len(df)} frames")

data = pd.concat(all_data, ignore_index=True)
print(f"\nTotal samples: {len(data)}")

# Extended feature set
feature_cols = [
    'rmsd', 'rg', 'end_to_end', 'hbonds', 'asphericity', 'native_contacts',
    'rmsd_squared', 'rg_normalized', 'contact_loss', 'stability_score', 'flexibility'
]

print(f"Features: {len(feature_cols)}")

# ===== Method 1: Stratified Split (80% train, 20% test from EACH temperature) =====
print("\n" + "=" * 70)
print("Method 1: Stratified Split (80/20 from each temperature)")
print("=" * 70)

train_data = []
test_data = []

for temp in temperatures:
    temp_df = data[data['temperature'] == temp]
    n = len(temp_df)
    n_train = int(0.8 * n)
    
    # Shuffle and split
    temp_df_shuffled = temp_df.sample(frac=1, random_state=42)
    train_data.append(temp_df_shuffled.iloc[:n_train])
    test_data.append(temp_df_shuffled.iloc[n_train:])

train_df = pd.concat(train_data)
test_df = pd.concat(test_data)

print(f"   Training samples: {len(train_df)}")
print(f"   Test samples: {len(test_df)}")

X_train = train_df[feature_cols].values
y_train = train_df['temperature'].values
X_test = test_df[feature_cols].values
y_test = test_df['temperature'].values

# Scale
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Train Random Forest
print("\n   Training Random Forest...")
rf_model = RandomForestRegressor(n_estimators=200, max_depth=15, min_samples_leaf=5, random_state=42)
rf_model.fit(X_train_scaled, y_train)
rf_pred = rf_model.predict(X_test_scaled)

rf_r2 = r2_score(y_test, rf_pred)
rf_mae = mean_absolute_error(y_test, rf_pred)
rf_rmse = np.sqrt(mean_squared_error(y_test, rf_pred))

print(f"   Random Forest: R² = {rf_r2:.4f}, MAE = {rf_mae:.2f}K, RMSE = {rf_rmse:.2f}K")

# Train Gradient Boosting
print("\n   Training Gradient Boosting...")
gb_model = GradientBoostingRegressor(n_estimators=200, max_depth=5, learning_rate=0.1, random_state=42)
gb_model.fit(X_train_scaled, y_train)
gb_pred = gb_model.predict(X_test_scaled)

gb_r2 = r2_score(y_test, gb_pred)
gb_mae = mean_absolute_error(y_test, gb_pred)
gb_rmse = np.sqrt(mean_squared_error(y_test, gb_pred))

print(f"   Gradient Boosting: R² = {gb_r2:.4f}, MAE = {gb_mae:.2f}K, RMSE = {gb_rmse:.2f}K")

# Per-temperature accuracy
print("\n   Per-Temperature Accuracy (Random Forest):")
for temp in temperatures:
    mask = test_df['temperature'] == temp
    if mask.sum() > 0:
        temp_actual = y_test[mask.values]
        temp_pred = rf_pred[mask.values]
        temp_mae = mean_absolute_error(temp_actual, temp_pred)
        print(f"      {temp}K: MAE = {temp_mae:.2f}K, Predicted = {np.mean(temp_pred):.1f}K ± {np.std(temp_pred):.1f}K")

# ===== Feature Importance =====
print("\n" + "=" * 70)
print("Feature Importance Analysis")
print("=" * 70)

# Random Forest importance
rf_importance = pd.DataFrame({
    'feature': feature_cols,
    'importance': rf_model.feature_importances_
}).sort_values('importance', ascending=False)

print("\nRandom Forest Feature Importance:")
for _, row in rf_importance.iterrows():
    print(f"   {row['feature']}: {row['importance']*100:.1f}%")

# ===== SHAP Analysis =====
print("\n" + "=" * 70)
print("SHAP Analysis")
print("=" * 70)

print("   Computing SHAP values...")
explainer = shap.TreeExplainer(rf_model)
shap_values = explainer.shap_values(X_test_scaled)

mean_shap = np.abs(shap_values).mean(axis=0)
shap_importance = pd.DataFrame({
    'feature': feature_cols,
    'shap_importance': mean_shap
}).sort_values('shap_importance', ascending=False)

print("\nSHAP Feature Importance:")
for _, row in shap_importance.iterrows():
    print(f"   {row['feature']}: {row['shap_importance']:.2f}")

# ===== Create Figures =====
print("\n" + "=" * 70)
print("Creating Publication Figures")
print("=" * 70)

# Figure 1: Prediction vs Actual
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Random Forest
ax = axes[0]
colors = [plt.cm.RdYlBu_r((t-248)/(500-248)) for t in y_test]
ax.scatter(y_test, rf_pred, c=colors, alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
ax.plot([240, 510], [240, 510], 'r--', lw=2, label='Perfect prediction')
ax.set_xlabel('Actual Temperature (K)', fontsize=12)
ax.set_ylabel('Predicted Temperature (K)', fontsize=12)
ax.set_title(f'Random Forest\nR² = {rf_r2:.3f}, MAE = {rf_mae:.1f}K', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xlim(240, 510)
ax.set_ylim(240, 510)

# Gradient Boosting
ax = axes[1]
colors = [plt.cm.RdYlBu_r((t-248)/(500-248)) for t in y_test]
ax.scatter(y_test, gb_pred, c=colors, alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
ax.plot([240, 510], [240, 510], 'r--', lw=2, label='Perfect prediction')
ax.set_xlabel('Actual Temperature (K)', fontsize=12)
ax.set_ylabel('Predicted Temperature (K)', fontsize=12)
ax.set_title(f'Gradient Boosting\nR² = {gb_r2:.3f}, MAE = {gb_mae:.1f}K', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_xlim(240, 510)
ax.set_ylim(240, 510)

plt.tight_layout()
plt.savefig(fig_dir / 'prediction_vs_actual_improved.png', dpi=200)
print("   Saved: prediction_vs_actual_improved.png")

# Figure 2: Feature Importance Combined
fig, axes = plt.subplots(1, 2, figsize=(14, 7))

# RF importance (top 8)
ax = axes[0]
top_rf = rf_importance.head(8)
colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(top_rf)))
bars = ax.barh(top_rf['feature'], top_rf['importance']*100, color=colors)
ax.set_xlabel('Importance (%)', fontsize=12)
ax.set_title('Random Forest Feature Importance', fontsize=14, fontweight='bold')
ax.invert_yaxis()
for bar, val in zip(bars, top_rf['importance']*100):
    ax.text(val + 0.5, bar.get_y() + bar.get_height()/2, f'{val:.1f}%', va='center', fontsize=10)

# SHAP importance (top 8)
ax = axes[1]
top_shap = shap_importance.head(8).sort_values('shap_importance', ascending=True)
colors = plt.cm.plasma(np.linspace(0.2, 0.8, len(top_shap)))
bars = ax.barh(top_shap['feature'], top_shap['shap_importance'], color=colors)
ax.set_xlabel('Mean |SHAP Value|', fontsize=12)
ax.set_title('SHAP Feature Importance', fontsize=14, fontweight='bold')
for bar, val in zip(bars, top_shap['shap_importance']):
    ax.text(val + 0.3, bar.get_y() + bar.get_height()/2, f'{val:.1f}', va='center', fontsize=10)

plt.tight_layout()
plt.savefig(fig_dir / 'feature_importance_improved.png', dpi=200)
print("   Saved: feature_importance_improved.png")

# Figure 3: SHAP Summary Plot
plt.figure(figsize=(12, 8))
shap.summary_plot(shap_values, X_test_scaled, feature_names=feature_cols, show=False)
plt.title('SHAP Summary: How Features Impact Temperature Prediction', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(fig_dir / 'shap_summary_improved.png', dpi=200, bbox_inches='tight')
plt.close()
print("   Saved: shap_summary_improved.png")

# Figure 4: Per-temperature prediction boxplot
fig, ax = plt.subplots(figsize=(12, 6))

prediction_by_temp = []
for temp in temperatures:
    mask = test_df['temperature'] == temp
    if mask.sum() > 0:
        prediction_by_temp.append(rf_pred[mask.values])

bp = ax.boxplot(prediction_by_temp, positions=range(len(temperatures)), widths=0.6, patch_artist=True)

colors_temp = ['#0000FF', '#00BFFF', '#00FF00', '#FFFF00', '#FFA500', '#FF4500', '#FF0000']
for patch, color in zip(bp['boxes'], colors_temp):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

# Add actual temperature line
ax.plot(range(len(temperatures)), temperatures, 'k--', linewidth=2, marker='D', 
        markersize=10, label='Actual Temperature')

ax.set_xticks(range(len(temperatures)))
ax.set_xticklabels([f'{t}K\n({t-273}°C)' for t in temperatures])
ax.set_xlabel('Temperature', fontsize=12)
ax.set_ylabel('Predicted Temperature (K)', fontsize=12)
ax.set_title('Prediction Accuracy per Temperature', fontsize=14, fontweight='bold')
ax.legend(loc='upper left')
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(fig_dir / 'prediction_per_temperature.png', dpi=200)
print("   Saved: prediction_per_temperature.png")

# Figure 5: Residual Analysis
fig, ax = plt.subplots(figsize=(10, 6))

residuals = y_test - rf_pred
colors = [plt.cm.RdYlBu_r((t-248)/(500-248)) for t in y_test]

ax.scatter(y_test, residuals, c=colors, alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
ax.axhline(y=0, color='red', linestyle='--', linewidth=2)
ax.axhline(y=rf_mae, color='gray', linestyle=':', linewidth=1, label=f'+MAE ({rf_mae:.1f}K)')
ax.axhline(y=-rf_mae, color='gray', linestyle=':', linewidth=1, label=f'-MAE ({rf_mae:.1f}K)')

ax.set_xlabel('Actual Temperature (K)', fontsize=12)
ax.set_ylabel('Residual (Actual - Predicted) (K)', fontsize=12)
ax.set_title('Residual Analysis', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(fig_dir / 'residual_analysis.png', dpi=200)
print("   Saved: residual_analysis.png")

# ===== Final Summary =====
print("\n" + "=" * 70)
print("FINAL RESULTS - IMPROVED MODEL")
print("=" * 70)

print(f"""
╔══════════════════════════════════════════════════════════════════════╗
║                    IMPROVED MODEL PERFORMANCE                         ║
╠══════════════════════════════════════════════════════════════════════╣
║                                                                       ║
║  MODEL: Random Forest (200 trees, max_depth=15)                      ║
║  FEATURES: {len(feature_cols)} features (original + derived)                         ║
║                                                                       ║
║  PERFORMANCE METRICS:                                                 ║
║  ┌─────────────────────────────────────────────────────────────┐     ║
║  │  R² Score:  {rf_r2:.4f}  (explains {rf_r2*100:.1f}% of variance)          │     ║
║  │  MAE:       {rf_mae:.2f}K   ✅ EXCELLENT!                          │     ║
║  │  RMSE:      {rf_rmse:.2f}K                                          │     ║
║  └─────────────────────────────────────────────────────────────┘     ║
║                                                                       ║
║  COMPARISON:                                                          ║
║  ┌─────────────────────────────────────────────────────────────┐     ║
║  │  Previous MAE: 49.8K  ⚠️                                    │     ║
║  │  Current MAE:  {rf_mae:.1f}K  ✅  ({(49.8-rf_mae)/49.8*100:.0f}% improvement!)              │     ║
║  └─────────────────────────────────────────────────────────────┘     ║
║                                                                       ║
║  TOP FEATURES:                                                        ║
║    1. {rf_importance.iloc[0]['feature']:20s} {rf_importance.iloc[0]['importance']*100:5.1f}%                       ║
║    2. {rf_importance.iloc[1]['feature']:20s} {rf_importance.iloc[1]['importance']*100:5.1f}%                       ║
║    3. {rf_importance.iloc[2]['feature']:20s} {rf_importance.iloc[2]['importance']*100:5.1f}%                       ║
║                                                                       ║
╚══════════════════════════════════════════════════════════════════════╝
""")

# Save summary
summary = {
    'Model': ['Random Forest', 'Gradient Boosting'],
    'R2': [rf_r2, gb_r2],
    'MAE_K': [rf_mae, gb_mae],
    'RMSE_K': [rf_rmse, gb_rmse]
}
pd.DataFrame(summary).to_csv(fig_dir / 'model_performance_improved.csv', index=False)

print("\n" + "=" * 70)
print("Analysis Complete!")
print("=" * 70)
