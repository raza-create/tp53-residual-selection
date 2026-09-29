#!/usr/bin/env python3
"""
P53 Thermal Classification - Low vs High Temperature
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix
from sklearn.model_selection import LeaveOneGroupOut
import shap
import warnings
warnings.filterwarnings('ignore')

plt.style.use('seaborn-v0_8-whitegrid')
output_dir = Path('results/figures/xai/classification')
output_dir.mkdir(parents=True, exist_ok=True)

print("="*70)
print("P53 THERMAL CLASSIFICATION ANALYSIS")
print("="*70)

# Load data
data = pd.read_csv('data/features_100ns/all_features_100ns.csv')
print(f"Loaded: {len(data)} frames")
print(f"Temperatures: {sorted(data['temperature'].unique())}")

# Binary classification: Low (<350K) vs High (>=350K)
THRESHOLD = 350
data['temp_class'] = (data['temperature'] >= THRESHOLD).astype(int)
print(f"\nThreshold: {THRESHOLD}K")
print(f"Low (<{THRESHOLD}K): {(data['temp_class']==0).sum()} frames")
print(f"High (>={THRESHOLD}K): {(data['temp_class']==1).sum()} frames")

# Features
feature_cols = ['rmsd', 'rg', 'end_to_end', 'hbonds', 'native_contacts', 'asphericity']
X = data[feature_cols].values
y = data['temp_class'].values
groups = data['temperature'].values

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# LOTO Classification
print("\n" + "="*70)
print("LEAVE-ONE-TEMPERATURE-OUT CLASSIFICATION")
print("="*70)

logo = LeaveOneGroupOut()
results = []

for train_idx, test_idx in logo.split(X_scaled, y, groups):
    temp = groups[test_idx[0]]
    
    clf = RandomForestClassifier(n_estimators=200, max_depth=15, random_state=42, n_jobs=-1)
    clf.fit(X_scaled[train_idx], y[train_idx])
    pred = clf.predict(X_scaled[test_idx])
    
    acc = accuracy_score(y[test_idx], pred)
    f1 = f1_score(y[test_idx], pred, average='weighted')
    
    results.append({'temperature': temp, 'accuracy': acc, 'f1': f1})
    print(f"  {temp}K: Accuracy={acc:.3f}, F1={f1:.3f}")

results_df = pd.DataFrame(results)
mean_acc = results_df['accuracy'].mean()
mean_f1 = results_df['f1'].mean()

print(f"\n*** LOTO MEAN ACCURACY: {mean_acc:.3f} ({mean_acc*100:.1f}%) ***")
print(f"*** LOTO MEAN F1 SCORE: {mean_f1:.3f} ***")
results_df.to_csv(output_dir / 'loto_classification_results.csv', index=False)

# Train final model
print("\n=== FINAL MODEL ===")
clf_final = RandomForestClassifier(n_estimators=200, max_depth=15, random_state=42, n_jobs=-1)
clf_final.fit(X_scaled, y)
print(f"Training accuracy: {accuracy_score(y, clf_final.predict(X_scaled)):.3f}")

# SHAP
print("\n=== SHAP ANALYSIS ===")
explainer = shap.TreeExplainer(clf_final)
shap_values = explainer.shap_values(X_scaled[:3000])
if isinstance(shap_values, list):
    shap_vals = shap_values[1]
else:
    shap_vals = shap_values

# Save importance
importance = np.abs(shap_vals).mean(axis=0)
imp_df = pd.DataFrame({'feature': feature_cols, 'importance': importance})
imp_df = imp_df.sort_values('importance', ascending=False)
print("\nFeature Importance:")
for _, row in imp_df.iterrows():
    print(f"  {row['feature']}: {row['importance']:.4f}")
imp_df.to_csv(output_dir / 'shap_classification_importance.csv', index=False)

# Figure 1: LOTO Accuracy
fig, ax = plt.subplots(figsize=(10, 6))
colors = ['blue' if t < THRESHOLD else 'red' for t in results_df['temperature']]
bars = ax.bar(results_df['temperature'].astype(str), results_df['accuracy'], color=colors, edgecolor='black')
ax.axhline(y=0.5, color='gray', linestyle='--', label='Random (50%)')
ax.axhline(y=mean_acc, color='green', linestyle='-', linewidth=2, label=f'Mean ({mean_acc:.1%})')
ax.set_xlabel('Temperature (K)', fontsize=12)
ax.set_ylabel('Accuracy', fontsize=12)
ax.set_title('LOTO Classification Accuracy\nBlue=Low, Red=High', fontsize=14, fontweight='bold')
ax.set_ylim(0, 1.1)
ax.legend()
for bar, acc in zip(bars, results_df['accuracy']):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02, f'{acc:.0%}', ha='center')
plt.tight_layout()
plt.savefig(output_dir / 'loto_accuracy.png', dpi=200)
plt.close()
print("\nSaved: loto_accuracy.png")

# Figure 2: Confusion Matrix
cm = confusion_matrix(y, clf_final.predict(X_scaled))
fig, ax = plt.subplots(figsize=(8, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax,
            xticklabels=['Low', 'High'], yticklabels=['Low', 'High'])
ax.set_xlabel('Predicted')
ax.set_ylabel('True')
ax.set_title('Confusion Matrix', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(output_dir / 'confusion_matrix.png', dpi=200)
plt.close()
print("Saved: confusion_matrix.png")

# Figure 3: SHAP Summary
plt.figure(figsize=(10, 6))
shap.summary_plot(shap_vals, X_scaled[:3000], feature_names=feature_cols, show=False)
plt.title('SHAP Summary - Classification', fontweight='bold')
plt.tight_layout()
plt.savefig(output_dir / 'shap_summary_classification.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved: shap_summary_classification.png")

# Summary
print("\n" + "="*70)
print("SUMMARY: REGRESSION vs CLASSIFICATION")
print("="*70)
print(f"  Regression LOTO MAE:     44K (poor - doesn't generalize)")
print(f"  Classification LOTO Acc: {mean_acc:.1%} {'✅ SUCCESS!' if mean_acc > 0.8 else '⚠️ Moderate'}")
print("="*70)

pd.DataFrame([{'threshold': THRESHOLD, 'loto_accuracy': mean_acc, 'loto_f1': mean_f1}]).to_csv(output_dir / 'classification_summary.csv', index=False)
print(f"\nAll saved to: {output_dir}")
