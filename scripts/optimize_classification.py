#!/usr/bin/env python3
"""
Find optimal classification threshold for p53 thermal analysis
"""
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import LeaveOneGroupOut
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("OPTIMIZING CLASSIFICATION THRESHOLD")
print("="*70)

# Load data
data = pd.read_csv('data/features_100ns/all_features_100ns.csv')
feature_cols = ['rmsd', 'rg', 'end_to_end', 'hbonds', 'native_contacts', 'asphericity']
X = data[feature_cols].values
groups = data['temperature'].values

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

temps = sorted(data['temperature'].unique())
print(f"Temperatures: {temps}")

# Test different thresholds
print("\n=== TESTING DIFFERENT THRESHOLDS ===\n")
results = []

for threshold in [273, 310, 350, 400]:
    y = (data['temperature'] >= threshold).astype(int)
    
    # Check class balance
    low_temps = [t for t in temps if t < threshold]
    high_temps = [t for t in temps if t >= threshold]
    
    if len(low_temps) < 1 or len(high_temps) < 1:
        continue
    
    # LOTO
    logo = LeaveOneGroupOut()
    loto_results = []
    
    for train_idx, test_idx in logo.split(X_scaled, y, groups):
        temp = groups[test_idx[0]]
        clf = RandomForestClassifier(n_estimators=200, max_depth=15, random_state=42, n_jobs=-1)
        clf.fit(X_scaled[train_idx], y[train_idx])
        pred = clf.predict(X_scaled[test_idx])
        acc = accuracy_score(y[test_idx], pred)
        loto_results.append({'temp': temp, 'acc': acc})
    
    mean_acc = np.mean([r['acc'] for r in loto_results])
    min_acc = np.min([r['acc'] for r in loto_results])
    
    results.append({
        'threshold': threshold,
        'low_temps': low_temps,
        'high_temps': high_temps,
        'mean_acc': mean_acc,
        'min_acc': min_acc,
        'per_temp': {r['temp']: r['acc'] for r in loto_results}
    })
    
    print(f"Threshold {threshold}K: Low={low_temps}, High={high_temps}")
    print(f"  Mean Accuracy: {mean_acc:.1%}")
    print(f"  Min Accuracy:  {min_acc:.1%}")
    print(f"  Per temperature: ", end="")
    for r in loto_results:
        print(f"{r['temp']}K={r['acc']:.0%} ", end="")
    print("\n")

# Find best
best = max(results, key=lambda x: x['min_acc'])  # Maximize worst-case
print("="*70)
print(f"BEST THRESHOLD: {best['threshold']}K")
print(f"  Mean Accuracy: {best['mean_acc']:.1%}")
print(f"  Min Accuracy:  {best['min_acc']:.1%}")
print("="*70)

# Now try 3-CLASS approach
print("\n" + "="*70)
print("ALTERNATIVE: 3-CLASS CLASSIFICATION")
print("="*70)
print("Low (248, 273) vs Medium (310, 350) vs High (400, 450)")

def assign_class(t):
    if t <= 273:
        return 0  # Low
    elif t <= 350:
        return 1  # Medium
    else:
        return 2  # High

data['temp_class_3'] = data['temperature'].apply(assign_class)
y3 = data['temp_class_3'].values

logo = LeaveOneGroupOut()
loto_3class = []

for train_idx, test_idx in logo.split(X_scaled, y3, groups):
    temp = groups[test_idx[0]]
    clf = RandomForestClassifier(n_estimators=200, max_depth=15, random_state=42, n_jobs=-1)
    clf.fit(X_scaled[train_idx], y3[train_idx])
    pred = clf.predict(X_scaled[test_idx])
    acc = accuracy_score(y3[test_idx], pred)
    loto_3class.append({'temp': temp, 'acc': acc})
    print(f"  {temp}K: Accuracy={acc:.1%}")

mean_3class = np.mean([r['acc'] for r in loto_3class])
min_3class = np.min([r['acc'] for r in loto_3class])
print(f"\n3-Class Mean Accuracy: {mean_3class:.1%}")
print(f"3-Class Min Accuracy:  {min_3class:.1%}")

# FINAL RECOMMENDATION
print("\n" + "="*70)
print("RECOMMENDATION")
print("="*70)
if best['min_acc'] > min_3class and best['min_acc'] > 0.7:
    print(f"USE BINARY: Threshold={best['threshold']}K (Min acc={best['min_acc']:.1%})")
elif min_3class > 0.7:
    print(f"USE 3-CLASS: Low/Medium/High (Min acc={min_3class:.1%})")
else:
    print("CONSIDER: Focus on EXTREME temperatures only (248K vs 450K)")
    # Test extreme only
    extreme_mask = data['temperature'].isin([248, 450])
    X_ext = X_scaled[extreme_mask]
    y_ext = (data.loc[extreme_mask, 'temperature'] >= 400).astype(int).values
    g_ext = data.loc[extreme_mask, 'temperature'].values
    
    logo = LeaveOneGroupOut()
    for train_idx, test_idx in logo.split(X_ext, y_ext, g_ext):
        temp = g_ext[test_idx[0]]
        clf = RandomForestClassifier(n_estimators=200, max_depth=15, random_state=42, n_jobs=-1)
        clf.fit(X_ext[train_idx], y_ext[train_idx])
        pred = clf.predict(X_ext[test_idx])
        acc = accuracy_score(y_ext[test_idx], pred)
        print(f"  Extreme {temp}K: {acc:.1%}")
