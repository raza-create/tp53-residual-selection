#!/usr/bin/env python3
"""
Analyze p53 MD simulations at different temperatures
Extract RMSD, Rg, and visualize structural changes
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from pathlib import Path

# Paths
sim_dir = Path("/home/raza/p53_thermal_xai/data/structures/long_simulations")
output_dir = Path("/home/raza/p53_thermal_xai/results")
output_dir.mkdir(parents=True, exist_ok=True)

# Temperature list
temperatures = [248, 273, 310, 350, 400, 450, 500]
temp_labels = ['-25°C', '0°C', '37°C', '77°C', '127°C', '177°C', '227°C']

# Reference structure
ref_gro = sim_dir.parent / "npt.gro"
ref_top = sim_dir.parent / "topol.top"

# Store results
results = {}

print("=" * 60)
print("Analyzing p53 MD Simulations at Different Temperatures")
print("=" * 60)

for temp in temperatures:
    print(f"\nAnalyzing {temp}K...")
    
    # Find trajectory file
    xtc_file = sim_dir / f"md_{temp}K.xtc"
    
    if not xtc_file.exists():
        print(f"  WARNING: {xtc_file} not found, skipping...")
        continue
    
    # Load trajectory
    u = mda.Universe(str(ref_gro), str(xtc_file))
    protein = u.select_atoms("protein")
    ca_atoms = u.select_atoms("protein and name CA")
    
    print(f"  Frames: {len(u.trajectory)}")
    print(f"  Protein atoms: {len(protein)}")
    print(f"  CA atoms: {len(ca_atoms)}")
    
    # Calculate RMSD (align to first frame)
    rmsd_values = []
    rg_values = []
    
    # Reference coordinates (first frame)
    u.trajectory[0]
    ref_positions = ca_atoms.positions.copy()
    
    for ts in u.trajectory:
        # RMSD of CA atoms
        current_pos = ca_atoms.positions
        rmsd = rms.rmsd(current_pos, ref_positions, superposition=True)
        rmsd_values.append(rmsd)
        
        # Radius of gyration
        rg = protein.radius_of_gyration()
        rg_values.append(rg)
    
    # Store results
    results[temp] = {
        'rmsd': np.array(rmsd_values),
        'rg': np.array(rg_values),
        'time': np.arange(len(rmsd_values)) * 5  # 5 ps per frame
    }
    
    print(f"  RMSD: {np.mean(rmsd_values):.2f} ± {np.std(rmsd_values):.2f} Å")
    print(f"  Rg: {np.mean(rg_values):.2f} ± {np.std(rg_values):.2f} Å")

print("\n" + "=" * 60)
print("Creating Visualizations...")
print("=" * 60)

# Color scheme (cold to hot)
colors = ['#0000FF', '#00BFFF', '#00FF00', '#FFFF00', '#FFA500', '#FF4500', '#FF0000']

# Figure 1: RMSD over time for all temperatures
fig, ax = plt.subplots(figsize=(12, 6))
for i, temp in enumerate(temperatures):
    if temp in results:
        time = results[temp]['time'] / 1000  # Convert to ns
        rmsd = results[temp]['rmsd']
        ax.plot(time, rmsd, color=colors[i], label=f'{temp}K ({temp_labels[i]})', alpha=0.8)

ax.set_xlabel('Time (ns)', fontsize=12)
ax.set_ylabel('RMSD (Å)', fontsize=12)
ax.set_title('p53 DNA-Binding Domain RMSD at Different Temperatures', fontsize=14)
ax.legend(loc='upper left')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(output_dir / 'figures' / 'rmsd_vs_time_all_temps.png', dpi=150)
print("  Saved: rmsd_vs_time_all_temps.png")

# Figure 2: Radius of Gyration over time
fig, ax = plt.subplots(figsize=(12, 6))
for i, temp in enumerate(temperatures):
    if temp in results:
        time = results[temp]['time'] / 1000
        rg = results[temp]['rg']
        ax.plot(time, rg, color=colors[i], label=f'{temp}K ({temp_labels[i]})', alpha=0.8)

ax.set_xlabel('Time (ns)', fontsize=12)
ax.set_ylabel('Radius of Gyration (Å)', fontsize=12)
ax.set_title('p53 DNA-Binding Domain Rg at Different Temperatures', fontsize=14)
ax.legend(loc='upper left')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(output_dir / 'figures' / 'rg_vs_time_all_temps.png', dpi=150)
print("  Saved: rg_vs_time_all_temps.png")

# Figure 3: Average RMSD vs Temperature (Bar chart)
fig, ax = plt.subplots(figsize=(10, 6))
temps_plot = []
rmsd_means = []
rmsd_stds = []
for temp in temperatures:
    if temp in results:
        temps_plot.append(temp)
        rmsd_means.append(np.mean(results[temp]['rmsd']))
        rmsd_stds.append(np.std(results[temp]['rmsd']))

bars = ax.bar(range(len(temps_plot)), rmsd_means, yerr=rmsd_stds, 
              color=colors[:len(temps_plot)], capsize=5, edgecolor='black')
ax.set_xticks(range(len(temps_plot)))
ax.set_xticklabels([f'{t}K\n({l})' for t, l in zip(temps_plot, temp_labels[:len(temps_plot)])])
ax.set_xlabel('Temperature', fontsize=12)
ax.set_ylabel('Average RMSD (Å)', fontsize=12)
ax.set_title('Average RMSD vs Temperature', fontsize=14)
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig(output_dir / 'figures' / 'rmsd_vs_temperature_bar.png', dpi=150)
print("  Saved: rmsd_vs_temperature_bar.png")

# Figure 4: Average Rg vs Temperature (Bar chart)
fig, ax = plt.subplots(figsize=(10, 6))
rg_means = []
rg_stds = []
for temp in temperatures:
    if temp in results:
        rg_means.append(np.mean(results[temp]['rg']))
        rg_stds.append(np.std(results[temp]['rg']))

bars = ax.bar(range(len(temps_plot)), rg_means, yerr=rg_stds,
              color=colors[:len(temps_plot)], capsize=5, edgecolor='black')
ax.set_xticks(range(len(temps_plot)))
ax.set_xticklabels([f'{t}K\n({l})' for t, l in zip(temps_plot, temp_labels[:len(temps_plot)])])
ax.set_xlabel('Temperature', fontsize=12)
ax.set_ylabel('Average Radius of Gyration (Å)', fontsize=12)
ax.set_title('Average Rg vs Temperature', fontsize=14)
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig(output_dir / 'figures' / 'rg_vs_temperature_bar.png', dpi=150)
print("  Saved: rg_vs_temperature_bar.png")

# Figure 5: Combined summary plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# RMSD trend
ax1 = axes[0]
ax1.errorbar(temps_plot, rmsd_means, yerr=rmsd_stds, marker='o', markersize=10,
             capsize=5, linewidth=2, color='red')
ax1.set_xlabel('Temperature (K)', fontsize=12)
ax1.set_ylabel('Average RMSD (Å)', fontsize=12)
ax1.set_title('RMSD Increases with Temperature', fontsize=14)
ax1.grid(True, alpha=0.3)

# Rg trend
ax2 = axes[1]
ax2.errorbar(temps_plot, rg_means, yerr=rg_stds, marker='s', markersize=10,
             capsize=5, linewidth=2, color='blue')
ax2.set_xlabel('Temperature (K)', fontsize=12)
ax2.set_ylabel('Average Rg (Å)', fontsize=12)
ax2.set_title('Rg Changes with Temperature', fontsize=14)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / 'figures' / 'temperature_summary.png', dpi=150)
print("  Saved: temperature_summary.png")

# Save summary table
summary_df = pd.DataFrame({
    'Temperature_K': temps_plot,
    'Temperature_C': [t - 273 for t in temps_plot],
    'RMSD_mean_A': rmsd_means,
    'RMSD_std_A': rmsd_stds,
    'Rg_mean_A': rg_means,
    'Rg_std_A': rg_stds
})
summary_df.to_csv(output_dir / 'tables' / 'temperature_analysis_summary.csv', index=False)
print("  Saved: temperature_analysis_summary.csv")

print("\n" + "=" * 60)
print("SUMMARY TABLE")
print("=" * 60)
print(summary_df.to_string(index=False))

print("\n" + "=" * 60)
print("Analysis Complete!")
print("=" * 60)
