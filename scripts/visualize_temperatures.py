#!/usr/bin/env python3
"""
Visualize protein structures at different temperatures
Including 310K (37°C) - Normal Body Temperature
"""

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

# Create output directory
os.makedirs('../results/figures', exist_ok=True)
os.makedirs('../results/tables', exist_ok=True)

# Temperature configurations - INCLUDING 310K (37°C)
temperatures = {
    '298K': {'path': '../data/structures/temp_298K', 'xtc': 'md.xtc', 'color': 'blue', 'label': '298K (25°C) - Room Temp'},
    '310K': {'path': '../data/structures', 'xtc': 'md_310K_short.xtc', 'color': 'green', 'label': '310K (37°C) - Body Temp'},
    '323K': {'path': '../data/structures/temp_323K', 'xtc': 'md.xtc', 'color': 'gold', 'label': '323K (50°C) - Fever'},
    '350K': {'path': '../data/structures/temp_350K', 'xtc': 'md.xtc', 'color': 'orange', 'label': '350K (77°C) - Hot'},
    '373K': {'path': '../data/structures/temp_373K', 'xtc': 'md.xtc', 'color': 'red', 'label': '373K (100°C) - Boiling'}
}

# Reference structure
ref_gro = '../data/structures/npt.gro'

# Store data for each temperature
all_data = {}

print("="*60)
print("LOADING TRAJECTORIES AT DIFFERENT TEMPERATURES")
print("="*60)

for temp_name, config in temperatures.items():
    print(f"\nLoading {temp_name} ({config['label']})...")
    xtc_file = f"{config['path']}/{config['xtc']}"
    
    u = mda.Universe(ref_gro, xtc_file)
    protein = u.select_atoms('protein')
    ca_atoms = u.select_atoms('protein and name CA')
    
    # Get last frame positions
    u.trajectory[-1]
    positions = ca_atoms.positions.copy()
    
    # Calculate radius of gyration
    rg = protein.radius_of_gyration()
    
    all_data[temp_name] = {
        'positions': positions,
        'rg': rg,
        'universe': u,
        'ca_atoms': ca_atoms,
        'color': config['color'],
        'label': config['label']
    }
    
    print(f"  Rg = {rg:.2f} Å")

#############################################
# FIGURE 1: 3D Structures at All Temperatures (2x3 Grid)
#############################################
print("\n" + "="*60)
print("CREATING 3D STRUCTURE COMPARISON")
print("="*60)

fig = plt.figure(figsize=(18, 12))

# Hotspot positions (residue number - 91 for index)
hotspots = {
    'R175': 175-91,
    'G245': 245-91,
    'R248': 248-91,
    'R273': 273-91,
    'R282': 282-91
}

temp_order = ['298K', '310K', '323K', '350K', '373K']

for idx, temp_name in enumerate(temp_order):
    data = all_data[temp_name]
    ax = fig.add_subplot(2, 3, idx+1, projection='3d')
    
    positions = data['positions']
    color = data['color']
    label = data['label']
    rg = data['rg']
    
    # Plot backbone
    ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], 
            color=color, alpha=0.4, linewidth=1)
    
    # Plot CA atoms
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
               c=color, s=20, alpha=0.7)
    
    # Mark hotspots with offset labels to avoid overlap
    offsets = [(8, 8, 5), (-10, 8, 3), (8, -10, 3), (-10, -8, 5), (8, 8, -5)]
    
    for i, (name, hs_idx) in enumerate(hotspots.items()):
        if 0 <= hs_idx < len(positions):
            ax.scatter(positions[hs_idx, 0], positions[hs_idx, 1], positions[hs_idx, 2],
                       c='black', s=150, marker='*', edgecolors='yellow', linewidths=1)
            
            # Offset labels to avoid overlap
            offset = offsets[i % len(offsets)]
            ax.text(positions[hs_idx, 0] + offset[0], 
                    positions[hs_idx, 1] + offset[1], 
                    positions[hs_idx, 2] + offset[2], 
                    name, fontsize=8, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
    
    ax.set_xlabel('X (Å)', fontsize=9)
    ax.set_ylabel('Y (Å)', fontsize=9)
    ax.set_zlabel('Z (Å)', fontsize=9)
    ax.set_title(f'{label}\nRg = {rg:.1f} Å', fontsize=11, fontweight='bold')

# Add explanation in 6th subplot
ax = fig.add_subplot(2, 3, 6)
ax.axis('off')
explanation = """
TEMPERATURE EFFECTS ON p53 STRUCTURE

- 298K (25°C): Room temperature - Stable structure
- 310K (37°C): Body temperature - Normal function
- 323K (50°C): Fever-like - Beginning instability
- 350K (77°C): High heat - Partial unfolding
- 373K (100°C): Boiling point - Denaturation

KEY OBSERVATIONS:
- Radius of Gyration (Rg) increases with temperature
- Higher Rg = more expanded/unfolded structure
- Hotspot mutations (★) affect stability

HOTSPOT POSITIONS:
★ R175 - Zinc binding region
★ G245 - Loop L3
★ R248 - DNA contact
★ R273 - DNA contact  
★ R282 - H2 helix
"""
ax.text(0.1, 0.9, explanation, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))

plt.suptitle('p53 DNA-Binding Domain Structure at Different Temperatures', 
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig('../results/figures/structures_all_temperatures.png', dpi=200, bbox_inches='tight')
print("Saved: structures_all_temperatures.png")
plt.close()

#############################################
# FIGURE 2: Overlay All Structures
#############################################
print("\nCreating structure overlay...")

fig = plt.figure(figsize=(14, 12))
ax = fig.add_subplot(111, projection='3d')

for temp_name in temp_order:
    data = all_data[temp_name]
    positions = data['positions']
    color = data['color']
    label = data['label']
    
    ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], 
            color=color, alpha=0.7, linewidth=2, label=label)

ax.set_xlabel('X (Å)', fontsize=12)
ax.set_ylabel('Y (Å)', fontsize=12)
ax.set_zlabel('Z (Å)', fontsize=12)
ax.set_title('Structure Overlay: All Temperatures\n(Blue=Cold → Green=Body Temp → Red=Hot)', 
             fontsize=14, fontweight='bold')
ax.legend(loc='upper left', fontsize=10)

plt.tight_layout()
plt.savefig('../results/figures/structure_overlay.png', dpi=200, bbox_inches='tight')
print("Saved: structure_overlay.png")
plt.close()

#############################################
# FIGURE 3: Radius of Gyration vs Temperature
#############################################
print("\nCreating Rg vs Temperature plot...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Bar chart
ax = axes[0]
temps = [298, 310, 323, 350, 373]
temp_labels = ['298K\n(25°C)\nRoom', '310K\n(37°C)\nBody', '323K\n(50°C)\nFever', 
               '350K\n(77°C)\nHot', '373K\n(100°C)\nBoiling']
rg_values = [all_data[f'{t}K']['rg'] for t in temps]
colors = [all_data[f'{t}K']['color'] for t in temps]

bars = ax.bar(temps, rg_values, color=colors, edgecolor='black', linewidth=2, width=15)
ax.set_xlabel('Temperature (K)', fontsize=12)
ax.set_ylabel('Radius of Gyration (Å)', fontsize=12)
ax.set_title('Protein Size vs Temperature', fontsize=14, fontweight='bold')
ax.set_xticks(temps)
ax.set_xticklabels(temp_labels, fontsize=9)

# Add value labels on bars
for bar, rg in zip(bars, rg_values):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, 
            f'{rg:.1f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

# Highlight body temperature
ax.axvline(x=310, color='green', linestyle='--', alpha=0.5, linewidth=2)
ax.text(310, ax.get_ylim()[1]*0.95, 'Body\nTemp', ha='center', fontsize=9, color='green')

ax.grid(True, alpha=0.3, axis='y')

# Line plot
ax = axes[1]
ax.plot(temps, rg_values, 'ko-', linewidth=2, markersize=10)
for i, (t, rg, c) in enumerate(zip(temps, rg_values, colors)):
    ax.scatter(t, rg, c=c, s=200, zorder=5, edgecolors='black', linewidths=2)
    
ax.set_xlabel('Temperature (K)', fontsize=12)
ax.set_ylabel('Radius of Gyration (Å)', fontsize=12)
ax.set_title('Rg Trend with Temperature', fontsize=14, fontweight='bold')
ax.axvline(x=310, color='green', linestyle='--', alpha=0.5, linewidth=2, label='Body Temp (37°C)')
ax.legend()
ax.grid(True, alpha=0.3)

plt.suptitle('p53 Structural Expansion with Increasing Temperature', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('../results/figures/rg_vs_temperature.png', dpi=200, bbox_inches='tight')
print("Saved: rg_vs_temperature.png")
plt.close()

#############################################
# FIGURE 4: 2D Projections
#############################################
print("\nCreating 2D projection plots...")

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# XY projection
ax = axes[0]
for temp_name in temp_order:
    data = all_data[temp_name]
    positions = data['positions']
    ax.plot(positions[:, 0], positions[:, 1], 
            color=data['color'], alpha=0.7, linewidth=1.5, label=data['label'])
ax.set_xlabel('X (Å)', fontsize=11)
ax.set_ylabel('Y (Å)', fontsize=11)
ax.set_title('XY Projection (Top View)', fontsize=12, fontweight='bold')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# XZ projection
ax = axes[1]
for temp_name in temp_order:
    data = all_data[temp_name]
    positions = data['positions']
    ax.plot(positions[:, 0], positions[:, 2], 
            color=data['color'], alpha=0.7, linewidth=1.5, label=data['label'])
ax.set_xlabel('X (Å)', fontsize=11)
ax.set_ylabel('Z (Å)', fontsize=11)
ax.set_title('XZ Projection (Front View)', fontsize=12, fontweight='bold')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# YZ projection
ax = axes[2]
for temp_name in temp_order:
    data = all_data[temp_name]
    positions = data['positions']
    ax.plot(positions[:, 1], positions[:, 2], 
            color=data['color'], alpha=0.7, linewidth=1.5, label=data['label'])
ax.set_xlabel('Y (Å)', fontsize=11)
ax.set_ylabel('Z (Å)', fontsize=11)
ax.set_title('YZ Projection (Side View)', fontsize=12, fontweight='bold')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.suptitle('p53 Structure Projections at Different Temperatures', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('../results/figures/structure_2D_projections.png', dpi=200, bbox_inches='tight')
print("Saved: structure_2D_projections.png")
plt.close()

#############################################
# FIGURE 5: Per-Residue Displacement
#############################################
print("\nCalculating per-residue displacement...")

ref_positions = all_data['310K']['positions']  # Use body temp as reference

fig, ax = plt.subplots(figsize=(14, 6))

residue_numbers = np.arange(91, 91 + len(ref_positions))

for temp_name in ['298K', '323K', '350K', '373K']:
    data = all_data[temp_name]
    positions = data['positions']
    
    # Calculate displacement from 310K (body temp)
    displacement = np.sqrt(np.sum((positions - ref_positions)**2, axis=1))
    
    ax.plot(residue_numbers, displacement, 
            color=data['color'], linewidth=1.5, label=data['label'], alpha=0.8)

# Mark hotspots
hotspot_residues = [175, 245, 248, 273, 282]
for hr in hotspot_residues:
    ax.axvline(x=hr, color='gray', linestyle='--', alpha=0.5)
    ax.text(hr, ax.get_ylim()[1]*0.95, f'{hr}', ha='center', fontsize=8, 
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))

ax.set_xlabel('Residue Number', fontsize=12)
ax.set_ylabel('Displacement from 310K Body Temp (Å)', fontsize=12)
ax.set_title('Per-Residue Structural Change Compared to Body Temperature (310K)', 
             fontsize=14, fontweight='bold')
ax.legend(loc='upper right', fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('../results/figures/residue_displacement.png', dpi=200, bbox_inches='tight')
print("Saved: residue_displacement.png")
plt.close()

#############################################
# Save PDB files for each temperature
#############################################
print("\nSaving PDB files for 3D visualization software...")

for temp_name, data in all_data.items():
    u = data['universe']
    protein = u.select_atoms('protein')
    u.trajectory[-1]
    
    pdb_file = f'../results/figures/p53_{temp_name}.pdb'
    protein.write(pdb_file)
    print(f"Saved: p53_{temp_name}.pdb")

#############################################
# Create Summary Table
#############################################
print("\n" + "="*60)
print("STRUCTURAL METRICS SUMMARY")
print("="*60)

print("\n{:<12} {:<25} {:<12}".format('Temperature', 'Condition', 'Rg (Å)'))
print("-"*50)
for temp_name in temp_order:
    data = all_data[temp_name]
    label = data['label'].split(' - ')[1] if ' - ' in data['label'] else ''
    print("{:<12} {:<25} {:<12.2f}".format(temp_name, label, data['rg']))

# Save summary CSV
summary_data = []
for temp_name in temp_order:
    data = all_data[temp_name]
    summary_data.append({
        'Temperature_K': int(temp_name.replace('K', '')),
        'Temperature_C': int(temp_name.replace('K', '')) - 273,
        'Condition': data['label'].split(' - ')[1] if ' - ' in data['label'] else '',
        'Rg_Angstrom': round(data['rg'], 2)
    })

df_summary = pd.DataFrame(summary_data)
df_summary.to_csv('../results/tables/temperature_summary.csv', index=False)
print("\nSaved: temperature_summary.csv")

print("\n" + "="*60)
print("VISUALIZATION COMPLETE!")
print("="*60)
print("\nGenerated files in results/figures/:")
print("  1. structures_all_temperatures.png  - 2x3 grid of 3D structures")
print("  2. structure_overlay.png            - All structures overlaid")
print("  3. rg_vs_temperature.png            - Radius of gyration trend")
print("  4. structure_2D_projections.png     - XY, XZ, YZ views")
print("  5. residue_displacement.png         - Per-residue structural change")
print("  6. p53_298K.pdb ... p53_373K.pdb    - PDB files for PyMOL/VMD")
