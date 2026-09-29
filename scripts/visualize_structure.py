#!/usr/bin/env python3
"""
Visualize protein structure and create plots
"""

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

# Load trajectory
print("Loading trajectory...")
u = mda.Universe(
    '../data/structures/npt.gro',
    '../data/structures/md_310K_short.xtc'
)

protein = u.select_atoms('protein')
ca_atoms = u.select_atoms('protein and name CA')  # Alpha carbons only

print(f"CA atoms: {len(ca_atoms)}")

# Create output directory
import os
os.makedirs('../results/figures', exist_ok=True)

#############################################
# FIGURE 1: 3D Structure at First Frame
#############################################
print("\nCreating 3D structure plot (first frame)...")

u.trajectory[0]  # Go to first frame
positions = ca_atoms.positions

fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Color by residue number
colors = np.arange(len(ca_atoms))

# Plot backbone as connected line
ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], 
        'b-', alpha=0.3, linewidth=1)

# Plot CA atoms as spheres
scatter = ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2],
                     c=colors, cmap='rainbow', s=50, alpha=0.8)

# Mark hotspot positions (approximate - residues 175, 245, 248, 249, 273, 282)
# Adjust indices based on starting residue 91
hotspots = [175-91, 245-91, 248-91, 249-91, 273-91, 282-91]
hotspot_names = ['R175', 'G245', 'R248', 'R249', 'R273', 'R282']

for i, (hs, name) in enumerate(zip(hotspots, hotspot_names)):
    if 0 <= hs < len(positions):
        ax.scatter(positions[hs, 0], positions[hs, 1], positions[hs, 2],
                   c='red', s=200, marker='*', edgecolors='black', linewidths=2)
        ax.text(positions[hs, 0], positions[hs, 1], positions[hs, 2], 
                f'  {name}', fontsize=10, fontweight='bold')

ax.set_xlabel('X (Å)', fontsize=12)
ax.set_ylabel('Y (Å)', fontsize=12)
ax.set_zlabel('Z (Å)', fontsize=12)
ax.set_title('p53 DNA-Binding Domain Structure (310K)\nHotspot Positions Marked with Stars', fontsize=14)

plt.colorbar(scatter, label='Residue Number', shrink=0.6)
plt.tight_layout()
plt.savefig('../results/figures/structure_3D_310K.png', dpi=150, bbox_inches='tight')
print("Saved: ../results/figures/structure_3D_310K.png")
plt.close()

#############################################
# FIGURE 2: Compare First vs Last Frame
#############################################
print("\nCreating structure comparison plot...")

fig = plt.figure(figsize=(16, 7))

# First frame
ax1 = fig.add_subplot(121, projection='3d')
u.trajectory[0]
pos_first = ca_atoms.positions.copy()
ax1.plot(pos_first[:, 0], pos_first[:, 1], pos_first[:, 2], 'b-', alpha=0.5, linewidth=1)
ax1.scatter(pos_first[:, 0], pos_first[:, 1], pos_first[:, 2], c='blue', s=30, alpha=0.7)
ax1.set_title('First Frame (t = 0 ps)', fontsize=14)
ax1.set_xlabel('X (Å)')
ax1.set_ylabel('Y (Å)')
ax1.set_zlabel('Z (Å)')

# Last frame
ax2 = fig.add_subplot(122, projection='3d')
u.trajectory[-1]
pos_last = ca_atoms.positions.copy()
ax2.plot(pos_last[:, 0], pos_last[:, 1], pos_last[:, 2], 'r-', alpha=0.5, linewidth=1)
ax2.scatter(pos_last[:, 0], pos_last[:, 1], pos_last[:, 2], c='red', s=30, alpha=0.7)
ax2.set_title(f'Last Frame (t = {u.trajectory[-1].time:.0f} ps)', fontsize=14)
ax2.set_xlabel('X (Å)')
ax2.set_ylabel('Y (Å)')
ax2.set_zlabel('Z (Å)')

plt.suptitle('p53 Structure: Start vs End of Simulation', fontsize=16)
plt.tight_layout()
plt.savefig('../results/figures/structure_comparison.png', dpi=150, bbox_inches='tight')
print("Saved: ../results/figures/structure_comparison.png")
plt.close()

#############################################
# FIGURE 3: RMSD and Rg Over Time
#############################################
print("\nCreating RMSD and Rg plot...")

# Load feature data
df = pd.read_csv('../data/features/features_310K.csv')

fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# RMSD plot
axes[0].plot(df['time_ps'], df['rmsd_nm'], 'b-', linewidth=2, marker='o')
axes[0].set_ylabel('RMSD (nm)', fontsize=12)
axes[0].set_title('Structural Stability Over Time (310K)', fontsize=14)
axes[0].axhline(y=df['rmsd_nm'].mean(), color='r', linestyle='--', label=f'Mean: {df["rmsd_nm"].mean():.3f} nm')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Rg plot
axes[1].plot(df['time_ps'], df['rg_nm'], 'g-', linewidth=2, marker='s')
axes[1].set_xlabel('Time (ps)', fontsize=12)
axes[1].set_ylabel('Radius of Gyration (nm)', fontsize=12)
axes[1].axhline(y=df['rg_nm'].mean(), color='r', linestyle='--', label=f'Mean: {df["rg_nm"].mean():.3f} nm')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('../results/figures/rmsd_rg_plot.png', dpi=150, bbox_inches='tight')
print("Saved: ../results/figures/rmsd_rg_plot.png")
plt.close()

#############################################
# FIGURE 4: Contact Map
#############################################
print("\nCreating contact map...")

u.trajectory[-1]  # Use last frame
positions = ca_atoms.positions

# Calculate distance matrix
n_atoms = len(positions)
dist_matrix = np.zeros((n_atoms, n_atoms))

for i in range(n_atoms):
    for j in range(n_atoms):
        dist_matrix[i, j] = np.linalg.norm(positions[i] - positions[j])

# Create contact map (distance < 8 Angstrom)
contact_map = dist_matrix < 8.0

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Distance matrix
im1 = axes[0].imshow(dist_matrix, cmap='viridis', origin='lower')
axes[0].set_xlabel('Residue Index', fontsize=12)
axes[0].set_ylabel('Residue Index', fontsize=12)
axes[0].set_title('Distance Matrix (Å)', fontsize=14)
plt.colorbar(im1, ax=axes[0], label='Distance (Å)')

# Contact map
im2 = axes[1].imshow(contact_map, cmap='Blues', origin='lower')
axes[1].set_xlabel('Residue Index', fontsize=12)
axes[1].set_ylabel('Residue Index', fontsize=12)
axes[1].set_title('Contact Map (< 8 Å)', fontsize=14)

# Mark hotspots
for hs in hotspots:
    if 0 <= hs < n_atoms:
        axes[1].axhline(y=hs, color='red', alpha=0.3, linewidth=1)
        axes[1].axvline(x=hs, color='red', alpha=0.3, linewidth=1)

plt.tight_layout()
plt.savefig('../results/figures/contact_map.png', dpi=150, bbox_inches='tight')
print("Saved: ../results/figures/contact_map.png")
plt.close()

#############################################
# FIGURE 5: Per-Residue RMSF (Flexibility)
#############################################
print("\nCalculating per-residue flexibility (RMSF)...")

from MDAnalysis.analysis import rms

# Calculate RMSF
rmsf_values = []
for i, atom in enumerate(ca_atoms):
    positions_over_time = []
    for ts in u.trajectory:
        positions_over_time.append(ca_atoms[i].position.copy())
    positions_over_time = np.array(positions_over_time)
    mean_pos = positions_over_time.mean(axis=0)
    rmsf = np.sqrt(((positions_over_time - mean_pos) ** 2).sum(axis=1).mean())
    rmsf_values.append(rmsf)

rmsf_values = np.array(rmsf_values)
residue_numbers = np.arange(91, 91 + len(rmsf_values))  # Starting from residue 91

fig, ax = plt.subplots(figsize=(14, 5))
ax.plot(residue_numbers, rmsf_values, 'b-', linewidth=1.5)
ax.fill_between(residue_numbers, 0, rmsf_values, alpha=0.3)

# Mark hotspots
hotspot_residues = [175, 245, 248, 249, 273, 282]
for hr in hotspot_residues:
    if 91 <= hr <= 289:
        idx = hr - 91
        ax.axvline(x=hr, color='red', linestyle='--', alpha=0.5)
        ax.scatter(hr, rmsf_values[idx], c='red', s=100, zorder=5, marker='*')
        ax.text(hr, rmsf_values[idx] + 0.2, f'{hr}', ha='center', fontsize=9)

ax.set_xlabel('Residue Number', fontsize=12)
ax.set_ylabel('RMSF (Å)', fontsize=12)
ax.set_title('Per-Residue Flexibility (RMSF) at 310K\nRed stars = Hotspot positions', fontsize=14)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('../results/figures/rmsf_per_residue.png', dpi=150, bbox_inches='tight')
print("Saved: ../results/figures/rmsf_per_residue.png")
plt.close()

#############################################
# Save PDB snapshots for external viewing
#############################################
print("\nSaving PDB snapshots for 3D visualization...")

# First frame
u.trajectory[0]
protein.write('../results/figures/p53_frame_first.pdb')
print("Saved: ../results/figures/p53_frame_first.pdb")

# Last frame
u.trajectory[-1]
protein.write('../results/figures/p53_frame_last.pdb')
print("Saved: ../results/figures/p53_frame_last.pdb")

print("\n" + "="*50)
print("VISUALIZATION COMPLETE!")
print("="*50)
print("\nGenerated files:")
print("  1. structure_3D_310K.png      - 3D structure with hotspots")
print("  2. structure_comparison.png   - First vs Last frame")
print("  3. rmsd_rg_plot.png          - RMSD and Rg over time")
print("  4. contact_map.png           - Distance and contact maps")
print("  5. rmsf_per_residue.png      - Per-residue flexibility")
print("  6. p53_frame_first.pdb       - PDB for 3D viewing")
print("  7. p53_frame_last.pdb        - PDB for 3D viewing")
print("\nTo view PDB files in 3D, use:")
print("  - PyMOL: pymol p53_frame_first.pdb")
print("  - VMD: vmd p53_frame_first.pdb")
print("  - Chimera: chimera p53_frame_first.pdb")
print("  - Online: https://www.rcsb.org/3d-view")
