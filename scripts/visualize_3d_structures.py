#!/usr/bin/env python3
"""
Create publication-quality 3D protein structure visualizations
showing p53 unfolding at different temperatures
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import MDAnalysis as mda
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Paths
sim_dir = Path("/home/raza/p53_thermal_xai/data/structures/long_simulations")
ref_gro = sim_dir.parent / "npt.gro"
output_dir = Path("/home/raza/p53_thermal_xai/results/figures")
output_dir.mkdir(parents=True, exist_ok=True)

# Temperatures
temperatures = [248, 273, 310, 350, 400, 450, 500]
temp_labels = ['-25°C', '0°C', '37°C', '77°C', '127°C', '177°C', '227°C']

# Hotspot residues (cancer mutation sites)
hotspots = {
    'R175': 175,
    'G245': 245,
    'R248': 248,
    'R249': 249,
    'R273': 273,
    'R282': 282
}

# Color scheme for temperatures (cold to hot)
temp_colors = ['#0000FF', '#00BFFF', '#00FF00', '#FFFF00', '#FFA500', '#FF4500', '#FF0000']

def get_secondary_structure_colors(n_residues):
    """Assign colors based on approximate secondary structure regions"""
    colors = []
    for i in range(n_residues):
        res_num = 91 + i  # Starting residue number
        
        # Beta sheets (approximate regions in p53 DBD)
        if res_num in range(112, 124) or res_num in range(128, 142) or \
           res_num in range(171, 182) or res_num in range(204, 214) or \
           res_num in range(232, 244) or res_num in range(253, 263):
            colors.append('#FFD700')  # Gold for beta sheets
        # Alpha helices
        elif res_num in range(163, 170) or res_num in range(278, 289):
            colors.append('#FF1493')  # Deep pink for helices
        else:
            colors.append('#00CED1')  # Dark cyan for loops
    
    return colors

def extract_final_frame(temp):
    """Extract final frame coordinates from trajectory"""
    xtc_file = sim_dir / f"md_{temp}K.xtc"
    
    if not xtc_file.exists():
        return None, None, None
    
    u = mda.Universe(str(ref_gro), str(xtc_file))
    ca_atoms = u.select_atoms("protein and name CA")
    protein = u.select_atoms("protein")
    
    # Go to last frame
    u.trajectory[-1]
    
    ca_positions = ca_atoms.positions.copy()
    all_positions = protein.positions.copy()
    
    # Get residue numbers
    residue_nums = [atom.resid for atom in ca_atoms]
    
    return ca_positions, all_positions, residue_nums

print("=" * 70)
print("Creating 3D Protein Structure Visualizations")
print("=" * 70)

# ===== Figure 1: Individual structures at each temperature =====
print("\n1. Creating individual temperature structures...")

fig = plt.figure(figsize=(20, 12))

for idx, temp in enumerate(temperatures):
    print(f"   Processing {temp}K...")
    
    ca_pos, all_pos, res_nums = extract_final_frame(temp)
    
    if ca_pos is None:
        continue
    
    ax = fig.add_subplot(2, 4, idx + 1, projection='3d')
    
    # Center the protein
    center = np.mean(ca_pos, axis=0)
    ca_pos_centered = ca_pos - center
    
    # Get secondary structure colors
    ss_colors = get_secondary_structure_colors(len(ca_pos))
    
    # Draw backbone as tube
    for i in range(len(ca_pos_centered) - 1):
        xs = [ca_pos_centered[i, 0], ca_pos_centered[i+1, 0]]
        ys = [ca_pos_centered[i, 1], ca_pos_centered[i+1, 1]]
        zs = [ca_pos_centered[i, 2], ca_pos_centered[i+1, 2]]
        ax.plot(xs, ys, zs, color=ss_colors[i], linewidth=3, alpha=0.8)
    
    # Draw CA atoms as small spheres
    ax.scatter(ca_pos_centered[:, 0], ca_pos_centered[:, 1], ca_pos_centered[:, 2],
               c=ss_colors, s=30, alpha=0.7, edgecolors='none')
    
    # Highlight hotspot residues as large spheres
    for name, res_num in hotspots.items():
        # Find index of this residue
        if res_num in res_nums:
            idx_res = res_nums.index(res_num)
            pos = ca_pos_centered[idx_res]
            ax.scatter([pos[0]], [pos[1]], [pos[2]], 
                      c='red', s=200, alpha=0.9, edgecolors='black', linewidth=1)
    
    # Set labels and title
    ax.set_xlabel('X (Å)', fontsize=8)
    ax.set_ylabel('Y (Å)', fontsize=8)
    ax.set_zlabel('Z (Å)', fontsize=8)
    ax.set_title(f'{temp}K ({temp_labels[idx]})', fontsize=12, fontweight='bold',
                 color=temp_colors[idx])
    
    # Set equal aspect ratio
    max_range = np.max(np.abs(ca_pos_centered)) * 1.1
    ax.set_xlim([-max_range, max_range])
    ax.set_ylim([-max_range, max_range])
    ax.set_zlim([-max_range, max_range])
    
    # Remove grid for cleaner look
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.grid(True, alpha=0.3)

# Add legend in the last subplot
ax = fig.add_subplot(2, 4, 8)
ax.axis('off')
ax.text(0.1, 0.9, 'Legend:', fontsize=14, fontweight='bold', transform=ax.transAxes)
ax.text(0.1, 0.75, '● Beta sheets', fontsize=11, color='#FFD700', transform=ax.transAxes)
ax.text(0.1, 0.60, '● Alpha helices', fontsize=11, color='#FF1493', transform=ax.transAxes)
ax.text(0.1, 0.45, '● Loops/coils', fontsize=11, color='#00CED1', transform=ax.transAxes)
ax.text(0.1, 0.30, '● Hotspot residues', fontsize=11, color='red', transform=ax.transAxes)
ax.text(0.1, 0.10, 'Hotspots: R175, G245, R248,\nR249, R273, R282', fontsize=9, transform=ax.transAxes)

plt.suptitle('p53 DNA-Binding Domain Structure at Different Temperatures', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(output_dir / 'p53_3d_all_temperatures.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print("   Saved: p53_3d_all_temperatures.png")

# ===== Figure 2: Overlay of all temperatures =====
print("\n2. Creating temperature overlay...")

fig = plt.figure(figsize=(14, 12))
ax = fig.add_subplot(111, projection='3d')

# Store first structure for alignment reference
ref_pos = None

for idx, temp in enumerate(temperatures):
    ca_pos, all_pos, res_nums = extract_final_frame(temp)
    
    if ca_pos is None:
        continue
    
    # Center the protein
    center = np.mean(ca_pos, axis=0)
    ca_pos_centered = ca_pos - center
    
    if ref_pos is None:
        ref_pos = ca_pos_centered.copy()
    
    # Draw backbone
    for i in range(len(ca_pos_centered) - 1):
        xs = [ca_pos_centered[i, 0], ca_pos_centered[i+1, 0]]
        ys = [ca_pos_centered[i, 1], ca_pos_centered[i+1, 1]]
        zs = [ca_pos_centered[i, 2], ca_pos_centered[i+1, 2]]
        ax.plot(xs, ys, zs, color=temp_colors[idx], linewidth=2.5, 
                alpha=0.6, label=f'{temp}K' if i == 0 else '')

ax.set_xlabel('X (Å)', fontsize=12)
ax.set_ylabel('Y (Å)', fontsize=12)
ax.set_zlabel('Z (Å)', fontsize=12)
ax.set_title('p53 Structure Overlay: Cold (Blue) → Hot (Red)', fontsize=14, fontweight='bold')
ax.legend(loc='upper left', fontsize=10)

max_range = 25
ax.set_xlim([-max_range, max_range])
ax.set_ylim([-max_range, max_range])
ax.set_zlim([-max_range, max_range])

plt.tight_layout()
plt.savefig(output_dir / 'p53_3d_overlay.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print("   Saved: p53_3d_overlay.png")

# ===== Figure 3: Side-by-side comparison (Cold vs Hot) =====
print("\n3. Creating cold vs hot comparison...")

fig = plt.figure(figsize=(16, 8))

# Cold structure (248K)
ax1 = fig.add_subplot(121, projection='3d')
ca_pos_cold, _, res_nums = extract_final_frame(248)
center = np.mean(ca_pos_cold, axis=0)
ca_pos_cold = ca_pos_cold - center
ss_colors = get_secondary_structure_colors(len(ca_pos_cold))

for i in range(len(ca_pos_cold) - 1):
    xs = [ca_pos_cold[i, 0], ca_pos_cold[i+1, 0]]
    ys = [ca_pos_cold[i, 1], ca_pos_cold[i+1, 1]]
    zs = [ca_pos_cold[i, 2], ca_pos_cold[i+1, 2]]
    ax1.plot(xs, ys, zs, color=ss_colors[i], linewidth=4, alpha=0.9)

ax1.scatter(ca_pos_cold[:, 0], ca_pos_cold[:, 1], ca_pos_cold[:, 2],
           c=ss_colors, s=50, alpha=0.8, edgecolors='none')

# Hotspots
for name, res_num in hotspots.items():
    if res_num in res_nums:
        idx_res = res_nums.index(res_num)
        pos = ca_pos_cold[idx_res]
        ax1.scatter([pos[0]], [pos[1]], [pos[2]], 
                   c='limegreen', s=300, alpha=0.95, edgecolors='black', linewidth=2)
        ax1.text(pos[0]+2, pos[1]+2, pos[2]+2, name, fontsize=8, color='black')

ax1.set_title('248K (-25°C) - STABLE', fontsize=14, fontweight='bold', color='blue')
ax1.set_xlabel('X (Å)')
ax1.set_ylabel('Y (Å)')
ax1.set_zlabel('Z (Å)')
max_range = 25
ax1.set_xlim([-max_range, max_range])
ax1.set_ylim([-max_range, max_range])
ax1.set_zlim([-max_range, max_range])

# Hot structure (500K)
ax2 = fig.add_subplot(122, projection='3d')
ca_pos_hot, _, res_nums = extract_final_frame(500)
center = np.mean(ca_pos_hot, axis=0)
ca_pos_hot = ca_pos_hot - center
ss_colors = get_secondary_structure_colors(len(ca_pos_hot))

for i in range(len(ca_pos_hot) - 1):
    xs = [ca_pos_hot[i, 0], ca_pos_hot[i+1, 0]]
    ys = [ca_pos_hot[i, 1], ca_pos_hot[i+1, 1]]
    zs = [ca_pos_hot[i, 2], ca_pos_hot[i+1, 2]]
    ax2.plot(xs, ys, zs, color=ss_colors[i], linewidth=4, alpha=0.9)

ax2.scatter(ca_pos_hot[:, 0], ca_pos_hot[:, 1], ca_pos_hot[:, 2],
           c=ss_colors, s=50, alpha=0.8, edgecolors='none')

# Hotspots
for name, res_num in hotspots.items():
    if res_num in res_nums:
        idx_res = res_nums.index(res_num)
        pos = ca_pos_hot[idx_res]
        ax2.scatter([pos[0]], [pos[1]], [pos[2]], 
                   c='red', s=300, alpha=0.95, edgecolors='black', linewidth=2)
        ax2.text(pos[0]+2, pos[1]+2, pos[2]+2, name, fontsize=8, color='black')

ax2.set_title('500K (227°C) - DENATURED', fontsize=14, fontweight='bold', color='red')
ax2.set_xlabel('X (Å)')
ax2.set_ylabel('Y (Å)')
ax2.set_zlabel('Z (Å)')
ax2.set_xlim([-max_range, max_range])
ax2.set_ylim([-max_range, max_range])
ax2.set_zlim([-max_range, max_range])

plt.suptitle('p53 Thermal Denaturation: Stable vs Unfolded', fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(output_dir / 'p53_3d_cold_vs_hot.png', dpi=200, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print("   Saved: p53_3d_cold_vs_hot.png")

# ===== Figure 4: Publication-style multi-view =====
print("\n4. Creating publication-style figure...")

fig = plt.figure(figsize=(18, 14))

# Create 3 rows: cold, medium, hot with 3 views each
view_angles = [(30, 45), (30, 135), (90, 0)]  # Different viewing angles
temps_to_show = [248, 350, 500]
temp_titles = ['248K (-25°C)\nStable', '350K (77°C)\nDestabilizing', '500K (227°C)\nDenatured']
row_colors = ['#0066CC', '#FF9900', '#CC0000']

for row, temp in enumerate(temps_to_show):
    ca_pos, _, res_nums = extract_final_frame(temp)
    center = np.mean(ca_pos, axis=0)
    ca_pos = ca_pos - center
    ss_colors = get_secondary_structure_colors(len(ca_pos))
    
    for col, (elev, azim) in enumerate(view_angles):
        ax = fig.add_subplot(3, 3, row * 3 + col + 1, projection='3d')
        
        # Draw backbone with thickness
        for i in range(len(ca_pos) - 1):
            xs = [ca_pos[i, 0], ca_pos[i+1, 0]]
            ys = [ca_pos[i, 1], ca_pos[i+1, 1]]
            zs = [ca_pos[i, 2], ca_pos[i+1, 2]]
            ax.plot(xs, ys, zs, color=ss_colors[i], linewidth=4, alpha=0.9)
        
        # CA atoms
        ax.scatter(ca_pos[:, 0], ca_pos[:, 1], ca_pos[:, 2],
                   c=ss_colors, s=40, alpha=0.7, edgecolors='none')
        
        # Hotspots as spheres
        for name, res_num in hotspots.items():
            if res_num in res_nums:
                idx_res = res_nums.index(res_num)
                pos = ca_pos[idx_res]
                ax.scatter([pos[0]], [pos[1]], [pos[2]], 
                          c='limegreen', s=250, alpha=0.95, edgecolors='darkgreen', linewidth=1.5)
        
        ax.view_init(elev=elev, azim=azim)
        
        max_range = 25
        ax.set_xlim([-max_range, max_range])
        ax.set_ylim([-max_range, max_range])
        ax.set_zlim([-max_range, max_range])
        
        # Remove axis labels for cleaner look
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        
        if col == 0:
            ax.set_ylabel(temp_titles[row], fontsize=12, fontweight='bold', 
                         color=row_colors[row], labelpad=20)
        
        if row == 0:
            view_names = ['Front View', 'Side View', 'Top View']
            ax.set_title(view_names[col], fontsize=11)
        
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.grid(False)

plt.suptitle('p53 DNA-Binding Domain: Temperature-Induced Structural Changes', 
             fontsize=16, fontweight='bold', y=0.98)
plt.tight_layout()
plt.savefig(output_dir / 'p53_3d_publication.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
print("   Saved: p53_3d_publication.png")

# ===== Save PDB files for PyMOL/Chimera =====
print("\n5. Saving PDB files for external visualization...")

pdb_dir = output_dir / "pdb_files"
pdb_dir.mkdir(exist_ok=True)

for temp in temperatures:
    xtc_file = sim_dir / f"md_{temp}K.xtc"
    
    if not xtc_file.exists():
        continue
    
    u = mda.Universe(str(ref_gro), str(xtc_file))
    protein = u.select_atoms("protein")
    
    # Save first frame
    u.trajectory[0]
    protein.write(str(pdb_dir / f"p53_{temp}K_first.pdb"))
    
    # Save last frame
    u.trajectory[-1]
    protein.write(str(pdb_dir / f"p53_{temp}K_last.pdb"))
    
    print(f"   Saved: p53_{temp}K_first.pdb, p53_{temp}K_last.pdb")

print("\n" + "=" * 70)
print("3D Visualization Complete!")
print("=" * 70)
print(f"\nFigures saved to: {output_dir}")
print(f"PDB files saved to: {pdb_dir}")
