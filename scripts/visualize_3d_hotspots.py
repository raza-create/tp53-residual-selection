#!/usr/bin/env python3
"""
Create 3D p53 Structure Visualizations with Cancer Hotspots
Publication-quality figures for high-impact journal
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

base_dir = Path("/home/raza/p53_thermal_xai")
struct_dir = base_dir / "data/structures"
output_dir = base_dir / "results/figures"

# Cancer hotspot residues (original numbering)
HOTSPOTS = {
    'R175': {'pos': 175, 'type': 'structural', 'color': 'red'},
    'G245': {'pos': 245, 'type': 'structural', 'color': 'orange'},
    'R248': {'pos': 248, 'type': 'contact', 'color': 'blue'},
    'R249': {'pos': 249, 'type': 'contact', 'color': 'blue'},
    'R273': {'pos': 273, 'type': 'contact', 'color': 'blue'},
    'R282': {'pos': 282, 'type': 'structural', 'color': 'red'},
}

def parse_pdb(pdb_file):
    """Parse PDB file and extract CA coordinates"""
    coords = []
    residues = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and ' CA ' in line:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                resname = line[17:20].strip()
                resnum = int(line[22:26])
                coords.append([x, y, z])
                residues.append({'name': resname, 'num': resnum})
    return np.array(coords), residues

# Find PDB files
pdb_files = {}
for temp in [248, 310, 450, 500]:
    # Check multiple possible locations
    candidates = [
        struct_dir / f"long_simulations/md_{temp}K.gro",
        struct_dir / f"p53_{temp}K.pdb",
        base_dir / f"results/figures/pdb_files/p53_{temp}K.pdb",
    ]
    for c in candidates:
        if c.exists():
            pdb_files[temp] = c
            break

# Use reference structure
ref_pdb = struct_dir / "2XWR.pdb"
if not ref_pdb.exists():
    ref_pdb = struct_dir / "p53_chainA.pdb"

print(f"Using reference: {ref_pdb}")
coords, residues = parse_pdb(ref_pdb)
print(f"Loaded {len(coords)} CA atoms")

# Create Figure 1: 3D structure with hotspots
fig = plt.figure(figsize=(14, 10))
ax = fig.add_subplot(111, projection='3d')

# Plot backbone
ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], 'gray', alpha=0.5, linewidth=1)

# Color by position (N to C terminal)
colors = plt.cm.viridis(np.linspace(0, 1, len(coords)))
ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=colors, s=20, alpha=0.6)

# Highlight hotspots
hotspot_coords = []
hotspot_labels = []
for name, info in HOTSPOTS.items():
    for i, res in enumerate(residues):
        if res['num'] == info['pos']:
            ax.scatter(coords[i, 0], coords[i, 1], coords[i, 2], 
                      c=info['color'], s=300, marker='o', edgecolors='black', linewidth=2)
            ax.text(coords[i, 0]+2, coords[i, 1]+2, coords[i, 2]+2, 
                   name, fontsize=12, fontweight='bold', color=info['color'])
            hotspot_coords.append(coords[i])
            hotspot_labels.append(name)
            break

ax.set_xlabel('X (Å)', fontsize=12)
ax.set_ylabel('Y (Å)', fontsize=12)
ax.set_zlabel('Z (Å)', fontsize=12)
ax.set_title('p53 DNA-Binding Domain - Cancer Hotspot Residues\n(Red=Structural, Blue=DNA Contact)', 
             fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig(output_dir / 'p53_3d_hotspots_labeled.png', dpi=300, bbox_inches='tight')
print(f"Saved: p53_3d_hotspots_labeled.png")
plt.close()

# Create Figure 2: 2D contact map with hotspots
fig, ax = plt.subplots(figsize=(10, 10))

# Calculate distance matrix
n = len(coords)
dist_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        dist_matrix[i, j] = np.linalg.norm(coords[i] - coords[j])

# Plot contact map (contacts < 8 Angstrom)
contact_map = dist_matrix < 8.0
ax.imshow(contact_map, cmap='Blues', origin='lower')

# Mark hotspots
hotspot_indices = []
for name, info in HOTSPOTS.items():
    for i, res in enumerate(residues):
        if res['num'] == info['pos']:
            hotspot_indices.append(i)
            ax.axhline(y=i, color='red', alpha=0.5, linewidth=1)
            ax.axvline(x=i, color='red', alpha=0.5, linewidth=1)
            ax.text(n+2, i, name, fontsize=10, color='red', fontweight='bold', va='center')
            break

ax.set_xlabel('Residue Index', fontsize=12)
ax.set_ylabel('Residue Index', fontsize=12)
ax.set_title('p53 Contact Map with Cancer Hotspot Positions', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(output_dir / 'p53_contact_map_hotspots.png', dpi=300, bbox_inches='tight')
print(f"Saved: p53_contact_map_hotspots.png")
plt.close()

# Create Figure 3: Mutation effect mapped on structure (1D projection)
fig, ax = plt.subplots(figsize=(14, 6))

# Residue positions along sequence
res_nums = [r['num'] for r in residues]
ax.fill_between(range(len(res_nums)), 0, 1, alpha=0.2, color='gray', label='p53 DBD')

# Mark hotspots with mutation effects
mutation_effects = {
    'R175H': -0.84,
    'G245S': 0.06,
    'R248W': -1.67,
    'R248Q': -0.62,
    'R273H': -1.09,
    'R282W': -1.60,
    'Y220C': 0.34,
}

for name, llr in mutation_effects.items():
    res_num = int(name[1:-1])
    for i, res in enumerate(residues):
        if res['num'] == res_num:
            color = 'red' if llr < -1 else ('orange' if llr < 0 else 'green')
            height = abs(llr) / 2
            ax.bar(i, height, color=color, width=3, alpha=0.8)
            ax.text(i, height + 0.1, name, rotation=90, fontsize=10, 
                   fontweight='bold', ha='center', va='bottom')
            break

ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='Deleterious threshold')
ax.set_xlabel('Residue Position', fontsize=12)
ax.set_ylabel('|Mutation Effect| (ESM-2 LLR)', fontsize=12)
ax.set_title('Cancer Mutation Effects Mapped on p53 Structure', fontsize=14, fontweight='bold')
ax.legend(loc='upper right')
ax.set_xlim(0, len(residues))
plt.tight_layout()
plt.savefig(output_dir / 'p53_mutation_effect_map.png', dpi=300, bbox_inches='tight')
print(f"Saved: p53_mutation_effect_map.png")
plt.close()

print("\n=== ALL 3D VISUALIZATIONS CREATED ===")
