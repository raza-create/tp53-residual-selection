#!/usr/bin/env python3
"""
p53 Figure - MODERN PROFESSIONAL DESIGN
Horizontal flow layout - Clean and elegant
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle, FancyBboxPatch, Ellipse, Wedge, Rectangle, FancyArrow
from matplotlib.collections import PatchCollection
import matplotlib.patheffects as path_effects
import numpy as np
from pathlib import Path

output_dir = Path('results/figures/publication')
output_dir.mkdir(parents=True, exist_ok=True)

print("="*70)
print("CREATING MODERN PROFESSIONAL p53 FIGURE")
print("="*70)

# Create figure - Wide format
fig, ax = plt.subplots(figsize=(18, 12))
ax.set_xlim(0, 180)
ax.set_ylim(0, 120)
ax.axis('off')

# Elegant gradient background
for i in range(100):
    y = i * 1.2
    alpha = 0.02
    color = (0.94 - i*0.001, 0.95 - i*0.001, 0.98 - i*0.0005)
    ax.axhspan(y, y + 1.2, facecolor=color, edgecolor='none')

# Frame
frame = FancyBboxPatch((3, 3), 174, 114,
                       boxstyle="round,pad=0.01,rounding_size=2",
                       facecolor='none', edgecolor='#3182CE', linewidth=3)
ax.add_patch(frame)

# ================================================================
# TITLE
# ================================================================
ax.text(90, 112, 'p53 Tumor Suppressor Protein', fontsize=20,
       fontweight='bold', color='#1A365D', ha='center',
       path_effects=[path_effects.withStroke(linewidth=3, foreground='white')])
ax.text(90, 106, '393 Amino Acids | The Guardian of the Genome', fontsize=12,
       color='#4A5568', ha='center', style='italic')

# ================================================================
# MAIN HORIZONTAL PROTEIN BAR (Central element)
# ================================================================
bar_y = 65
bar_height = 20
bar_start = 15
bar_end = 165

# Define domains with exact positions scaled to bar
total_length = 393
scale = (bar_end - bar_start) / total_length

domains = [
    ('N', 1, 1, '#48BB78', True),           # N-terminus (circle)
    ('TAD1', 1, 40, '#F56565', False),       # Red
    ('TAD2', 41, 61, '#ED8936', False),      # Orange
    ('PRD', 62, 93, '#ECC94B', False),       # Yellow
    ('DBD', 94, 292, '#4299E1', False),      # Blue (main)
    ('Linker', 293, 323, '#A0AEC0', False),  # Gray
    ('OD', 324, 355, '#9F7AEA', False),      # Purple
    ('CTD', 356, 393, '#667EEA', False),     # Indigo
    ('C', 393, 393, '#F56565', True),        # C-terminus (circle)
]

# Draw shadow for entire bar
shadow = FancyBboxPatch((bar_start + 2, bar_y - bar_height/2 - 2), 
                        bar_end - bar_start, bar_height,
                        boxstyle="round,pad=0.01,rounding_size=3",
                        facecolor='#2D3748', alpha=0.15, edgecolor='none')
ax.add_patch(shadow)

# Draw each domain
for name, start, end, color, is_circle in domains:
    x_start = bar_start + (start - 1) * scale
    x_end = bar_start + end * scale
    width = x_end - x_start
    
    if is_circle:
        # Terminal circles
        cx = x_start if name == 'N' else x_end
        
        # Glow
        for i in range(3, 0, -1):
            ax.add_patch(Circle((cx, bar_y), 8 + i*2, facecolor=color,
                               alpha=0.1, edgecolor='none', zorder=10))
        
        # Main circle
        ax.add_patch(Circle((cx, bar_y), 8, facecolor=color,
                           edgecolor='#2D3748', linewidth=2, zorder=11))
        
        # Highlight
        ax.add_patch(Circle((cx - 2, bar_y + 2), 3, facecolor='white',
                           alpha=0.5, edgecolor='none', zorder=12))
        
        # Label
        ax.text(cx, bar_y, name, fontsize=14, fontweight='bold', color='white',
               ha='center', va='center', zorder=13)
    else:
        # Domain rectangles with 3D effect
        # Top highlight
        highlight = FancyBboxPatch((x_start, bar_y), width, bar_height/2,
                                   boxstyle="round,pad=0.01,rounding_size=2",
                                   facecolor='white', alpha=0.2, edgecolor='none', zorder=3)
        ax.add_patch(highlight)
        
        # Main rectangle
        rect = FancyBboxPatch((x_start, bar_y - bar_height/2), width, bar_height,
                              boxstyle="round,pad=0.01,rounding_size=2",
                              facecolor=color, edgecolor='#2D3748', linewidth=1.5, zorder=2)
        ax.add_patch(rect)
        
        # Domain label
        if name == 'DBD':
            ax.text(x_start + width/2, bar_y + 2, name, fontsize=16, fontweight='bold',
                   color='white', ha='center', va='center', zorder=4)
            ax.text(x_start + width/2, bar_y - 5, f'({start}-{end})', fontsize=10,
                   color='white', ha='center', va='center', zorder=4)
        elif width > 15:
            ax.text(x_start + width/2, bar_y, name, fontsize=11, fontweight='bold',
                   color='white', ha='center', va='center', zorder=4)

# ================================================================
# RESIDUE SCALE (Below bar)
# ================================================================
scale_y = bar_y - bar_height/2 - 8

ax.plot([bar_start, bar_end], [scale_y, scale_y], color='#4A5568', lw=1.5)

tick_positions = [1, 40, 61, 93, 94, 200, 292, 323, 355, 393]
for pos in tick_positions:
    x = bar_start + (pos - 1) * scale
    ax.plot([x, x], [scale_y - 2, scale_y + 2], color='#4A5568', lw=1.5)
    ax.text(x, scale_y - 5, str(pos), fontsize=8, color='#4A5568', ha='center', va='top')

# ================================================================
# DBD EXPANDED VIEW (Below main bar)
# ================================================================
dbd_detail_y = 28
dbd_x_center = bar_start + (193 - 1) * scale  # Center of DBD

# Connector lines from DBD to detail
ax.plot([dbd_x_center - 30, dbd_x_center - 40], [bar_y - bar_height/2, dbd_detail_y + 18],
       color='#4299E1', lw=2, linestyle='--', alpha=0.5)
ax.plot([dbd_x_center + 30, dbd_x_center + 40], [bar_y - bar_height/2, dbd_detail_y + 18],
       color='#4299E1', lw=2, linestyle='--', alpha=0.5)

# DBD detail box
dbd_detail_width = 100
dbd_detail_height = 35
dbd_detail_x = dbd_x_center - dbd_detail_width/2

detail_box = FancyBboxPatch((dbd_detail_x, dbd_detail_y), dbd_detail_width, dbd_detail_height,
                            boxstyle="round,pad=0.02,rounding_size=3",
                            facecolor='#EBF8FF', edgecolor='#3182CE', linewidth=2)
ax.add_patch(detail_box)

ax.text(dbd_x_center, dbd_detail_y + dbd_detail_height - 5, 
       'DNA-Binding Domain (DBD) - Key Features', fontsize=11, fontweight='bold',
       color='#2C5282', ha='center')

# Key structural features inside detail box
features = [
    ('L2 Loop', '167-195', 'Zn²⁺ binding', '#E53E3E'),
    ('L3 Loop', '237-250', 'DNA contact', '#DD6B20'),
    ('H2 Helix', '278-286', 'Recognition', '#805AD5'),
]

feat_x_start = dbd_detail_x + 8
for i, (name, residues, function, color) in enumerate(features):
    fx = feat_x_start + i * 33
    fy = dbd_detail_y + 12
    
    # Feature box
    ax.add_patch(FancyBboxPatch((fx, fy), 28, 18,
                                boxstyle="round,pad=0.02,rounding_size=1.5",
                                facecolor=color, edgecolor='#2D3748', linewidth=1))
    
    ax.text(fx + 14, fy + 12, name, fontsize=9, fontweight='bold', color='white', ha='center')
    ax.text(fx + 14, fy + 6, f'({residues})', fontsize=7, color='white', ha='center')
    ax.text(fx + 14, fy + 1, function, fontsize=6, color='white', ha='center', style='italic')

# ================================================================
# CANCER HOTSPOTS (Above DBD)
# ================================================================
hotspot_y = bar_y + bar_height/2 + 8

hotspots = [
    (175, 'R175H', 4.8),
    (245, 'G245S', 2.7),
    (248, 'R248Q', 6.4),
    (273, 'R273H', 4.5),
    (282, 'R282W', 1.8),
]

# Hotspot title
ax.text(bar_start + (193 - 1) * scale, hotspot_y + 20, 'Cancer Hotspot Mutations',
       fontsize=11, fontweight='bold', color='#C53030', ha='center')

for pos, name, freq in hotspots:
    hx = bar_start + (pos - 1) * scale
    
    # Connector line
    ax.plot([hx, hx], [bar_y + bar_height/2, hotspot_y], color='#C53030', lw=1.5, alpha=0.6)
    
    # Hotspot circle
    ax.add_patch(Circle((hx, hotspot_y + 8), 4, facecolor='#FC8181',
                       edgecolor='#C53030', linewidth=1.5))
    
    # Label
    ax.text(hx, hotspot_y + 16, f'{name}', fontsize=8, fontweight='bold',
           color='#C53030', ha='center', rotation=45)

# ================================================================
# ZINC COORDINATION (Left side)
# ================================================================
zn_x, zn_y = 20, 28

# Zinc box
zn_box = FancyBboxPatch((zn_x - 5, zn_y - 5), 40, 35,
                        boxstyle="round,pad=0.02,rounding_size=2",
                        facecolor='#FFFFF0', edgecolor='#D69E2E', linewidth=2)
ax.add_patch(zn_box)

ax.text(zn_x + 15, zn_y + 25, 'Zinc Coordination', fontsize=10, fontweight='bold',
       color='#975A16', ha='center')

# Zinc sphere
ax.add_patch(Circle((zn_x + 15, zn_y + 10), 6, facecolor='#F6E05E',
                   edgecolor='#B7791F', linewidth=2))
ax.text(zn_x + 15, zn_y + 10, 'Zn²⁺', fontsize=8, fontweight='bold',
       color='#744210', ha='center', va='center')

# Coordinating residues
residues = ['C176', 'H179', 'C238', 'C242']
positions = [(zn_x + 5, zn_y + 2), (zn_x + 25, zn_y + 2), 
             (zn_x + 5, zn_y + 18), (zn_x + 25, zn_y + 18)]

for (rx, ry), res in zip(positions, residues):
    ax.text(rx, ry, res, fontsize=7, color='#975A16', ha='center')
    # Line to zinc
    ax.plot([rx + 3, zn_x + 15], [ry + 2, zn_y + 10], color='#D69E2E', lw=1, alpha=0.5)

# ================================================================
# DNA BINDING (Right side)
# ================================================================
dna_x, dna_y = 155, 28

# DNA box
dna_box = FancyBboxPatch((dna_x - 15, dna_y - 5), 35, 35,
                         boxstyle="round,pad=0.02,rounding_size=2",
                         facecolor='#E6FFFA', edgecolor='#319795', linewidth=2)
ax.add_patch(dna_box)

ax.text(dna_x + 2, dna_y + 25, 'DNA Binding', fontsize=10, fontweight='bold',
       color='#234E52', ha='center')

# Simple DNA helix
for i in range(5):
    offset = i * 5
    # Strand 1
    ax.add_patch(Circle((dna_x - 5 + offset*0.3, dna_y + 5 + np.sin(i)*3), 2,
                       facecolor='#38B2AC', edgecolor='#285E61', linewidth=1))
    # Strand 2  
    ax.add_patch(Circle((dna_x + 5 + offset*0.3, dna_y + 5 - np.sin(i)*3), 2,
                       facecolor='#4FD1C5', edgecolor='#285E61', linewidth=1))
    # Connection
    if i % 2 == 0:
        ax.plot([dna_x - 5 + offset*0.3, dna_x + 5 + offset*0.3], 
               [dna_y + 5 + np.sin(i)*3, dna_y + 5 - np.sin(i)*3],
               color='#81E6D9', lw=1.5)

# ================================================================
# FUNCTION LABELS (Below scale)
# ================================================================
func_y = scale_y - 15

functions = [
    (bar_start + 20 * scale, 'Transcription\nActivation', '#F56565'),
    (bar_start + 77 * scale, 'Proline-Rich\nSignaling', '#ECC94B'),
    (bar_start + 193 * scale, 'DNA Binding\n& Zn²⁺ Coordination', '#4299E1'),
    (bar_start + 340 * scale, 'Tetramerization', '#9F7AEA'),
]

for fx, text, color in functions:
    ax.text(fx, func_y, text, fontsize=8, color=color, ha='center', va='top',
           fontweight='bold')

# ================================================================
# LEGEND (Bottom right)
# ================================================================
legend_x, legend_y = 145, 8

ax.add_patch(Circle((legend_x, legend_y + 3), 3, facecolor='#FC8181', edgecolor='#C53030', lw=1.5))
ax.text(legend_x + 6, legend_y + 3, 'Mutation hotspot', fontsize=8, va='center')

ax.add_patch(Circle((legend_x, legend_y - 4), 3, facecolor='#F6E05E', edgecolor='#B7791F', lw=1.5))
ax.text(legend_x + 6, legend_y - 4, 'Zinc ion', fontsize=8, va='center')

# ================================================================
# SAVE
# ================================================================
plt.tight_layout()
plt.savefig(output_dir / 'Figure1_p53_Modern.png', dpi=300,
            bbox_inches='tight', facecolor='white')
plt.savefig(output_dir / 'Figure1_p53_Modern.pdf', bbox_inches='tight')
plt.savefig(output_dir / 'Figure1_p53_Modern.tiff', dpi=300,
            bbox_inches='tight', facecolor='white')
plt.close()

print("\n" + "="*70)
print("MODERN PROFESSIONAL FIGURE CREATED")
print("="*70)
print(f"\nSaved:")
print(f"  ✓ {output_dir}/Figure1_p53_Modern.png")
print(f"  ✓ {output_dir}/Figure1_p53_Modern.pdf") 
print(f"  ✓ {output_dir}/Figure1_p53_Modern.tiff")
print("\nDESIGN FEATURES:")
print("  ✓ Clean horizontal layout")
print("  ✓ Color-coded domains with smooth transitions")
print("  ✓ N and C termini as distinct spheres")
print("  ✓ DBD expanded view with key features")
print("  ✓ Cancer hotspots clearly marked above")
print("  ✓ Zinc coordination diagram")
print("  ✓ DNA binding illustration")
print("  ✓ Residue scale")
print("  ✓ Functional annotations")
