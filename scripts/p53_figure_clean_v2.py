#!/usr/bin/env python3
"""
Professional p53 Figure - Clean Version 2
Clear domain arrangement, no confusion
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle, FancyBboxPatch, Polygon, Ellipse, Wedge, Rectangle
import matplotlib.patheffects as path_effects
import numpy as np
from pathlib import Path

output_dir = Path('results/figures/publication')
output_dir.mkdir(parents=True, exist_ok=True)

print("="*70)
print("CREATING CLEAN p53 FIGURE V2")
print("="*70)

# Create figure
fig, ax = plt.subplots(figsize=(14, 16))
ax.set_xlim(0, 140)
ax.set_ylim(0, 160)
ax.axis('off')

# Background
background = FancyBboxPatch((2, 2), 136, 156,
                            boxstyle="round,pad=0.01,rounding_size=3",
                            facecolor='#F5F7FA', edgecolor='#4A90A4', linewidth=4)
ax.add_patch(background)

# ================================================================
# COLOR SCHEME
# ================================================================
colors = {
    'n_term': '#2ECC71',      # Green
    'tad1': '#E74C3C',        # Red
    'tad2': '#F39C12',        # Orange
    'prd': '#F1C40F',         # Yellow
    'dbd': '#3498DB',         # Blue
    'dbd_dark': '#1A5276',    # Dark blue
    'linker': '#95A5A6',      # Gray
    'od': '#9B59B6',          # Purple
    'ctd': '#7F8C8D',         # Dark gray
    'c_term': '#E74C3C',      # Red
    'hotspot': '#C0392B',     # Dark red
    'zinc': '#F4D03F',        # Gold
    'dna_backbone': '#1ABC9C',# Teal
    'dna_base': '#16A085',    # Dark teal
}

# ================================================================
# N-TERMINUS (Top left corner)
# ================================================================
n_x, n_y = 22, 142

# Green sphere with glow
for i in range(4, 0, -1):
    glow = Circle((n_x, n_y), 7 + i*2, facecolor=colors['n_term'], 
                 alpha=0.1, edgecolor='none')
    ax.add_patch(glow)

n_sphere = Circle((n_x, n_y), 7, facecolor=colors['n_term'], 
                  edgecolor='#1D8348', linewidth=3)
ax.add_patch(n_sphere)

# White highlight
n_highlight = Circle((n_x - 2, n_y + 2), 2.5, facecolor='white', 
                     alpha=0.6, edgecolor='none')
ax.add_patch(n_highlight)

# "N" label inside sphere
ax.text(n_x, n_y, 'N', fontsize=18, fontweight='bold', color='white',
       ha='center', va='center')

# ================================================================
# TAD1 - Right of N-terminus
# ================================================================
# TAD1 box
tad1_box = FancyBboxPatch((n_x + 12, n_y - 5), 20, 10,
                          boxstyle="round,pad=0.02,rounding_size=2",
                          facecolor=colors['tad1'], edgecolor='#922B21', linewidth=2)
ax.add_patch(tad1_box)
ax.text(n_x + 22, n_y, 'TAD1', fontsize=12, fontweight='bold', 
       color='white', ha='center', va='center')

# Residue range above TAD1
ax.text(n_x + 22, n_y + 10, '( 1 - 40 )', fontsize=10, color='#5D6D7E', ha='center')

# ================================================================
# GRADIENT BAR from TAD1 to PRD
# ================================================================
# Starting point (below TAD1)
bar_start_x, bar_start_y = n_x + 22, n_y - 10

# Draw gradient bar going diagonally down-right
n_segments = 25
for i in range(n_segments):
    t = i / n_segments
    
    # Calculate position (curved path)
    x = bar_start_x + t * 30
    y = bar_start_y - t * 40
    
    # Gradient from red to orange to yellow
    if t < 0.5:
        r = 0.91
        g = 0.30 + t * 0.6
        b = 0.24
    else:
        r = 0.91 + (t-0.5) * 0.18
        g = 0.60 + (t-0.5) * 0.4
        b = 0.24 - (t-0.5) * 0.15
    
    circle = Circle((x, y), 4, facecolor=(r, g, b), edgecolor='none', zorder=2)
    ax.add_patch(circle)

# ================================================================
# TAD2 - Along the gradient bar
# ================================================================
tad2_x, tad2_y = bar_start_x + 8, bar_start_y - 18

# TAD2 box (to the left of bar)
tad2_box = FancyBboxPatch((tad2_x - 25, tad2_y - 5), 20, 10,
                          boxstyle="round,pad=0.02,rounding_size=2",
                          facecolor=colors['tad2'], edgecolor='#B7950B', linewidth=2)
ax.add_patch(tad2_box)
ax.text(tad2_x - 15, tad2_y, 'TAD2', fontsize=11, fontweight='bold',
       color='white', ha='center', va='center')

# Residue range
ax.text(tad2_x - 15, tad2_y + 10, '( 41 - 61 )', fontsize=9, color='#5D6D7E', ha='center')

# ================================================================
# PRD - End of gradient bar
# ================================================================
prd_x, prd_y = bar_start_x + 25, bar_start_y - 35

# PRD box (to the left of bar)
prd_box = FancyBboxPatch((prd_x - 28, prd_y - 4), 16, 8,
                         boxstyle="round,pad=0.02,rounding_size=2",
                         facecolor=colors['prd'], edgecolor='#B7950B', linewidth=2)
ax.add_patch(prd_box)
ax.text(prd_x - 20, prd_y, 'PRD', fontsize=10, fontweight='bold',
       color='#1a1a1a', ha='center', va='center')

# Residue range
ax.text(prd_x - 20, prd_y + 9, '( 62 - 93 )', fontsize=9, color='#5D6D7E', ha='center')

# ================================================================
# DBD - LARGE CENTRAL SPHERE
# ================================================================
dbd_x, dbd_y = 70, 75
dbd_r = 30

# Shadow
shadow = Ellipse((dbd_x + 4, dbd_y - 4), dbd_r*2.2, dbd_r*2.2,
                facecolor='#2C3E50', alpha=0.2, edgecolor='none', zorder=1)
ax.add_patch(shadow)

# Main sphere
dbd_outer = Circle((dbd_x, dbd_y), dbd_r, facecolor=colors['dbd'],
                   edgecolor=colors['dbd_dark'], linewidth=5, zorder=5)
ax.add_patch(dbd_outer)

# Inner gradient effect
dbd_inner = Circle((dbd_x - 8, dbd_y + 8), dbd_r * 0.6,
                   facecolor='#5DADE2', alpha=0.4, edgecolor='none', zorder=6)
ax.add_patch(dbd_inner)

# Highlight
highlight = Wedge((dbd_x - 10, dbd_y + 12), dbd_r * 0.4, 30, 150,
                  facecolor='white', alpha=0.35, edgecolor='none', zorder=7)
ax.add_patch(highlight)

# DBD label
ax.text(dbd_x, dbd_y + 3, 'DBD', fontsize=32, fontweight='bold',
       color='white', ha='center', va='center', zorder=8)
ax.text(dbd_x, dbd_y - 10, '( 94 - 292 )', fontsize=12, color='white',
       ha='center', va='center', zorder=8)

# ================================================================
# CANCER HOTSPOTS ON DBD
# ================================================================
hotspot_positions = [
    (dbd_x - 20, dbd_y + 18),
    (dbd_x + 20, dbd_y + 18),
    (dbd_x - 24, dbd_y - 3),
    (dbd_x + 24, dbd_y - 3),
    (dbd_x - 15, dbd_y - 20),
    (dbd_x + 15, dbd_y - 20),
    (dbd_x, dbd_y + 25),
]

for hx, hy in hotspot_positions:
    # Glow
    glow = Circle((hx, hy), 4.5, facecolor=colors['hotspot'], alpha=0.3, 
                 edgecolor='none', zorder=8)
    ax.add_patch(glow)
    # Sphere
    hotspot = Circle((hx, hy), 3, facecolor=colors['hotspot'],
                    edgecolor='#7B241C', linewidth=1.5, zorder=9)
    ax.add_patch(hotspot)

# ================================================================
# ZINC BINDING SITES
# ================================================================
zinc_positions = [
    (dbd_x - 12, dbd_y + 8),
    (dbd_x + 5, dbd_y - 8),
]

for zx, zy in zinc_positions:
    # Glow
    for i in range(3, 0, -1):
        glow = Circle((zx, zy), 4 + i*1.5, facecolor=colors['zinc'],
                     alpha=0.12, edgecolor='none', zorder=8)
        ax.add_patch(glow)
    # Sphere
    zinc = Circle((zx, zy), 4, facecolor=colors['zinc'],
                 edgecolor='#B7950B', linewidth=2, zorder=9)
    ax.add_patch(zinc)

# Zn2+ label
ax.text(dbd_x - 18, dbd_y + 14, 'Zn²⁺', fontsize=10, fontweight='bold',
       color='#9A7D0A', zorder=10)

# ================================================================
# OD - HELIX (Right side of DBD)
# ================================================================
od_start_x, od_start_y = dbd_x + dbd_r + 8, dbd_y + 5

# Draw helix as overlapping ellipses
for i in range(8):
    ell_x = od_start_x + i * 3
    ell_y = od_start_y - i * 5
    
    # Shadow
    shadow_ell = Ellipse((ell_x + 1, ell_y - 1), 10, 16,
                        angle=-25, facecolor='#4A235A', alpha=0.3,
                        edgecolor='none', zorder=3)
    ax.add_patch(shadow_ell)
    
    # Main ellipse
    ell = Ellipse((ell_x, ell_y), 10, 16, angle=-25,
                 facecolor=colors['od'], edgecolor='#6C3483',
                 linewidth=2, zorder=4)
    ax.add_patch(ell)
    
    # Highlight
    ell_high = Ellipse((ell_x - 1.5, ell_y + 3), 4, 7, angle=-25,
                      facecolor='white', alpha=0.35, edgecolor='none', zorder=5)
    ax.add_patch(ell_high)

# OD label box (below helix)
od_label_x, od_label_y = od_start_x + 15, od_start_y - 50

od_box = FancyBboxPatch((od_label_x - 10, od_label_y - 5), 20, 10,
                        boxstyle="round,pad=0.02,rounding_size=2",
                        facecolor=colors['od'], edgecolor='#5B2C6F', linewidth=2)
ax.add_patch(od_box)
ax.text(od_label_x, od_label_y, 'OD', fontsize=12, fontweight='bold',
       color='white', ha='center', va='center')

# Residue range below OD
ax.text(od_label_x, od_label_y - 12, '( 324 - 355 )', fontsize=10, 
       color='#5D6D7E', ha='center')

# ================================================================
# LINKER (Gray bar from OD to CTD)
# ================================================================
linker_start_x = od_label_x + 5
linker_start_y = od_label_y + 10

ax.plot([linker_start_x, linker_start_x + 10], 
       [linker_start_y, linker_start_y + 25],
       color=colors['linker'], linewidth=8, solid_capstyle='round', zorder=2)

# ================================================================
# CTD (Top right area)
# ================================================================
ctd_x, ctd_y = linker_start_x + 8, linker_start_y + 32

# CTD box
ctd_box = FancyBboxPatch((ctd_x - 10, ctd_y - 5), 20, 10,
                         boxstyle="round,pad=0.02,rounding_size=2",
                         facecolor=colors['ctd'], edgecolor='#5D6D7E', linewidth=2)
ax.add_patch(ctd_box)
ax.text(ctd_x, ctd_y, 'CTD', fontsize=11, fontweight='bold',
       color='white', ha='center', va='center')

# Residue range
ax.text(ctd_x, ctd_y + 10, '( 356 - 393 )', fontsize=10, color='#5D6D7E', ha='center')

# ================================================================
# C-TERMINUS (Red sphere, right of CTD)
# ================================================================
c_x, c_y = ctd_x + 18, ctd_y + 5

# Red glow
for i in range(4, 0, -1):
    glow = Circle((c_x, c_y), 6 + i*2, facecolor=colors['c_term'],
                 alpha=0.1, edgecolor='none')
    ax.add_patch(glow)

# Red sphere
c_sphere = Circle((c_x, c_y), 6, facecolor=colors['c_term'],
                 edgecolor='#922B21', linewidth=3)
ax.add_patch(c_sphere)

# Highlight
c_highlight = Circle((c_x - 2, c_y + 2), 2, facecolor='white',
                    alpha=0.6, edgecolor='none')
ax.add_patch(c_highlight)

# "C" label
ax.text(c_x, c_y, 'C', fontsize=16, fontweight='bold', color='white',
       ha='center', va='center')

# ================================================================
# DNA DOUBLE HELIX (Bottom)
# ================================================================
dna_x, dna_y = 70, 18
dna_length = 50

# Arrow from DBD to DNA
ax.annotate('', xy=(dna_x, dna_y + 12), xytext=(dbd_x, dbd_y - dbd_r - 3),
           arrowprops=dict(arrowstyle='->', color='#2C3E50', lw=4))

# Draw DNA helix
n_points = 50
t = np.linspace(0, 5 * np.pi, n_points)

x_positions = dna_x - dna_length/2 + (dna_length/n_points) * np.arange(n_points)
y1 = dna_y + 5 * np.sin(t)
y2 = dna_y + 5 * np.sin(t + np.pi)

# Backbone strands
ax.plot(x_positions, y1, color=colors['dna_backbone'], linewidth=5, zorder=3)
ax.plot(x_positions, y2, color=colors['dna_backbone'], linewidth=5, zorder=3)

# Base pairs and nucleotides
for i in range(0, n_points, 2):
    if abs(np.sin(t[i])) < 0.85:
        # Base pair
        ax.plot([x_positions[i], x_positions[i]], [y1[i], y2[i]], 
               color='#7DCEA0', linewidth=2, zorder=2)
        
        # Nucleotides
        c1 = Circle((x_positions[i], y1[i]), 2, facecolor=colors['dna_backbone'],
                   edgecolor='#148F77', linewidth=1, zorder=4)
        c2 = Circle((x_positions[i], y2[i]), 2, facecolor='#27AE60',
                   edgecolor='#1E8449', linewidth=1, zorder=4)
        ax.add_patch(c1)
        ax.add_patch(c2)

# DNA label
ax.text(dna_x, dna_y - 12, 'Target DNA', fontsize=14, fontweight='bold',
       color='#1E8449', ha='center')

# ================================================================
# LEGEND: Cancer Hotspots (Left side)
# ================================================================
leg1_x, leg1_y = 8, 80
leg1_w, leg1_h = 30, 40

leg1_box = FancyBboxPatch((leg1_x, leg1_y), leg1_w, leg1_h,
                          boxstyle="round,pad=0.02,rounding_size=2",
                          facecolor='white', edgecolor=colors['hotspot'],
                          linewidth=3)
ax.add_patch(leg1_box)

ax.text(leg1_x + leg1_w/2, leg1_y + leg1_h - 5, 'Cancer Hotspots',
       fontsize=12, fontweight='bold', color=colors['hotspot'], ha='center')

hotspots_list = ['R175H', 'R248', 'R273', 'R282', 'G245', 'R249']
for i, name in enumerate(hotspots_list):
    hy = leg1_y + leg1_h - 12 - i * 5
    dot = Circle((leg1_x + 6, hy), 2, facecolor=colors['hotspot'],
                edgecolor='#7B241C', linewidth=1)
    ax.add_patch(dot)
    ax.text(leg1_x + 11, hy, name, fontsize=10, color='#1a1a1a', va='center')

# ================================================================
# LEGEND: Zinc Binding Sites (Left side, below hotspots)
# ================================================================
leg2_x, leg2_y = 8, 55
leg2_w, leg2_h = 30, 18

leg2_box = FancyBboxPatch((leg2_x, leg2_y), leg2_w, leg2_h,
                          boxstyle="round,pad=0.02,rounding_size=2",
                          facecolor='white', edgecolor=colors['zinc'],
                          linewidth=3)
ax.add_patch(leg2_box)

ax.text(leg2_x + leg2_w/2, leg2_y + leg2_h - 5, 'Zinc Binding Sites',
       fontsize=11, fontweight='bold', color='#9A7D0A', ha='center')

zn_dot = Circle((leg2_x + 6, leg2_y + 6), 3, facecolor=colors['zinc'],
               edgecolor='#B7950B', linewidth=1.5)
ax.add_patch(zn_dot)
ax.text(leg2_x + 12, leg2_y + 6, 'Structural Zn²⁺', fontsize=10,
       color='#1a1a1a', va='center')

# ================================================================
# SAVE
# ================================================================
plt.savefig(output_dir / 'Figure1_p53_Clean_v2.png', dpi=300,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure1_p53_Clean_v2.pdf',
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure1_p53_Clean_v2.tiff', dpi=300,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()

print("\n" + "="*70)
print("CLEAN p53 FIGURE V2 CREATED")
print("="*70)
print(f"\nSaved:")
print(f"  ✓ {output_dir}/Figure1_p53_Clean_v2.png")
print(f"  ✓ {output_dir}/Figure1_p53_Clean_v2.pdf")
print(f"  ✓ {output_dir}/Figure1_p53_Clean_v2.tiff")
print("\nLayout:")
print("  • N (green sphere) → TAD1 (red box) → gradient bar")
print("  • TAD2 (orange) and PRD (yellow) along the bar")
print("  • DBD (large blue sphere) in center")
print("  • OD (purple helix) → Linker → CTD (gray) → C (red sphere)")
print("  • DNA double helix at bottom")
print("  • Clear legends on left side")
