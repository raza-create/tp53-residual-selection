#!/usr/bin/env python3
"""
Professional p53 Figure - Similar style to reference but different colors
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle, FancyBboxPatch, Polygon, Ellipse, Wedge, Arc
from matplotlib.collections import PatchCollection
import matplotlib.patheffects as path_effects
import numpy as np
from pathlib import Path

output_dir = Path('results/figures/publication')
output_dir.mkdir(parents=True, exist_ok=True)

print("="*70)
print("CREATING PROFESSIONAL p53 FIGURE - NEW STYLE")
print("="*70)

# Create figure
fig, ax = plt.subplots(figsize=(14, 16))
ax.set_xlim(0, 140)
ax.set_ylim(0, 160)
ax.axis('off')

# Background gradient effect (light lavender)
background = FancyBboxPatch((2, 2), 136, 156,
                            boxstyle="round,pad=0.01,rounding_size=3",
                            facecolor='#F0F4F8', edgecolor='#5E60CE', linewidth=4)
ax.add_patch(background)

# ================================================================
# COLOR SCHEME (Different from reference - Teal/Coral/Gold)
# ================================================================
colors = {
    'n_term': '#06D6A0',      # Mint green
    'tad1': '#EF476F',        # Coral pink
    'tad2': '#FFD166',        # Golden yellow
    'prd': '#F4A261',         # Sandy orange
    'dbd': '#118AB2',         # Ocean blue
    'dbd_dark': '#073B4C',    # Dark teal
    'linker': '#8D99AE',      # Gray
    'od': '#9D4EDD',          # Purple
    'ctd': '#7B8CDE',         # Periwinkle
    'c_term': '#EF476F',      # Coral pink
    'hotspot': '#E63946',     # Red
    'zinc': '#FFD700',        # Gold
    'dna1': '#2EC4B6',        # Teal
    'dna2': '#E71D36',        # Red
}

# ================================================================
# N-TERMINUS (Top left - Mint green sphere with glow)
# ================================================================
n_x, n_y = 25, 140

# Glow effect
for i in range(5, 0, -1):
    glow = Circle((n_x, n_y), 6 + i*1.5, facecolor=colors['n_term'], 
                 alpha=0.08, edgecolor='none')
    ax.add_patch(glow)

# Main sphere
n_sphere = Circle((n_x, n_y), 6, facecolor=colors['n_term'], 
                  edgecolor='#048A6E', linewidth=2.5)
ax.add_patch(n_sphere)

# Highlight
n_highlight = Circle((n_x - 2, n_y + 2), 2, facecolor='white', 
                     alpha=0.6, edgecolor='none')
ax.add_patch(n_highlight)

# N label
ax.text(n_x, n_y, 'N', fontsize=16, fontweight='bold', color='white',
       ha='center', va='center')

# Residue range
ax.text(n_x + 10, n_y + 3, '( 1 - 40 )', fontsize=10, color='#495057')

# TAD1 label box
tad1_box = FancyBboxPatch((n_x + 8, n_y - 8), 18, 8,
                          boxstyle="round,pad=0.02,rounding_size=1",
                          facecolor=colors['tad1'], edgecolor='none')
ax.add_patch(tad1_box)
ax.text(n_x + 17, n_y - 4, 'TAD1', fontsize=10, fontweight='bold', 
       color='white', ha='center', va='center')

# ================================================================
# GRADIENT BAR (N-term to DBD)
# ================================================================
# Draw curved gradient bar from TAD1 to PRD to DBD
bar_points_1 = [
    (n_x + 3, n_y - 5),
    (n_x + 10, n_y - 20),
    (n_x + 20, n_y - 35),
    (n_x + 30, n_y - 45),
]

# Draw gradient segments
for i in range(len(bar_points_1) - 1):
    x1, y1 = bar_points_1[i]
    x2, y2 = bar_points_1[i + 1]
    
    # Gradient color
    t = i / (len(bar_points_1) - 1)
    r = (1-t) * 0.94 + t * 0.96  # coral to orange
    g = (1-t) * 0.28 + t * 0.63
    b = (1-t) * 0.44 + t * 0.38
    
    ax.plot([x1, x2], [y1, y2], color=(r, g, b), linewidth=10, 
           solid_capstyle='round', zorder=2)
    # Shadow
    ax.plot([x1+1, x2+1], [y1-1, y2-1], color='#2C3E50', linewidth=10,
           solid_capstyle='round', alpha=0.2, zorder=1)

# TAD2 label
ax.text(n_x + 22, n_y - 22, '( 41 - 61 )', fontsize=9, color='#495057')
tad2_box = FancyBboxPatch((n_x + 15, n_y - 32), 18, 8,
                          boxstyle="round,pad=0.02,rounding_size=1",
                          facecolor=colors['tad2'], edgecolor='none')
ax.add_patch(tad2_box)
ax.text(n_x + 24, n_y - 28, 'TAD2', fontsize=10, fontweight='bold',
       color='#1a1a1a', ha='center', va='center')

# PRD label
ax.text(n_x + 32, n_y - 40, '( 62 - 93 )', fontsize=9, color='#495057')
prd_box = FancyBboxPatch((n_x + 28, n_y - 50), 14, 7,
                         boxstyle="round,pad=0.02,rounding_size=1",
                         facecolor=colors['prd'], edgecolor='none')
ax.add_patch(prd_box)
ax.text(n_x + 35, n_y - 46.5, 'PRD', fontsize=9, fontweight='bold',
       color='white', ha='center', va='center')

# ================================================================
# DBD - LARGE CENTRAL SPHERE
# ================================================================
dbd_x, dbd_y = 70, 80
dbd_r = 28

# Shadow
shadow = Ellipse((dbd_x + 3, dbd_y - 3), dbd_r*2.1, dbd_r*2.1,
                facecolor='#2C3E50', alpha=0.25, edgecolor='none', zorder=1)
ax.add_patch(shadow)

# Main sphere with gradient effect (dark edge)
dbd_outer = Circle((dbd_x, dbd_y), dbd_r, facecolor=colors['dbd'],
                   edgecolor=colors['dbd_dark'], linewidth=4, zorder=5)
ax.add_patch(dbd_outer)

# Inner lighter area
dbd_inner = Circle((dbd_x - 5, dbd_y + 5), dbd_r * 0.7, 
                   facecolor='#1A9DB8', alpha=0.5, edgecolor='none', zorder=6)
ax.add_patch(dbd_inner)

# Highlight arc
highlight_wedge = Wedge((dbd_x - 8, dbd_y + 10), dbd_r * 0.5, 20, 160,
                        facecolor='white', alpha=0.3, edgecolor='none', zorder=7)
ax.add_patch(highlight_wedge)

# DBD label
ax.text(dbd_x, dbd_y + 5, 'DBD', fontsize=28, fontweight='bold', 
       color='white', ha='center', va='center', zorder=8)
ax.text(dbd_x, dbd_y - 8, '( 94 - 292 )', fontsize=12, color='white',
       ha='center', va='center', zorder=8)

# ================================================================
# CANCER HOTSPOTS ON DBD (Red spheres)
# ================================================================
hotspot_positions = [
    (dbd_x - 18, dbd_y + 20),   # R175H
    (dbd_x + 18, dbd_y + 20),   # R248
    (dbd_x - 22, dbd_y),        # G245
    (dbd_x + 22, dbd_y),        # R273
    (dbd_x - 15, dbd_y - 18),   # R282
    (dbd_x + 15, dbd_y - 18),   # R249
    (dbd_x, dbd_y + 24),        # extra
]

for hx, hy in hotspot_positions:
    # Glow
    glow = Circle((hx, hy), 4, facecolor=colors['hotspot'], alpha=0.3, 
                 edgecolor='none', zorder=8)
    ax.add_patch(glow)
    # Sphere
    hotspot = Circle((hx, hy), 2.5, facecolor=colors['hotspot'],
                    edgecolor='#9D0208', linewidth=1, zorder=9)
    ax.add_patch(hotspot)

# ================================================================
# ZINC BINDING SITES (Gold spheres with glow)
# ================================================================
zinc_positions = [
    (dbd_x - 10, dbd_y + 10),
    (dbd_x + 8, dbd_y - 5),
]

for zx, zy in zinc_positions:
    # Glow
    for i in range(3, 0, -1):
        glow = Circle((zx, zy), 3 + i*1.2, facecolor=colors['zinc'],
                     alpha=0.15, edgecolor='none', zorder=8)
        ax.add_patch(glow)
    # Sphere
    zinc = Circle((zx, zy), 3, facecolor=colors['zinc'],
                 edgecolor='#B8860B', linewidth=1.5, zorder=9)
    ax.add_patch(zinc)

# Zn2+ label
ax.text(dbd_x - 15, dbd_y + 15, 'Zn²⁺', fontsize=9, fontweight='bold',
       color='#B8860B', zorder=10)

# ================================================================
# LINKER AND OD (Right side)
# ================================================================
# Linker (gray curved line)
linker_points = [(dbd_x + 25, dbd_y - 5), (dbd_x + 35, dbd_y - 15),
                (dbd_x + 40, dbd_y - 25), (dbd_x + 38, dbd_y - 35)]

for i in range(len(linker_points) - 1):
    x1, y1 = linker_points[i]
    x2, y2 = linker_points[i + 1]
    ax.plot([x1, x2], [y1, y2], color=colors['linker'], linewidth=8,
           solid_capstyle='round', zorder=3)

# OD - Helix representation (Purple ellipses)
od_x, od_y = 115, 55
helix_angle = -30

for i in range(7):
    offset = i * 5
    ell_x = od_x - offset * 0.3
    ell_y = od_y - offset
    
    # Shadow
    shadow_ell = Ellipse((ell_x + 1, ell_y - 1), 8, 14,
                        angle=helix_angle, facecolor='#2C3E50', 
                        alpha=0.2, edgecolor='none', zorder=3)
    ax.add_patch(shadow_ell)
    
    # Main ellipse
    ell = Ellipse((ell_x, ell_y), 8, 14, angle=helix_angle,
                 facecolor=colors['od'], edgecolor='#7B2CBF',
                 linewidth=1.5, alpha=0.9, zorder=4)
    ax.add_patch(ell)
    
    # Highlight
    ell_high = Ellipse((ell_x - 1, ell_y + 2), 3, 6, angle=helix_angle,
                      facecolor='white', alpha=0.3, edgecolor='none', zorder=5)
    ax.add_patch(ell_high)

# OD label
od_label_box = FancyBboxPatch((od_x + 5, od_y - 25), 12, 8,
                              boxstyle="round,pad=0.02,rounding_size=1",
                              facecolor=colors['od'], edgecolor='none')
ax.add_patch(od_label_box)
ax.text(od_x + 11, od_y - 21, 'OD', fontsize=10, fontweight='bold',
       color='white', ha='center', va='center')
ax.text(od_x + 11, od_y - 32, '( 324 - 355 )', fontsize=9, color='#495057',
       ha='center')

# ================================================================
# CTD AND C-TERMINUS (Top right)
# ================================================================
# Linker to CTD
ax.plot([od_x - 5, od_x + 5], [od_y + 20, od_y + 35], color=colors['linker'],
       linewidth=6, solid_capstyle='round', zorder=2)

# CTD label
ctd_box = FancyBboxPatch((od_x - 2, od_y + 38), 14, 7,
                         boxstyle="round,pad=0.02,rounding_size=1",
                         facecolor=colors['ctd'], edgecolor='none')
ax.add_patch(ctd_box)
ax.text(od_x + 5, od_y + 41.5, 'CTD', fontsize=9, fontweight='bold',
       color='white', ha='center', va='center')
ax.text(od_x + 5, od_y + 50, '( 356 - 393 )', fontsize=9, color='#495057',
       ha='center')

# C-terminus sphere
c_x, c_y = od_x + 15, od_y + 55

# Glow
for i in range(5, 0, -1):
    glow = Circle((c_x, c_y), 5 + i*1.2, facecolor=colors['c_term'],
                 alpha=0.08, edgecolor='none')
    ax.add_patch(glow)

c_sphere = Circle((c_x, c_y), 5, facecolor=colors['c_term'],
                 edgecolor='#C1121F', linewidth=2)
ax.add_patch(c_sphere)
c_highlight = Circle((c_x - 1.5, c_y + 1.5), 1.8, facecolor='white',
                    alpha=0.6, edgecolor='none')
ax.add_patch(c_highlight)
ax.text(c_x, c_y, 'C', fontsize=14, fontweight='bold', color='white',
       ha='center', va='center')

# ================================================================
# DNA DOUBLE HELIX (Bottom)
# ================================================================
dna_x, dna_y = 70, 22
dna_length = 45

# Arrow from DBD to DNA
ax.annotate('', xy=(dna_x, dna_y + 15), xytext=(dbd_x, dbd_y - dbd_r - 2),
           arrowprops=dict(arrowstyle='->', color='#2C3E50', lw=3))

# Draw DNA helix
n_points = 40
t = np.linspace(0, 4 * np.pi, n_points)

# Strand 1
x1 = dna_x - dna_length/2 + (dna_length/n_points) * np.arange(n_points)
y1 = dna_y + 5 * np.sin(t)

# Strand 2
x2 = x1
y2 = dna_y + 5 * np.sin(t + np.pi)

# Draw strands
ax.plot(x1, y1, color=colors['dna1'], linewidth=4, zorder=3)
ax.plot(x2, y2, color=colors['dna1'], linewidth=4, zorder=3)

# Base pairs and nucleotides
for i in range(0, n_points, 2):
    if abs(np.sin(t[i])) < 0.8:
        # Base pair line
        ax.plot([x1[i], x2[i]], [y1[i], y2[i]], color='#6C757D', 
               linewidth=2, zorder=2)
        
        # Nucleotide circles
        c1 = Circle((x1[i], y1[i]), 1.5, facecolor=colors['dna1'],
                   edgecolor='#1D9A8C', linewidth=1, zorder=4)
        c2 = Circle((x2[i], y2[i]), 1.5, facecolor=colors['dna2'],
                   edgecolor='#B71C33', linewidth=1, zorder=4)
        ax.add_patch(c1)
        ax.add_patch(c2)

# DNA label
ax.text(dna_x, dna_y - 12, 'Target DNA', fontsize=14, fontweight='bold',
       color='#2C3E50', ha='center')

# ================================================================
# LEGEND BOXES
# ================================================================
# Cancer Hotspots legend (left)
legend1_x, legend1_y = 8, 85
legend1_w, legend1_h = 28, 35

leg1_box = FancyBboxPatch((legend1_x, legend1_y), legend1_w, legend1_h,
                          boxstyle="round,pad=0.02,rounding_size=2",
                          facecolor='white', edgecolor=colors['hotspot'],
                          linewidth=2.5)
ax.add_patch(leg1_box)

ax.text(legend1_x + legend1_w/2, legend1_y + legend1_h - 4,
       'Cancer Hotspots', fontsize=11, fontweight='bold',
       color=colors['hotspot'], ha='center')

hotspot_names = ['R175H', 'R248', 'R273', 'R282', 'G245', 'R249']
for i, name in enumerate(hotspot_names):
    hy = legend1_y + legend1_h - 10 - i * 4.5
    dot = Circle((legend1_x + 5, hy), 1.5, facecolor=colors['hotspot'],
                edgecolor='#9D0208', linewidth=0.8)
    ax.add_patch(dot)
    ax.text(legend1_x + 9, hy, name, fontsize=9, color='#1a1a1a', va='center')

# Zinc Binding Sites legend
legend2_x, legend2_y = 8, 58
legend2_w, legend2_h = 28, 20

leg2_box = FancyBboxPatch((legend2_x, legend2_y), legend2_w, legend2_h,
                          boxstyle="round,pad=0.02,rounding_size=2",
                          facecolor='white', edgecolor=colors['zinc'],
                          linewidth=2.5)
ax.add_patch(leg2_box)

ax.text(legend2_x + legend2_w/2, legend2_y + legend2_h - 4,
       'Zinc Binding Sites', fontsize=10, fontweight='bold',
       color='#B8860B', ha='center')

zn_dot = Circle((legend2_x + 5, legend2_y + 7), 2.5, facecolor=colors['zinc'],
               edgecolor='#B8860B', linewidth=1)
ax.add_patch(zn_dot)
ax.text(legend2_x + 10, legend2_y + 7, 'Structural Zn²⁺', fontsize=9,
       color='#1a1a1a', va='center')

# ================================================================
# TITLE (Optional - can be removed for journal)
# ================================================================
ax.text(70, 153, 'Human p53 Tumor Suppressor', fontsize=16, 
       fontweight='bold', color='#2C3E50', ha='center')
ax.text(70, 148, 'Domain Architecture (393 amino acids)', fontsize=11,
       color='#6C757D', ha='center', style='italic')

# ================================================================
# SAVE
# ================================================================
plt.savefig(output_dir / 'Figure1_p53_NewStyle.png', dpi=300,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure1_p53_NewStyle.pdf',
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure1_p53_NewStyle.tiff', dpi=300,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()

print("\n" + "="*70)
print("NEW STYLE p53 FIGURE CREATED")
print("="*70)
print(f"\nSaved:")
print(f"  ✓ {output_dir}/Figure1_p53_NewStyle.png")
print(f"  ✓ {output_dir}/Figure1_p53_NewStyle.pdf")
print(f"  ✓ {output_dir}/Figure1_p53_NewStyle.tiff")
print("\nColor scheme: Teal/Coral/Gold/Purple")
print("Style: Similar layout, different colors")
