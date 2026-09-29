#!/usr/bin/env python3
"""
Professional p53 Figure - Unique Clean Design
Horizontal layout with proper text alignment
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Circle, Ellipse, Polygon, Rectangle, Arc
from matplotlib.lines import Line2D
import matplotlib.patheffects as path_effects
import numpy as np
from pathlib import Path as FilePath

output_dir = FilePath('results/figures/publication')
output_dir.mkdir(parents=True, exist_ok=True)

print("="*70)
print("CREATING CLEAN PROFESSIONAL p53 FIGURE")
print("="*70)

# Figure setup
fig = plt.figure(figsize=(16, 12))
ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])
ax.set_xlim(0, 160)
ax.set_ylim(0, 120)
ax.axis('off')
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

# ================================================================
# COLOR SCHEME
# ================================================================
colors = {
    'TAD1': '#E63946',      # Red
    'TAD2': '#F4A261',      # Orange
    'PRD': '#E9C46A',       # Yellow
    'DBD': '#2A9D8F',       # Teal
    'Linker': '#8D99AE',    # Gray
    'OD': '#6A4C93',        # Purple
    'CTD': '#457B9D',       # Blue-gray
    'DNA': '#1D3557',       # Dark blue
    'hotspot': '#E63946',   # Red
    'zinc': '#FFD700',      # Gold
}

# ================================================================
# TITLE SECTION (Top)
# ================================================================
ax.text(80, 115, 'p53 Tumor Suppressor Protein Architecture', 
       fontsize=20, fontweight='bold', ha='center', color='#1D3557')
ax.text(80, 110, 'Full-length: 393 amino acids | Gene: TP53 | UniProt: P04637',
       fontsize=11, ha='center', color='#6B7280', style='italic')

# ================================================================
# MAIN PROTEIN STRUCTURE (Horizontal bar with 3D domains)
# ================================================================
y_protein = 70
bar_height = 18

# Domain definitions with positions scaled to figure
domains = [
    # (name, start_x, width, color, residues, label_position)
    ('TAD1', 10, 12, colors['TAD1'], '1-40', 'above'),
    ('TAD2', 22, 7, colors['TAD2'], '41-61', 'below'),
    ('PRD', 29, 10, colors['PRD'], '62-93', 'above'),
    ('DBD', 39, 55, colors['DBD'], '94-292', 'center'),
    ('Linker', 94, 9, colors['Linker'], '293-323', 'below'),
    ('OD', 103, 10, colors['OD'], '324-355', 'above'),
    ('CTD', 113, 12, colors['CTD'], '356-393', 'below'),
]

def draw_3d_rounded_box(ax, x, y, w, h, color, depth=3):
    """Draw a 3D rounded box with shadow and highlights"""
    from matplotlib.colors import to_rgb, to_hex
    
    rgb = to_rgb(color)
    dark = tuple(max(0, c - 0.2) for c in rgb)
    light = tuple(min(1, c + 0.15) for c in rgb)
    
    # Shadow
    shadow = FancyBboxPatch((x + depth*0.7, y - depth*0.7), w, h,
                            boxstyle="round,pad=0.02,rounding_size=1.5",
                            facecolor='#2C3E50', alpha=0.25, edgecolor='none', zorder=1)
    ax.add_patch(shadow)
    
    # Right side (dark)
    side_r = Polygon([(x+w, y), (x+w+depth, y-depth), (x+w+depth, y+h-depth), (x+w, y+h)],
                    facecolor=dark, edgecolor='none', zorder=2)
    ax.add_patch(side_r)
    
    # Top side (light)
    side_t = Polygon([(x, y+h), (x+depth, y+h-depth), (x+w+depth, y+h-depth), (x+w, y+h)],
                    facecolor=light, edgecolor='none', zorder=2)
    ax.add_patch(side_t)
    
    # Main face
    main = FancyBboxPatch((x, y), w, h,
                          boxstyle="round,pad=0.02,rounding_size=1.5",
                          facecolor=color, edgecolor='#1a1a1a', linewidth=1.5, zorder=3)
    ax.add_patch(main)
    
    # Highlight gradient (top part)
    highlight = FancyBboxPatch((x+1, y+h*0.65), w-2, h*0.3,
                               boxstyle="round,pad=0.01,rounding_size=1",
                               facecolor='white', alpha=0.2, edgecolor='none', zorder=4)
    ax.add_patch(highlight)

# Draw all domains
for name, x, w, color, residues, label_pos in domains:
    draw_3d_rounded_box(ax, x, y_protein, w, bar_height, color)
    
    # Domain name on the box
    if name == 'DBD':
        ax.text(x + w/2, y_protein + bar_height/2 + 2, name, ha='center', va='center',
               fontsize=16, fontweight='bold', color='white', zorder=5)
        ax.text(x + w/2, y_protein + bar_height/2 - 4, f'({residues})', ha='center', va='center',
               fontsize=10, color='white', zorder=5)
    elif w > 8:
        ax.text(x + w/2, y_protein + bar_height/2, name, ha='center', va='center',
               fontsize=11, fontweight='bold', color='white', zorder=5)

# ================================================================
# DOMAIN LABELS (Above and Below - No Overlap)
# ================================================================
label_above_y = y_protein + bar_height + 8
label_below_y = y_protein - 12

for name, x, w, color, residues, label_pos in domains:
    if name == 'DBD':
        continue  # Already labeled on the box
    
    cx = x + w/2
    
    if label_pos == 'above':
        # Connector line
        ax.plot([cx, cx], [y_protein + bar_height, label_above_y - 3], 
               color=color, lw=1.5, linestyle='-', alpha=0.7)
        # Label box
        bbox = dict(boxstyle='round,pad=0.3', facecolor=color, edgecolor='none', alpha=0.9)
        ax.text(cx, label_above_y, f'{name}\n({residues})', ha='center', va='bottom',
               fontsize=9, fontweight='bold', color='white', bbox=bbox)
    else:
        # Connector line
        ax.plot([cx, cx], [y_protein, label_below_y + 6], 
               color=color, lw=1.5, linestyle='-', alpha=0.7)
        # Label box
        bbox = dict(boxstyle='round,pad=0.3', facecolor=color, edgecolor='none', alpha=0.9)
        ax.text(cx, label_below_y, f'{name}\n({residues})', ha='center', va='top',
               fontsize=9, fontweight='bold', color='white', bbox=bbox)

# ================================================================
# N and C TERMINUS MARKERS
# ================================================================
# N-terminus (left)
n_circle = Circle((6, y_protein + bar_height/2), 4, facecolor='#2ECC71', 
                  edgecolor='#1a1a1a', linewidth=2, zorder=6)
ax.add_patch(n_circle)
ax.text(6, y_protein + bar_height/2, 'N', ha='center', va='center',
       fontsize=12, fontweight='bold', color='white', zorder=7)
ax.plot([10, 10], [y_protein, y_protein + bar_height], color='#1a1a1a', lw=2)

# C-terminus (right)
c_circle = Circle((129, y_protein + bar_height/2), 4, facecolor='#E74C3C',
                  edgecolor='#1a1a1a', linewidth=2, zorder=6)
ax.add_patch(c_circle)
ax.text(129, y_protein + bar_height/2, 'C', ha='center', va='center',
       fontsize=12, fontweight='bold', color='white', zorder=7)
ax.plot([125, 125], [y_protein, y_protein + bar_height], color='#1a1a1a', lw=2)

# ================================================================
# DBD DETAIL SECTION (Below main bar)
# ================================================================
detail_y = 30
detail_x = 39
detail_w = 55
detail_h = 25

# DBD detail box
detail_box = FancyBboxPatch((detail_x, detail_y), detail_w, detail_h,
                            boxstyle="round,pad=0.02,rounding_size=2",
                            facecolor='#E8F4F3', edgecolor=colors['DBD'], 
                            linewidth=2, zorder=1)
ax.add_patch(detail_box)

ax.text(detail_x + detail_w/2, detail_y + detail_h - 3, 
       'DBD Structural Features (94-292)', ha='center', va='top',
       fontsize=12, fontweight='bold', color=colors['DBD'])

# Connector from main DBD to detail
ax.annotate('', xy=(detail_x + detail_w/2, detail_y + detail_h),
           xytext=(detail_x + detail_w/2, y_protein - 2),
           arrowprops=dict(arrowstyle='->', color=colors['DBD'], lw=2,
                          connectionstyle='arc3,rad=0'))

# Key features in DBD detail box
features = [
    ('L2 Loop', '167-195', 'Zn²⁺ binding', colors['hotspot']),
    ('L3 Loop', '237-250', 'DNA contact', colors['hotspot']),
    ('H2 Helix', '278-286', 'DNA recognition', colors['OD']),
]

feat_start_y = detail_y + detail_h - 10
for i, (feat_name, residues, function, color) in enumerate(features):
    fx = detail_x + 5 + i * 18
    
    # Feature box
    feat_box = FancyBboxPatch((fx, feat_start_y - 8), 16, 12,
                              boxstyle="round,pad=0.02,rounding_size=0.8",
                              facecolor=color, edgecolor='#1a1a1a', 
                              linewidth=1, alpha=0.85, zorder=2)
    ax.add_patch(feat_box)
    
    ax.text(fx + 8, feat_start_y - 2, feat_name, ha='center', va='center',
           fontsize=9, fontweight='bold', color='white', zorder=3)
    ax.text(fx + 8, feat_start_y - 5.5, f'({residues})', ha='center', va='center',
           fontsize=7, color='white', zorder=3)

# ================================================================
# CANCER HOTSPOTS PANEL (Left side)
# ================================================================
hotspot_x = 5
hotspot_y = 5
hotspot_w = 28
hotspot_h = 38

# Panel box
hotspot_box = FancyBboxPatch((hotspot_x, hotspot_y), hotspot_w, hotspot_h,
                             boxstyle="round,pad=0.02,rounding_size=1.5",
                             facecolor='white', edgecolor=colors['hotspot'],
                             linewidth=2.5, zorder=1)
ax.add_patch(hotspot_box)

# Title
ax.text(hotspot_x + hotspot_w/2, hotspot_y + hotspot_h - 3, 
       'Cancer Mutation Hotspots', ha='center', va='top',
       fontsize=11, fontweight='bold', color=colors['hotspot'])

# Hotspot list
hotspots = [
    ('R175H', '4.8%', 'Structural'),
    ('G245S', '2.7%', 'Structural'),
    ('R248Q/W', '6.4%', 'DNA Contact'),
    ('R249S', '2.5%', 'Structural'),
    ('R273H/C', '4.5%', 'DNA Contact'),
    ('R282W', '1.8%', 'Structural'),
]

for i, (mut, freq, mut_type) in enumerate(hotspots):
    hy = hotspot_y + hotspot_h - 9 - i * 5
    
    # Colored dot
    dot_color = colors['hotspot'] if mut_type == 'Structural' else '#3498DB'
    dot = Circle((hotspot_x + 4, hy), 1.5, facecolor=dot_color, 
                edgecolor='#1a1a1a', linewidth=0.8, zorder=2)
    ax.add_patch(dot)
    
    # Mutation name and frequency
    ax.text(hotspot_x + 7, hy, f'{mut}', ha='left', va='center',
           fontsize=9, fontweight='bold', color='#1a1a1a')
    ax.text(hotspot_x + 18, hy, f'({freq})', ha='left', va='center',
           fontsize=8, color='#6B7280')

# Legend for mutation types
ax.plot([hotspot_x + 3], [hotspot_y + 3], 'o', color=colors['hotspot'], markersize=6)
ax.text(hotspot_x + 5.5, hotspot_y + 3, 'Structural', fontsize=7, va='center', color='#1a1a1a')
ax.plot([hotspot_x + 15], [hotspot_y + 3], 'o', color='#3498DB', markersize=6)
ax.text(hotspot_x + 17.5, hotspot_y + 3, 'DNA Contact', fontsize=7, va='center', color='#1a1a1a')

# ================================================================
# ZINC BINDING PANEL (Right side)
# ================================================================
zinc_x = 127
zinc_y = 5
zinc_w = 28
zinc_h = 25

# Panel box
zinc_box = FancyBboxPatch((zinc_x, zinc_y), zinc_w, zinc_h,
                          boxstyle="round,pad=0.02,rounding_size=1.5",
                          facecolor='white', edgecolor=colors['zinc'],
                          linewidth=2.5, zorder=1)
ax.add_patch(zinc_box)

# Title
ax.text(zinc_x + zinc_w/2, zinc_y + zinc_h - 3,
       'Zinc Coordination', ha='center', va='top',
       fontsize=11, fontweight='bold', color='#B8860B')

# Zinc info
zn_circle = Circle((zinc_x + 7, zinc_y + zinc_h - 11), 3, facecolor=colors['zinc'],
                   edgecolor='#B8860B', linewidth=1.5, zorder=2)
ax.add_patch(zn_circle)
ax.text(zinc_x + 7, zinc_y + zinc_h - 11, 'Zn²⁺', ha='center', va='center',
       fontsize=8, fontweight='bold', color='#1a1a1a', zorder=3)

# Coordinating residues
coord_residues = ['Cys176', 'His179', 'Cys238', 'Cys242']
for i, res in enumerate(coord_residues):
    ax.text(zinc_x + 14, zinc_y + zinc_h - 8 - i * 4, f'• {res}',
           fontsize=9, color='#1a1a1a', va='center')

# ================================================================
# FUNCTION ANNOTATIONS (Bottom)
# ================================================================
func_y = 8
functions = [
    (16, 'Transcriptional\nActivation', colors['TAD1']),
    (34, 'Proline-Rich\nSignaling', colors['PRD']),
    (66, 'Sequence-Specific\nDNA Binding', colors['DBD']),
    (108, 'Tetramer\nFormation', colors['OD']),
]

for fx, text, color in functions:
    bbox = dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor=color, linewidth=1.5)
    ax.text(fx, func_y, text, ha='center', va='center', fontsize=8,
           color='#1a1a1a', bbox=bbox)

# ================================================================
# SCALE BAR
# ================================================================
scale_y = y_protein - 25
ax.plot([10, 125], [scale_y, scale_y], color='#1a1a1a', lw=2)

# Tick marks for key residue positions
scale_ticks = [
    (10, '1'), (22, '41'), (29, '62'), (39, '94'),
    (94, '292'), (103, '324'), (113, '356'), (125, '393')
]

for sx, label in scale_ticks:
    ax.plot([sx, sx], [scale_y - 1.5, scale_y + 1.5], color='#1a1a1a', lw=1.5)
    ax.text(sx, scale_y - 4, label, ha='center', va='top', fontsize=8, color='#4B5563')

ax.text(67.5, scale_y - 8, 'Residue Position', ha='center', fontsize=9, 
       color='#4B5563', style='italic')

# ================================================================
# FRAME
# ================================================================
frame = FancyBboxPatch((1, 1), 158, 118,
                       boxstyle="round,pad=0.01,rounding_size=1",
                       facecolor='none', edgecolor='#CBD5E1', linewidth=3)
ax.add_patch(frame)

# ================================================================
# SAVE
# ================================================================
plt.savefig(output_dir / 'Figure1_p53_Clean_Professional.png', dpi=300, 
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure1_p53_Clean_Professional.pdf', 
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure1_p53_Clean_Professional.tiff', dpi=300,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()

print("\n" + "="*70)
print("CLEAN PROFESSIONAL FIGURE CREATED")
print("="*70)
print(f"\nSaved to {output_dir}:")
print("  ✓ Figure1_p53_Clean_Professional.png")
print("  ✓ Figure1_p53_Clean_Professional.pdf")
print("  ✓ Figure1_p53_Clean_Professional.tiff")
print("\nFigure Features:")
print("  ✓ Horizontal 3D domain layout (1-393)")
print("  ✓ All domains clearly labeled (no overlap)")
print("  ✓ N/C terminus markers")
print("  ✓ DBD detail panel with L2, L3, H2")
print("  ✓ Cancer hotspots panel with frequencies")
print("  ✓ Zinc coordination panel")
print("  ✓ Functional annotations")
print("  ✓ Residue scale bar")
