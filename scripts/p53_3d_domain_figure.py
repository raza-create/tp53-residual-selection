#!/usr/bin/env python3
"""
Professional 3D-style Full-Length p53 Domain Figure (1-393)
Single panel, high-impact journal quality
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, Polygon, PathPatch
from matplotlib.path import Path
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as path_effects
import numpy as np
from pathlib import Path as FilePath

output_dir = FilePath('results/figures/publication')
output_dir.mkdir(parents=True, exist_ok=True)

# Publication settings
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 11,
    'axes.linewidth': 0,
    'figure.dpi': 300,
})

print("="*70)
print("CREATING 3D-STYLE FULL-LENGTH p53 FIGURE")
print("="*70)

def draw_3d_domain(ax, x, y, width, height, color, label, label_below=True, 
                   depth=8, highlight=False):
    """Draw a 3D-looking domain box with gradient and shadow"""
    
    # Convert color to RGB if string
    if isinstance(color, str):
        from matplotlib.colors import to_rgb
        base_color = to_rgb(color)
    else:
        base_color = color
    
    # Create darker and lighter versions
    dark_color = tuple(max(0, c - 0.25) for c in base_color)
    light_color = tuple(min(1, c + 0.15) for c in base_color)
    
    # Shadow (back face)
    shadow = FancyBboxPatch((x + depth/2, y - depth/2), width, height,
                            boxstyle="round,pad=0.02,rounding_size=3",
                            facecolor='#2C3E50', edgecolor='none', alpha=0.4,
                            zorder=1)
    ax.add_patch(shadow)
    
    # Side face (right)
    side_points = [
        (x + width, y),
        (x + width + depth, y - depth),
        (x + width + depth, y + height - depth),
        (x + width, y + height),
    ]
    side = Polygon(side_points, facecolor=dark_color, edgecolor='#1a1a1a', 
                   linewidth=1, zorder=2)
    ax.add_patch(side)
    
    # Top face
    top_points = [
        (x, y + height),
        (x + depth, y + height - depth),
        (x + width + depth, y + height - depth),
        (x + width, y + height),
    ]
    top = Polygon(top_points, facecolor=light_color, edgecolor='#1a1a1a',
                  linewidth=1, zorder=2)
    ax.add_patch(top)
    
    # Front face (main)
    front = FancyBboxPatch((x, y), width, height,
                           boxstyle="round,pad=0.02,rounding_size=3",
                           facecolor=color, edgecolor='#1a1a1a', linewidth=2,
                           zorder=3)
    ax.add_patch(front)
    
    # Gradient overlay for 3D effect
    gradient = FancyBboxPatch((x, y + height*0.6), width, height*0.4,
                              boxstyle="round,pad=0.02,rounding_size=3",
                              facecolor='white', edgecolor='none', alpha=0.15,
                              zorder=4)
    ax.add_patch(gradient)
    
    # Highlight border if specified
    if highlight:
        highlight_box = FancyBboxPatch((x - 4, y - 4), width + depth + 8, height + 8,
                                       boxstyle="round,pad=0.02,rounding_size=5",
                                       facecolor='none', edgecolor='#E74C3C', 
                                       linewidth=4, linestyle='--', zorder=0)
        ax.add_patch(highlight_box)
    
    return x + width/2, y + height/2  # Return center for label


# Create figure
fig, ax = plt.subplots(figsize=(18, 10))
ax.set_xlim(-20, 480)
ax.set_ylim(-60, 120)
ax.axis('off')
ax.set_aspect('equal')

# Title with shadow effect
title = ax.text(230, 105, 'Full-Length Human p53 Tumor Suppressor Protein', 
               fontsize=18, fontweight='bold', ha='center', va='bottom',
               color='#1a1a1a')
title.set_path_effects([path_effects.withStroke(linewidth=3, foreground='white')])

subtitle = ax.text(230, 95, '393 Amino Acids | UniProt: P04637 | Gene: TP53', 
                  fontsize=12, ha='center', va='bottom', color='#5D6D7E', style='italic')

# ================================================================
# DOMAIN DEFINITIONS (1-393)
# ================================================================
y_main = 35
bar_height = 35
depth = 6

domains = [
    # (name, start, end, color, full_name)
    ('TAD1', 1, 40, '#E74C3C', 'Transactivation\nDomain 1'),
    ('TAD2', 41, 61, '#EC7063', 'TAD2'),
    ('PRR', 62, 93, '#F5B041', 'Proline-Rich\nRegion'),
    ('DBD', 94, 292, '#2ECC71', 'DNA-Binding Domain'),
    ('Linker', 293, 323, '#85C1E9', 'Linker'),
    ('TD', 324, 355, '#3498DB', 'Tetramerization\nDomain'),
    ('CTD', 356, 393, '#9B59B6', 'C-Terminal\nDomain'),
]

# Draw domains
domain_centers = {}
for name, start, end, color, full_name in domains:
    width = end - start + 1
    highlight = (name == 'DBD')  # Highlight DBD
    
    cx, cy = draw_3d_domain(ax, start, y_main, width, bar_height, color, name,
                           highlight=highlight, depth=depth)
    domain_centers[name] = (cx, cy)
    
    # Domain label on the bar
    if name == 'DBD':
        ax.text(cx, cy + 5, 'DNA-Binding Domain', ha='center', va='center',
               fontsize=13, fontweight='bold', color='white', zorder=5)
        ax.text(cx, cy - 8, '(94-292)', ha='center', va='center',
               fontsize=11, color='white', zorder=5)
    elif name in ['TAD1', 'PRR', 'TD']:
        ax.text(cx, cy, name, ha='center', va='center',
               fontsize=11, fontweight='bold', color='white', zorder=5)
    else:
        ax.text(cx, cy, name, ha='center', va='center',
               fontsize=10, fontweight='bold', color='white', zorder=5)

# "This Study" annotation for DBD
ax.annotate('', xy=(193, y_main + bar_height + 15), xytext=(193, y_main + bar_height + 35),
           arrowprops=dict(arrowstyle='->', color='#E74C3C', lw=3))
ax.text(193, y_main + bar_height + 40, '★ This Study', ha='center', va='bottom',
       fontsize=14, fontweight='bold', color='#E74C3C')

# ================================================================
# RESIDUE NUMBER SCALE
# ================================================================
scale_y = y_main - 15

# Main scale line with 3D effect
ax.plot([1, 393 + depth], [scale_y, scale_y], color='#2C3E50', lw=3, zorder=1)
ax.plot([1, 393 + depth], [scale_y + 1, scale_y + 1], color='#5D6D7E', lw=1, zorder=2)

# Tick marks and numbers
tick_positions = [1, 40, 61, 93, 94, 150, 200, 250, 292, 323, 355, 393]
for pos in tick_positions:
    ax.plot([pos, pos], [scale_y - 4, scale_y + 4], color='#2C3E50', lw=2)
    ax.text(pos, scale_y - 8, str(pos), ha='center', va='top', fontsize=9,
           fontweight='bold', color='#2C3E50')

# ================================================================
# DOMAIN ANNOTATIONS (Below)
# ================================================================
annotation_y = -35

annotations = [
    (20, 'Transcriptional\nActivation', '#E74C3C'),
    (77, 'SH3 Binding\nMotifs', '#F5B041'),
    (150, 'Sequence-Specific\nDNA Binding', '#2ECC71'),
    (230, 'Zn²⁺\nCoordination', '#2ECC71'),
    (340, 'Oligo-\nmerization', '#3498DB'),
    (375, 'Regulatory\nSignals', '#9B59B6'),
]

for x, text, color in annotations:
    # Connector line
    ax.plot([x, x], [y_main - 2, annotation_y + 20], color=color, lw=1.5, 
           linestyle=':', alpha=0.7)
    
    # Text box
    bbox = dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=color,
                linewidth=2, alpha=0.95)
    ax.text(x, annotation_y, text, ha='center', va='top', fontsize=9,
           fontweight='bold', color='#2C3E50', bbox=bbox)

# ================================================================
# KEY STRUCTURAL FEATURES IN DBD
# ================================================================
# Mark L2, L3 loops and key residues within DBD
feature_y = y_main + bar_height + 5

# L2 loop marker
l2_start, l2_end = 167, 195
ax.annotate('', xy=(l2_start, feature_y), xytext=(l2_end, feature_y),
           arrowprops=dict(arrowstyle='<->', color='#E74C3C', lw=2))
ax.text((l2_start + l2_end)/2, feature_y + 5, 'L2 Loop\n(Zn²⁺)', ha='center', va='bottom',
       fontsize=8, fontweight='bold', color='#E74C3C')

# L3 loop marker  
l3_start, l3_end = 237, 250
ax.annotate('', xy=(l3_start, feature_y), xytext=(l3_end, feature_y),
           arrowprops=dict(arrowstyle='<->', color='#C0392B', lw=2))
ax.text((l3_start + l3_end)/2, feature_y + 5, 'L3 Loop\n(DNA)', ha='center', va='bottom',
       fontsize=8, fontweight='bold', color='#C0392B')

# H2 helix marker
h2_start, h2_end = 278, 286
ax.annotate('', xy=(h2_start, feature_y), xytext=(h2_end, feature_y),
           arrowprops=dict(arrowstyle='<->', color='#1ABC9C', lw=2))
ax.text((h2_start + h2_end)/2, feature_y + 5, 'H2', ha='center', va='bottom',
       fontsize=8, fontweight='bold', color='#1ABC9C')

# ================================================================
# CANCER HOTSPOT MUTATIONS
# ================================================================
hotspots = [
    (175, 'R175H', 4.8, 'Structural'),
    (245, 'G245S', 2.7, 'Structural'),
    (248, 'R248Q', 6.4, 'Contact'),
    (273, 'R273H', 4.5, 'Contact'),
    (282, 'R282W', 1.8, 'Structural'),
]

hotspot_y = y_main + 3
for i, (pos, name, freq, mut_type) in enumerate(hotspots):
    color = '#FF006E' if mut_type == 'Structural' else '#00B4D8'
    
    # 3D sphere effect
    circle_outer = Circle((pos, hotspot_y), 6, facecolor=color, 
                          edgecolor='#1a1a1a', linewidth=1.5, zorder=6)
    ax.add_patch(circle_outer)
    
    # Highlight
    circle_highlight = Circle((pos - 1.5, hotspot_y + 1.5), 2, facecolor='white',
                              edgecolor='none', alpha=0.4, zorder=7)
    ax.add_patch(circle_highlight)
    
    # Label with connector
    label_y = hotspot_y + 18 + (i % 2) * 12
    ax.plot([pos, pos], [hotspot_y + 6, label_y - 3], color=color, lw=1.5, zorder=5)
    ax.text(pos, label_y, f'{name}\n({freq}%)', ha='center', va='bottom',
           fontsize=8, fontweight='bold', color=color)

# ================================================================
# LEGEND
# ================================================================
legend_x = 10
legend_y = -50

# Domain colors legend
ax.text(legend_x, legend_y + 8, 'Domains:', fontsize=10, fontweight='bold', color='#2C3E50')
legend_items = [
    ('TAD', '#E74C3C'),
    ('PRR', '#F5B041'),
    ('DBD', '#2ECC71'),
    ('TD', '#3498DB'),
    ('CTD', '#9B59B6'),
]
for i, (name, color) in enumerate(legend_items):
    x = legend_x + 55 + i * 45
    rect = FancyBboxPatch((x, legend_y), 15, 10, boxstyle="round,pad=0.02,rounding_size=2",
                          facecolor=color, edgecolor='#1a1a1a', linewidth=1)
    ax.add_patch(rect)
    ax.text(x + 18, legend_y + 5, name, fontsize=9, va='center', fontweight='bold')

# Mutation type legend
ax.text(legend_x + 300, legend_y + 8, 'Mutations:', fontsize=10, fontweight='bold', color='#2C3E50')
mut_legend = [('Structural', '#FF006E'), ('DNA Contact', '#00B4D8')]
for i, (name, color) in enumerate(mut_legend):
    x = legend_x + 375 + i * 85
    circle = Circle((x, legend_y + 5), 5, facecolor=color, edgecolor='#1a1a1a', linewidth=1)
    ax.add_patch(circle)
    ax.text(x + 10, legend_y + 5, name, fontsize=9, va='center', fontweight='bold')

# ================================================================
# SAVE FIGURE
# ================================================================
plt.tight_layout()
plt.savefig(output_dir / 'Figure1_p53_3D_Domains.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure1_p53_3D_Domains.pdf', bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure1_p53_3D_Domains.tiff', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

print("\n" + "="*70)
print("3D DOMAIN FIGURE CREATED SUCCESSFULLY")
print("="*70)
print(f"\nSaved:")
print(f"  - {output_dir}/Figure1_p53_3D_Domains.png")
print(f"  - {output_dir}/Figure1_p53_3D_Domains.pdf")
print(f"  - {output_dir}/Figure1_p53_3D_Domains.tiff")
print("\nFigure shows:")
print("  ✓ Full-length p53 (1-393 residues)")
print("  ✓ All domains with 3D effect")
print("  ✓ DBD highlighted (This Study)")
print("  ✓ L2, L3 loops and H2 helix marked")
print("  ✓ Cancer hotspot mutations")
print("  ✓ Residue scale")
print("  ✓ Functional annotations")
