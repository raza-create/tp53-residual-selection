#!/usr/bin/env python3
"""
Composite Introduction Figure for p53 Paper
Combines PyMOL 3D structure with domain schematic
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle, FancyArrowPatch
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from PIL import Image
import numpy as np
from pathlib import Path

output_dir = Path('results/figures/publication')
output_dir.mkdir(parents=True, exist_ok=True)

print("="*70)
print("CREATING COMPOSITE INTRODUCTION FIGURE")
print("="*70)

# Create figure with GridSpec for precise layout
fig = plt.figure(figsize=(14, 18))

# ================================================================
# PANEL A: Domain Schematic (Top)
# ================================================================
ax_a = fig.add_axes([0.05, 0.72, 0.9, 0.25])
ax_a.set_xlim(0, 400)
ax_a.set_ylim(0, 100)
ax_a.axis('off')

# Title
ax_a.text(200, 95, 'A', fontsize=24, fontweight='bold', ha='center', va='top')
ax_a.text(200, 88, 'p53 Domain Architecture (393 amino acids)', 
         fontsize=14, fontweight='bold', ha='center', color='#1a1a1a')

# Domain bar settings
y_bar = 50
bar_h = 22

# Domain data: (name, start, end, color, short_label)
domains = [
    ('TAD1', 1, 40, '#DC3545', 'TAD1'),
    ('TAD2', 41, 61, '#FD7E14', 'TAD2'),
    ('PRD', 62, 93, '#FFC107', 'PRD'),
    ('DBD', 94, 292, '#28A745', 'DNA-Binding Domain (DBD)'),
    ('Linker', 293, 323, '#6C757D', ''),
    ('OD', 324, 355, '#6F42C1', 'OD'),
    ('CTD', 356, 393, '#17A2B8', 'CTD'),
]

# Draw domains
for name, start, end, color, label in domains:
    # Main rectangle
    rect = FancyBboxPatch((start, y_bar - bar_h/2), end-start, bar_h,
                          boxstyle="round,pad=0.01,rounding_size=2",
                          facecolor=color, edgecolor='#343a40', linewidth=1.5)
    ax_a.add_patch(rect)
    
    # Label on bar
    if end - start > 30:
        ax_a.text((start+end)/2, y_bar, label, ha='center', va='center',
                 fontsize=10 if name != 'DBD' else 11, fontweight='bold', color='white')

# Residue labels below
residue_labels = [(1, '1'), (40, '40'), (61, '61'), (93, '93'), (94, '94'), 
                  (292, '292'), (323, '323'), (355, '355'), (393, '393')]
for pos, label in residue_labels:
    ax_a.text(pos, y_bar - bar_h/2 - 8, label, ha='center', fontsize=8, color='#495057')
    ax_a.plot([pos, pos], [y_bar - bar_h/2, y_bar - bar_h/2 - 4], 'k-', lw=1)

# N and C terminus
ax_a.text(1, y_bar + bar_h/2 + 8, 'N', fontsize=14, fontweight='bold', 
         ha='center', color='#28A745')
ax_a.text(393, y_bar + bar_h/2 + 8, 'C', fontsize=14, fontweight='bold',
         ha='center', color='#DC3545')

# Domain function labels (above)
func_labels = [
    (20, 'Transactivation', '#DC3545'),
    (77, 'Proline-rich', '#FFC107'),
    (193, 'DNA Binding & Zn²⁺', '#28A745'),
    (340, 'Tetramerization', '#6F42C1'),
    (375, 'Regulation', '#17A2B8'),
]
for x, text, color in func_labels:
    ax_a.annotate(text, xy=(x, y_bar + bar_h/2), xytext=(x, y_bar + bar_h/2 + 20),
                 ha='center', fontsize=9, color=color,
                 arrowprops=dict(arrowstyle='->', color=color, lw=1.2))

# Highlight DBD (This Study)
highlight = FancyBboxPatch((92, y_bar - bar_h/2 - 3), 202, bar_h + 6,
                           boxstyle="round,pad=0.01,rounding_size=3",
                           facecolor='none', edgecolor='#DC3545', 
                           linewidth=3, linestyle='--')
ax_a.add_patch(highlight)
ax_a.text(193, y_bar - bar_h/2 - 12, '★ This Study', ha='center', 
         fontsize=11, fontweight='bold', color='#DC3545')

# Legend
legend_x = 10
legend_y = 15
ax_a.text(legend_x, legend_y + 8, 'Domains:', fontsize=10, fontweight='bold')
legend_items = [('TAD', '#DC3545'), ('PRD', '#FFC107'), ('DBD', '#28A745'), 
                ('OD', '#6F42C1'), ('CTD', '#17A2B8')]
for i, (name, color) in enumerate(legend_items):
    x = legend_x + 50 + i * 55
    rect = Rectangle((x, legend_y), 15, 10, facecolor=color, edgecolor='#343a40', lw=1)
    ax_a.add_patch(rect)
    ax_a.text(x + 18, legend_y + 5, name, fontsize=9, va='center')

# ================================================================
# PANEL B: 3D Structure (Middle - PyMOL image)
# ================================================================
ax_b = fig.add_axes([0.05, 0.28, 0.55, 0.42])
ax_b.axis('off')
ax_b.text(0.02, 0.98, 'B', fontsize=24, fontweight='bold', transform=ax_b.transAxes, va='top')
ax_b.set_title('3D Structure of p53 DBD (PDB: 2XWR)', fontsize=12, fontweight='bold', loc='left', pad=10)

# Load PyMOL image
pymol_paths = [
    'results/figures/pymol/Fig1B_p53_hotspots.png',
    'results/figures/pymol/Fig1A_p53_cartoon.png',
    'results/figures/publication/Figure1_Complete_Annotated.png',
]

img_loaded = False
for img_path in pymol_paths:
    if Path(img_path).exists():
        try:
            img = Image.open(img_path)
            ax_b.imshow(img)
            img_loaded = True
            print(f"Loaded: {img_path}")
            break
        except:
            continue

if not img_loaded:
    ax_b.text(0.5, 0.5, '3D Structure\n(PyMOL)', ha='center', va='center',
             transform=ax_b.transAxes, fontsize=14,
             bbox=dict(boxstyle='round', facecolor='#e9ecef'))

# Structure annotation (right side of image)
ax_b.text(1.02, 0.85, 'Key Features:', fontsize=10, fontweight='bold',
         transform=ax_b.transAxes, va='top')
features = [
    ('• β-sandwich core', '#FD7E14'),
    ('• L2 loop (Zn²⁺)', '#DC3545'),
    ('• L3 loop (DNA)', '#DC3545'),
    ('• H2 helix', '#17A2B8'),
    ('• Zinc ion (Zn²⁺)', '#FFC107'),
]
for i, (text, color) in enumerate(features):
    ax_b.text(1.02, 0.78 - i*0.08, text, fontsize=9, color=color,
             transform=ax_b.transAxes, va='top')

# ================================================================
# PANEL C: Cancer Hotspots (Right side)
# ================================================================
ax_c = fig.add_axes([0.62, 0.50, 0.35, 0.20])
ax_c.axis('off')
ax_c.text(0.02, 0.98, 'C', fontsize=24, fontweight='bold', transform=ax_c.transAxes, va='top')
ax_c.set_title('Cancer Mutation Hotspots', fontsize=12, fontweight='bold', loc='left', pad=5)

# Hotspot data
hotspots = [
    ('R248Q/W', 6.4, 'Contact', '#3498DB'),
    ('R175H', 4.8, 'Structural', '#E74C3C'),
    ('R273H/C', 4.5, 'Contact', '#3498DB'),
    ('G245S', 2.7, 'Structural', '#E74C3C'),
    ('R249S', 2.5, 'Structural', '#E74C3C'),
    ('R282W', 1.8, 'Structural', '#E74C3C'),
]

# Create horizontal bar chart
y_pos = np.arange(len(hotspots))
freqs = [h[1] for h in hotspots]
colors = [h[3] for h in hotspots]
names = [h[0] for h in hotspots]

bars = ax_c.barh(y_pos, freqs, color=colors, edgecolor='#343a40', linewidth=1, height=0.7)
ax_c.set_yticks(y_pos)
ax_c.set_yticklabels(names, fontsize=10, fontweight='bold')
ax_c.set_xlabel('Mutation Frequency (%)', fontsize=10)
ax_c.set_xlim(0, 8)
ax_c.invert_yaxis()

# Add frequency labels
for bar, freq in zip(bars, freqs):
    ax_c.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
             f'{freq}%', va='center', fontsize=9)

# Legend
ax_c.plot([], [], 's', color='#E74C3C', markersize=10, label='Structural')
ax_c.plot([], [], 's', color='#3498DB', markersize=10, label='DNA Contact')
ax_c.legend(loc='lower right', fontsize=8, frameon=True)

ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)

# ================================================================
# PANEL D: Zinc Coordination (Right side bottom)
# ================================================================
ax_d = fig.add_axes([0.62, 0.28, 0.35, 0.18])
ax_d.set_xlim(0, 100)
ax_d.set_ylim(0, 100)
ax_d.axis('off')
ax_d.text(0.02, 0.98, 'D', fontsize=24, fontweight='bold', transform=ax_d.transAxes, va='top')
ax_d.text(5, 92, 'Zinc Coordination Site', fontsize=12, fontweight='bold')

# Draw zinc coordination diagram
zn_x, zn_y = 50, 50

# Zinc ion (center)
zn_circle = Circle((zn_x, zn_y), 8, facecolor='#FFC107', edgecolor='#B8860B', linewidth=2)
ax_d.add_patch(zn_circle)
ax_d.text(zn_x, zn_y, 'Zn²⁺', ha='center', va='center', fontsize=11, fontweight='bold')

# Coordinating residues
coord_residues = [
    (zn_x - 25, zn_y + 20, 'Cys176'),
    (zn_x + 25, zn_y + 20, 'His179'),
    (zn_x - 25, zn_y - 20, 'Cys238'),
    (zn_x + 25, zn_y - 20, 'Cys242'),
]

for cx, cy, name in coord_residues:
    # Residue box
    res_box = FancyBboxPatch((cx - 12, cy - 6), 24, 12,
                             boxstyle="round,pad=0.02,rounding_size=1",
                             facecolor='#28A745', edgecolor='#1e7e34', linewidth=1.5)
    ax_d.add_patch(res_box)
    ax_d.text(cx, cy, name, ha='center', va='center', fontsize=9, 
             fontweight='bold', color='white')
    
    # Coordination bond
    ax_d.annotate('', xy=(zn_x + (cx-zn_x)*0.35, zn_y + (cy-zn_y)*0.35),
                 xytext=(cx - (cx-zn_x)*0.15, cy - (cy-zn_y)*0.15),
                 arrowprops=dict(arrowstyle='-', color='#6C757D', lw=2, 
                                connectionstyle='arc3'))

ax_d.text(50, 10, 'Tetrahedral coordination essential\nfor DNA-binding loop stability',
         ha='center', fontsize=9, style='italic', color='#495057')

# ================================================================
# PANEL E: Study Focus Summary (Bottom)
# ================================================================
ax_e = fig.add_axes([0.05, 0.03, 0.9, 0.22])
ax_e.set_xlim(0, 100)
ax_e.set_ylim(0, 100)
ax_e.axis('off')
ax_e.text(0.02, 0.98, 'E', fontsize=24, fontweight='bold', transform=ax_e.transAxes, va='top')
ax_e.text(50, 92, 'Study Overview: Thermal Stability Analysis of p53 DBD', 
         fontsize=12, fontweight='bold', ha='center')

# Workflow boxes
workflow = [
    (8, 'MD Simulations\n6 temperatures\n248K - 450K\n100 ns each', '#3498DB'),
    (30, 'Feature\nExtraction\nRMSD, Rg, H-bonds\nNative contacts', '#28A745'),
    (52, 'Machine\nLearning\nRandom Forest\nR² = 0.994', '#6F42C1'),
    (74, 'SHAP\nAnalysis\nFeature importance\nExplainability', '#E74C3C'),
]

for x, text, color in workflow:
    box = FancyBboxPatch((x, 25), 18, 55,
                         boxstyle="round,pad=0.02,rounding_size=2",
                         facecolor=color, edgecolor='#343a40', linewidth=2, alpha=0.9)
    ax_e.add_patch(box)
    ax_e.text(x + 9, 52, text, ha='center', va='center', fontsize=9,
             fontweight='bold', color='white')

# Arrows between boxes
for i in range(3):
    x_start = 8 + i * 22 + 18
    x_end = 8 + (i+1) * 22
    ax_e.annotate('', xy=(x_end, 52), xytext=(x_start, 52),
                 arrowprops=dict(arrowstyle='->', color='#343a40', lw=2.5))

# Key finding box
key_box = FancyBboxPatch((5, 2), 90, 18,
                         boxstyle="round,pad=0.02,rounding_size=2",
                         facecolor='#FFF3CD', edgecolor='#FFC107', linewidth=3)
ax_e.add_patch(key_box)
ax_e.text(50, 11, 'Key Finding: Native contacts (79.6% SHAP importance) dominate p53 thermal response',
         ha='center', va='center', fontsize=11, fontweight='bold', color='#856404')

# ================================================================
# SAVE
# ================================================================
plt.savefig(output_dir / 'Figure1_Introduction_Composite.png', dpi=300, 
            bbox_inches='tight', facecolor='white')
plt.savefig(output_dir / 'Figure1_Introduction_Composite.pdf', bbox_inches='tight')
plt.savefig(output_dir / 'Figure1_Introduction_Composite.tiff', dpi=300, 
            bbox_inches='tight', facecolor='white')
plt.close()

print("\n" + "="*70)
print("COMPOSITE INTRODUCTION FIGURE CREATED")
print("="*70)
print(f"\nSaved:")
print(f"  ✓ {output_dir}/Figure1_Introduction_Composite.png")
print(f"  ✓ {output_dir}/Figure1_Introduction_Composite.pdf")
print(f"  ✓ {output_dir}/Figure1_Introduction_Composite.tiff")
print("\nPanels:")
print("  A: Domain architecture (1-393)")
print("  B: 3D PyMOL structure")
print("  C: Cancer hotspots bar chart")
print("  D: Zinc coordination diagram")
print("  E: Study workflow + key finding")
