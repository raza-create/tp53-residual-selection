#!/usr/bin/env python3
"""
Complete Full-Length p53 Figure (1-393 residues)
For high-impact journal introduction
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrowPatch, Polygon
from matplotlib.collections import PatchCollection
import numpy as np
from pathlib import Path
from PIL import Image

output_dir = Path('results/figures/publication')
output_dir.mkdir(parents=True, exist_ok=True)

# Publication settings
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.linewidth': 1.5,
    'figure.dpi': 300,
})

print("="*70)
print("CREATING FULL-LENGTH p53 FIGURE (1-393)")
print("="*70)

# Create figure
fig = plt.figure(figsize=(16, 20))

# ================================================================
# PANEL A: Full-length p53 domain organization (1-393)
# ================================================================
ax_a = fig.add_axes([0.05, 0.75, 0.9, 0.22])
ax_a.set_xlim(0, 430)
ax_a.set_ylim(0, 100)
ax_a.axis('off')
ax_a.text(0.5, 0.98, 'A', transform=ax_a.transAxes, fontsize=20, fontweight='bold', va='top')
ax_a.set_title('Full-Length p53 Domain Organization (393 residues)', fontsize=14, fontweight='bold', pad=10)

# Scale: 1 residue = 1 unit
y_main = 50
bar_height = 25

# ================================================================
# COMPLETE DOMAIN DEFINITIONS (1-393)
# ================================================================
domains = [
    # (name, start, end, color, label, label_pos)
    ('TAD1', 1, 40, '#E74C3C', 'TAD1\n(1-40)', 'below'),
    ('TAD2', 41, 61, '#FF6B6B', 'TAD2\n(41-61)', 'above'),
    ('PRR', 62, 93, '#F39C12', 'PRR\n(62-93)', 'below'),
    ('DBD', 94, 292, '#27AE60', 'DNA-Binding Domain\n(94-292)', 'center'),
    ('Linker', 293, 323, '#BDC3C7', 'Linker\n(293-323)', 'above'),
    ('OD', 324, 355, '#3498DB', 'OD\n(324-355)', 'below'),
    ('CTD', 356, 393, '#9B59B6', 'CTD\n(356-393)', 'above'),
]

# Draw domains
for name, start, end, color, label, label_pos in domains:
    width = end - start + 1
    
    # Main domain box
    rect = FancyBboxPatch((start, y_main - bar_height/2), width, bar_height,
                          boxstyle="round,pad=0.01,rounding_size=3",
                          facecolor=color, edgecolor='black', linewidth=2)
    ax_a.add_patch(rect)
    
    # Labels
    if label_pos == 'center':
        ax_a.text(start + width/2, y_main, label, ha='center', va='center',
                 fontsize=11, fontweight='bold', color='white')
    elif label_pos == 'above':
        ax_a.text(start + width/2, y_main + bar_height/2 + 5, label, ha='center', va='bottom',
                 fontsize=9, fontweight='bold')
    else:  # below
        ax_a.text(start + width/2, y_main - bar_height/2 - 5, label, ha='center', va='top',
                 fontsize=9, fontweight='bold')

# Highlight DBD (This Study)
highlight = Rectangle((90, y_main - bar_height/2 - 8), 206, bar_height + 16,
                      fill=False, edgecolor='red', linewidth=3, linestyle='--')
ax_a.add_patch(highlight)
ax_a.text(193, y_main + bar_height/2 + 12, '★ This Study', ha='center', va='bottom',
         fontsize=12, fontweight='bold', color='red')

# Scale bar with tick marks
ax_a.plot([1, 393], [y_main + bar_height/2 + 25, y_main + bar_height/2 + 25], 'k-', lw=2)
for pos in [1, 50, 100, 150, 200, 250, 300, 350, 393]:
    ax_a.plot([pos, pos], [y_main + bar_height/2 + 23, y_main + bar_height/2 + 27], 'k-', lw=1.5)
    ax_a.text(pos, y_main + bar_height/2 + 29, str(pos), ha='center', va='bottom', fontsize=8)

# Domain function annotations
functions = [
    (20, 'Transcriptional\nactivation'),
    (77, 'SH3 domain\nbinding'),
    (193, 'DNA binding\n& Zn coordination'),
    (340, 'Tetramerization'),
    (375, 'Regulation'),
]
for x, text in functions:
    ax_a.text(x, y_main - bar_height/2 - 25, text, ha='center', va='top', fontsize=8, style='italic')

# ================================================================
# PANEL B: DBD Secondary Structure Detail (94-292)
# ================================================================
ax_b = fig.add_axes([0.05, 0.42, 0.9, 0.30])
ax_b.set_xlim(90, 296)
ax_b.set_ylim(0, 100)
ax_b.axis('off')
ax_b.text(0.0, 0.98, 'B', transform=ax_b.transAxes, fontsize=20, fontweight='bold', va='top')
ax_b.set_title('DNA-Binding Domain Secondary Structure (94-292)', fontsize=14, fontweight='bold', pad=10)

y_ss = 60
ss_height = 20

# Complete secondary structure of DBD
# Based on PDB 2XWR annotations
secondary_structure = [
    # (type, start, end, color, label)
    # N-terminal
    ('coil', 94, 102, '#BDC3C7', ''),
    
    # Beta sandwich - Sheet 1
    ('strand', 103, 109, '#F39C12', 'S1'),
    ('coil', 110, 111, '#BDC3C7', ''),
    ('strand', 112, 117, '#F39C12', 'S2'),
    ('coil', 118, 123, '#BDC3C7', ''),
    ('strand', 124, 127, '#F39C12', 'S3'),
    ('strand', 128, 132, '#F39C12', 'S4'),
    
    # L1 loop
    ('coil', 133, 141, '#A569BD', 'L1'),
    
    # Continue beta
    ('strand', 142, 146, '#F39C12', 'S5'),
    ('coil', 147, 150, '#BDC3C7', ''),
    ('strand', 151, 155, '#F39C12', 'S6'),
    ('coil', 156, 162, '#BDC3C7', ''),
    ('strand', 163, 166, '#F39C12', 'S7'),
    
    # L2 loop with H1 helix (zinc binding region)
    ('coil', 167, 176, '#E74C3C', 'L2a'),
    ('helix', 177, 181, '#1ABC9C', 'H1'),
    ('coil', 182, 195, '#E74C3C', 'L2b'),
    
    # Continue beta
    ('strand', 196, 198, '#F39C12', 'S8'),
    ('strand', 199, 205, '#F39C12', 'S9'),
    ('coil', 206, 212, '#BDC3C7', ''),
    ('strand', 213, 220, '#F39C12', 'S10'),
    ('coil', 221, 228, '#BDC3C7', ''),
    ('strand', 229, 232, '#F39C12', 'S11'),
    ('coil', 233, 236, '#BDC3C7', ''),
    
    # L3 loop (DNA binding)
    ('coil', 237, 250, '#C0392B', 'L3'),
    
    # Continue beta
    ('strand', 251, 254, '#F39C12', 'S12'),
    ('coil', 255, 260, '#BDC3C7', ''),
    ('strand', 261, 265, '#F39C12', 'S13'),
    ('coil', 266, 277, '#BDC3C7', ''),
    
    # H2 helix (DNA binding)
    ('helix', 278, 286, '#16A085', 'H2'),
    
    # C-terminal
    ('coil', 287, 292, '#BDC3C7', ''),
]

# Draw secondary structure elements
for ss_type, start, end, color, label in secondary_structure:
    width = end - start + 1
    
    if ss_type == 'helix':
        # Draw helix as cylinder/rectangle with rounded ends
        rect = FancyBboxPatch((start, y_ss - ss_height/2), width, ss_height,
                              boxstyle="round,pad=0.01,rounding_size=2",
                              facecolor=color, edgecolor='black', linewidth=1.5)
        ax_b.add_patch(rect)
        # Add helix pattern
        for i in range(start, end, 2):
            ax_b.plot([i, i+1], [y_ss - ss_height/4, y_ss + ss_height/4], 
                     color='white', lw=1, alpha=0.5)
    elif ss_type == 'strand':
        # Draw strand as arrow
        arrow_width = ss_height * 0.8
        points = [
            (start, y_ss - arrow_width/2),
            (end - 2, y_ss - arrow_width/2),
            (end - 2, y_ss - ss_height/2),
            (end + 1, y_ss),
            (end - 2, y_ss + ss_height/2),
            (end - 2, y_ss + arrow_width/2),
            (start, y_ss + arrow_width/2),
        ]
        arrow = Polygon(points, facecolor=color, edgecolor='black', linewidth=1)
        ax_b.add_patch(arrow)
    else:  # coil
        rect = Rectangle((start, y_ss - ss_height/4), width, ss_height/2,
                         facecolor=color, edgecolor='none')
        ax_b.add_patch(rect)
    
    # Add labels for named elements
    if label and label not in ['L2a', 'L2b']:
        y_offset = ss_height/2 + 8 if 'S' in label else ss_height/2 + 5
        ax_b.text(start + width/2, y_ss + y_offset, label, ha='center', va='bottom',
                 fontsize=8, fontweight='bold')

# Special labels for L2 (combined)
ax_b.annotate('L2 Loop\n(167-195)\nZn binding', xy=(180, y_ss - ss_height/2 - 2),
             xytext=(180, y_ss - ss_height/2 - 25),
             ha='center', fontsize=9, fontweight='bold', color='#E74C3C',
             arrowprops=dict(arrowstyle='->', color='#E74C3C', lw=1.5))

ax_b.annotate('L3 Loop\n(237-250)\nDNA contact', xy=(243, y_ss - ss_height/2 - 2),
             xytext=(243, y_ss - ss_height/2 - 25),
             ha='center', fontsize=9, fontweight='bold', color='#C0392B',
             arrowprops=dict(arrowstyle='->', color='#C0392B', lw=1.5))

# Zinc coordinating residues markers
zn_residues = [176, 179, 238, 242]
for res in zn_residues:
    ax_b.plot(res, y_ss + ss_height/2 + 3, 'o', color='#F4D03F', markersize=8, 
             markeredgecolor='black', markeredgewidth=1)
ax_b.text(177, y_ss + ss_height/2 + 12, 'Zn²⁺ ligands\n(C176,H179,C238,C242)', 
         ha='center', fontsize=8, fontweight='bold', color='#F4D03F')

# Cancer hotspot markers
hotspots = [
    (175, 'R175H', '#FF006E', 4.8),
    (245, 'G245S', '#FF006E', 2.7),
    (248, 'R248Q', '#00D4FF', 6.4),
    (249, 'R249S', '#FF006E', 2.5),
    (273, 'R273H', '#00D4FF', 4.5),
    (282, 'R282W', '#FF006E', 1.8),
]

for res, name, color, freq in hotspots:
    ax_b.plot(res, y_ss - ss_height/2 - 3, 's', color=color, markersize=10,
             markeredgecolor='black', markeredgewidth=1.5)
    ax_b.text(res, y_ss - ss_height/2 - 8, f'{name}\n({freq}%)', ha='center', va='top',
             fontsize=7, fontweight='bold')

# Scale bar
ax_b.plot([94, 292], [y_ss + ss_height/2 + 22, y_ss + ss_height/2 + 22], 'k-', lw=1.5)
for pos in [94, 120, 150, 180, 210, 240, 270, 292]:
    ax_b.plot([pos, pos], [y_ss + ss_height/2 + 20, y_ss + ss_height/2 + 24], 'k-', lw=1)
    ax_b.text(pos, y_ss + ss_height/2 + 26, str(pos), ha='center', va='bottom', fontsize=7)

# Legend for Panel B
legend_y = 8
legend_items = [
    ('β-strand', '#F39C12'),
    ('α-helix', '#1ABC9C'),
    ('L1 loop', '#A569BD'),
    ('L2 loop (Zn)', '#E74C3C'),
    ('L3 loop (DNA)', '#C0392B'),
    ('Coil', '#BDC3C7'),
]
for i, (name, color) in enumerate(legend_items):
    x = 100 + i * 33
    rect = Rectangle((x, legend_y), 8, 6, facecolor=color, edgecolor='black', linewidth=1)
    ax_b.add_patch(rect)
    ax_b.text(x + 10, legend_y + 3, name, fontsize=8, va='center')

# Hotspot legend
ax_b.plot(280, legend_y + 3, 's', color='#FF006E', markersize=8, markeredgecolor='black')
ax_b.text(283, legend_y + 3, 'Structural mut.', fontsize=8, va='center')
ax_b.plot(280, legend_y - 5, 's', color='#00D4FF', markersize=8, markeredgecolor='black')
ax_b.text(283, legend_y - 5, 'Contact mut.', fontsize=8, va='center')

# ================================================================
# PANEL C: 3D Structure (Load PyMOL image)
# ================================================================
ax_c = fig.add_axes([0.05, 0.02, 0.55, 0.38])
ax_c.text(0.0, 0.98, 'C', transform=ax_c.transAxes, fontsize=20, fontweight='bold', va='top')
ax_c.set_title('3D Structure of p53 DBD (PDB: 2XWR)', fontsize=14, fontweight='bold', pad=10)

# Try to load PyMOL figure
pymol_img_path = Path('results/figures/publication/Figure1_Complete_Annotated.png')
if not pymol_img_path.exists():
    pymol_img_path = Path('results/figures/pymol/Fig1A_p53_cartoon.png')

if pymol_img_path.exists():
    try:
        img = Image.open(pymol_img_path)
        ax_c.imshow(img)
    except Exception as e:
        ax_c.text(0.5, 0.5, f'3D Structure\n(See PyMOL figure)', ha='center', va='center',
                 transform=ax_c.transAxes, fontsize=12)
else:
    ax_c.text(0.5, 0.5, '3D Structure\n(PyMOL rendering)', ha='center', va='center',
             transform=ax_c.transAxes, fontsize=12,
             bbox=dict(boxstyle='round', facecolor='lightgray'))
ax_c.axis('off')

# ================================================================
# PANEL D: Key Information Table
# ================================================================
ax_d = fig.add_axes([0.62, 0.02, 0.35, 0.38])
ax_d.axis('off')
ax_d.text(0.0, 0.98, 'D', transform=ax_d.transAxes, fontsize=20, fontweight='bold', va='top')
ax_d.set_title('p53 Domain Summary', fontsize=14, fontweight='bold', pad=10)

# Create table data
table_data = [
    ['Domain', 'Residues', 'Function'],
    ['TAD1', '1-40', 'Transcriptional activation'],
    ['TAD2', '41-61', 'Secondary transactivation'],
    ['PRR', '62-93', 'Proline-rich, SH3 binding'],
    ['DBD', '94-292', 'DNA binding, Zn coordination'],
    ['  - L2 loop', '167-195', 'Zinc binding, thermal sensitive'],
    ['  - L3 loop', '237-250', 'DNA contact, most sensitive'],
    ['  - H2 helix', '278-286', 'DNA recognition helix'],
    ['Linker', '293-323', 'Flexible connection'],
    ['OD', '324-355', 'Tetramerization'],
    ['CTD', '356-393', 'Regulatory, post-translational mod.'],
]

# Draw table
cell_height = 0.075
cell_widths = [0.25, 0.2, 0.55]
y_start = 0.88

for i, row in enumerate(table_data):
    y = y_start - i * cell_height
    x = 0.02
    
    for j, (cell, width) in enumerate(zip(row, cell_widths)):
        # Header styling
        if i == 0:
            facecolor = '#2C3E50'
            textcolor = 'white'
            fontweight = 'bold'
        elif 'L2' in cell or 'L3' in cell or 'H2' in cell:
            facecolor = '#E8F6F3'
            textcolor = 'black'
            fontweight = 'normal'
        elif 'DBD' in cell:
            facecolor = '#D5F5E3'
            textcolor = 'black'
            fontweight = 'bold'
        else:
            facecolor = 'white'
            textcolor = 'black'
            fontweight = 'normal'
        
        rect = FancyBboxPatch((x, y - cell_height + 0.01), width - 0.01, cell_height - 0.005,
                              boxstyle="round,pad=0.01,rounding_size=0.01",
                              facecolor=facecolor, edgecolor='gray', linewidth=0.5,
                              transform=ax_d.transAxes)
        ax_d.add_patch(rect)
        ax_d.text(x + width/2, y - cell_height/2 + 0.005, cell, 
                 ha='center', va='center', fontsize=8, fontweight=fontweight,
                 color=textcolor, transform=ax_d.transAxes)
        x += width

# Cancer hotspot summary
ax_d.text(0.02, 0.08, 'Cancer Hotspot Mutations:', fontsize=10, fontweight='bold',
         transform=ax_d.transAxes)
hotspot_text = """• R175H (4.8%) - L2 loop, structural
- G245S (2.7%) - L3 loop, structural  
- R248Q/W (6.4%) - L3 loop, DNA contact
- R249S (2.5%) - L3 loop, structural
- R273H/C (4.5%) - pre-H2, DNA contact
- R282W (1.8%) - H2 helix, structural"""
ax_d.text(0.02, 0.02, hotspot_text, fontsize=8, transform=ax_d.transAxes,
         verticalalignment='bottom', family='monospace')

# Save figure
plt.savefig(output_dir / 'Figure1_FullLength_p53.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure1_FullLength_p53.pdf', bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure1_FullLength_p53.tiff', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

print("\n" + "="*70)
print("FIGURE CREATED SUCCESSFULLY")
print("="*70)
print(f"\nSaved:")
print(f"  - {output_dir}/Figure1_FullLength_p53.png")
print(f"  - {output_dir}/Figure1_FullLength_p53.pdf")
print(f"  - {output_dir}/Figure1_FullLength_p53.tiff")
print("\nFigure contains:")
print("  Panel A: Full-length p53 domain organization (1-393)")
print("  Panel B: DBD secondary structure detail (94-292)")
print("  Panel C: 3D structure (PyMOL)")
print("  Panel D: Domain summary table + cancer hotspots")
