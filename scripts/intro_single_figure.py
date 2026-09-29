#!/usr/bin/env python3
"""
Create single-panel introduction figure: p53 domain organization
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrow
import numpy as np
from pathlib import Path

output_dir = Path('results/figures/publication')
output_dir.mkdir(parents=True, exist_ok=True)

# Publication style
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 11,
    'axes.linewidth': 1.5,
    'figure.dpi': 300,
})

print("Creating Introduction Figure...")

# Single figure
fig, ax = plt.subplots(figsize=(12, 6))
ax.set_xlim(0, 420)
ax.set_ylim(0, 100)
ax.axis('off')

# ================================================================
# p53 Full-length protein bar
# ================================================================
y_main = 70
bar_height = 18

# Domain definitions (residue positions)
domains = [
    ('TAD1', 1, 40, '#FF6B6B', 'TAD1\n(1-40)'),
    ('TAD2', 41, 61, '#FF9999', 'TAD2'),
    ('PRR', 64, 92, '#FFD93D', 'PRR\n(64-92)'),
    ('DBD', 94, 292, '#6BCB77', 'DNA-Binding Domain\n(94-292)'),
    ('OD', 325, 356, '#4D96FF', 'OD\n(325-356)'),
    ('CTD', 363, 393, '#9B59B6', 'CTD'),
]

# Draw domains
for name, start, end, color, label in domains:
    width = end - start
    rect = FancyBboxPatch((start, y_main - bar_height/2), width, bar_height,
                          boxstyle="round,pad=0.01,rounding_size=2",
                          facecolor=color, edgecolor='black', linewidth=2)
    ax.add_patch(rect)
    
    # Label positioning
    if name == 'DBD':
        ax.text(start + width/2, y_main, label, ha='center', va='center',
               fontsize=12, fontweight='bold', color='white')
    elif name in ['TAD2', 'CTD']:
        ax.text(start + width/2, y_main + bar_height/2 + 8, label, ha='center', va='bottom',
               fontsize=9, fontweight='bold')
    else:
        ax.text(start + width/2, y_main - bar_height/2 - 3, label, ha='center', va='top',
               fontsize=9, fontweight='bold')

# Highlight DBD with red dashed box
highlight = Rectangle((90, y_main - bar_height/2 - 5), 206, bar_height + 10,
                      fill=False, edgecolor='red', linewidth=3, linestyle='--')
ax.add_patch(highlight)

# "This Study" label with arrow
ax.annotate('This Study', xy=(193, y_main - bar_height/2 - 5), xytext=(193, y_main - 35),
           fontsize=12, fontweight='bold', color='red', ha='center',
           arrowprops=dict(arrowstyle='->', color='red', lw=2))

# Residue number markers
markers = [(1, '1'), (94, '94'), (292, '292'), (393, '393')]
for pos, label in markers:
    ax.plot([pos, pos], [y_main + bar_height/2, y_main + bar_height/2 + 5], 'k-', lw=1.5)
    ax.text(pos, y_main + bar_height/2 + 8, label, ha='center', va='bottom', fontsize=10)

# Scale bar
ax.plot([1, 393], [y_main + bar_height/2 + 5, y_main + bar_height/2 + 5], 'k-', lw=1.5)

# ================================================================
# Cancer hotspot mutations on DBD
# ================================================================
hotspots = [
    (175, 'R175H', 4.8, '#E74C3C'),  # Structural
    (245, 'G245S', 2.7, '#E74C3C'),  # Structural
    (248, 'R248Q', 6.4, '#3498DB'),  # Contact
    (273, 'R273H', 4.5, '#3498DB'),  # Contact
    (282, 'R282W', 1.8, '#E74C3C'),  # Structural
]

# Draw hotspot markers below DBD
y_hotspot = y_main - bar_height/2 - 18
for pos, name, freq, color in hotspots:
    # Vertical line from DBD to hotspot marker
    ax.plot([pos, pos], [y_main - bar_height/2, y_hotspot + 3], 
           color=color, lw=1.5, linestyle=':')
    
    # Circle marker
    circle = Circle((pos, y_hotspot), 4, facecolor=color, edgecolor='black', linewidth=1.5)
    ax.add_patch(circle)
    
    # Label
    ax.text(pos, y_hotspot - 8, f'{name}\n({freq}%)', ha='center', va='top', 
           fontsize=8, fontweight='bold')

# Hotspot legend
ax.text(320, y_hotspot + 5, 'Mutation Type:', fontsize=10, fontweight='bold')
structural_circle = Circle((340, y_hotspot - 8), 4, facecolor='#E74C3C', edgecolor='black')
contact_circle = Circle((340, y_hotspot - 20), 4, facecolor='#3498DB', edgecolor='black')
ax.add_patch(structural_circle)
ax.add_patch(contact_circle)
ax.text(350, y_hotspot - 8, 'Structural', fontsize=9, va='center')
ax.text(350, y_hotspot - 20, 'DNA Contact', fontsize=9, va='center')

# ================================================================
# Title
# ================================================================
ax.text(200, 97, 'p53 Domain Organization and Cancer Mutation Hotspots', 
       ha='center', va='top', fontsize=14, fontweight='bold')

# ================================================================
# Domain legend at bottom
# ================================================================
legend_y = 8
legend_items = [
    ('TAD', '#FF6B6B', 'Transactivation Domain'),
    ('PRR', '#FFD93D', 'Proline-Rich Region'),
    ('DBD', '#6BCB77', 'DNA-Binding Domain'),
    ('OD', '#4D96FF', 'Oligomerization Domain'),
    ('CTD', '#9B59B6', 'C-Terminal Domain'),
]

x_start = 30
for i, (abbr, color, full_name) in enumerate(legend_items):
    x_pos = x_start + i * 80
    rect = Rectangle((x_pos, legend_y - 3), 12, 8, facecolor=color, edgecolor='black', linewidth=1)
    ax.add_patch(rect)
    ax.text(x_pos + 15, legend_y + 1, f'{abbr}: {full_name}', fontsize=8, va='center')

plt.tight_layout()
plt.savefig(output_dir / 'Figure1_Introduction.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure1_Introduction.pdf', bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure1_Introduction.tiff', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

print("="*60)
print("INTRODUCTION FIGURE CREATED")
print("="*60)
print(f"\nSaved:")
print(f"  - {output_dir}/Figure1_Introduction.png")
print(f"  - {output_dir}/Figure1_Introduction.pdf")
print(f"  - {output_dir}/Figure1_Introduction.tiff")
print("\nFigure shows:")
print("  - p53 full-length domain organization (393 aa)")
print("  - DNA-Binding Domain highlighted (this study)")
print("  - Cancer mutation hotspots with frequencies")
print("  - Structural vs Contact mutation classification")
