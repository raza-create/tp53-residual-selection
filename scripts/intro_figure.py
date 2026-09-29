#!/usr/bin/env python3
"""
Create publication-quality introduction figure for p53 thermal analysis paper
Figure 1: Overview of p53 structure, cancer mutations, and study design
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch, Circle
from matplotlib.collections import PatchCollection
import numpy as np
from pathlib import Path
from PIL import Image

# Output directory
output_dir = Path('results/figures/publication')
output_dir.mkdir(parents=True, exist_ok=True)

# Set publication style
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.2,
})

print("="*70)
print("CREATING INTRODUCTION FIGURE")
print("="*70)

# Create figure with subplots
fig = plt.figure(figsize=(14, 16))

# ================================================================
# Panel A: p53 Domain Organization
# ================================================================
ax_a = fig.add_axes([0.05, 0.78, 0.9, 0.18])
ax_a.set_xlim(0, 100)
ax_a.set_ylim(0, 10)
ax_a.axis('off')
ax_a.set_title('A) p53 Domain Organization', fontsize=14, fontweight='bold', loc='left', x=-0.02)

# Full-length p53 bar (393 aa)
# Scale: 100 units = 393 aa
def aa_to_scale(aa):
    return aa * 100 / 393

# Domain colors
colors = {
    'TAD1': '#FF6B6B',
    'TAD2': '#FF8E8E', 
    'PRR': '#FFE66D',
    'DBD': '#4ECDC4',
    'TD': '#45B7D1',
    'CTD': '#96CEB4'
}

# Draw domains
domains = [
    ('TAD1', 1, 40, colors['TAD1'], 'Transactivation\nDomain 1'),
    ('TAD2', 40, 61, colors['TAD2'], 'TAD2'),
    ('PRR', 64, 92, colors['PRR'], 'Proline-Rich\nRegion'),
    ('DBD', 94, 292, colors['DBD'], 'DNA-Binding Domain\n(94-292)'),
    ('TD', 324, 355, colors['TD'], 'Tetramerization\nDomain'),
    ('CTD', 363, 393, colors['CTD'], 'C-Terminal\nDomain'),
]

y_center = 5
bar_height = 2.5

for name, start, end, color, label in domains:
    x_start = aa_to_scale(start)
    width = aa_to_scale(end - start)
    rect = FancyBboxPatch((x_start, y_center - bar_height/2), width, bar_height,
                          boxstyle="round,pad=0.02,rounding_size=0.3",
                          facecolor=color, edgecolor='black', linewidth=1.5)
    ax_a.add_patch(rect)
    
    # Add label
    if name == 'DBD':
        ax_a.text(x_start + width/2, y_center, label, ha='center', va='center',
                 fontsize=10, fontweight='bold', color='white')
        # Highlight DBD with box
        highlight = Rectangle((x_start-1, y_center - bar_height/2 - 1), width+2, bar_height + 2,
                              fill=False, edgecolor='red', linewidth=2, linestyle='--')
        ax_a.add_patch(highlight)
        ax_a.text(x_start + width/2, y_center - bar_height/2 - 1.8, 'This Study',
                 ha='center', va='top', fontsize=10, fontweight='bold', color='red')
    elif name in ['TAD1', 'PRR']:
        ax_a.text(x_start + width/2, y_center + bar_height/2 + 0.5, label, ha='center', va='bottom',
                 fontsize=8, fontweight='bold')
    else:
        ax_a.text(x_start + width/2, y_center - bar_height/2 - 0.3, label, ha='center', va='top',
                 fontsize=8, fontweight='bold')

# Add residue numbers
ax_a.text(aa_to_scale(1), y_center + bar_height/2 + 2.2, '1', ha='center', fontsize=8)
ax_a.text(aa_to_scale(393), y_center + bar_height/2 + 2.2, '393', ha='center', fontsize=8)
ax_a.text(aa_to_scale(94), y_center + bar_height/2 + 2.2, '94', ha='center', fontsize=8)
ax_a.text(aa_to_scale(292), y_center + bar_height/2 + 2.2, '292', ha='center', fontsize=8)

# Add scale bar
ax_a.plot([0, aa_to_scale(393)], [y_center + bar_height/2 + 1.8, y_center + bar_height/2 + 1.8], 
         'k-', linewidth=1)

# ================================================================
# Panel B: 3D Structure with Hotspots (load PyMOL image)
# ================================================================
ax_b = fig.add_axes([0.05, 0.42, 0.42, 0.32])
ax_b.set_title('B) p53 DBD Structure with Cancer Hotspots', fontsize=14, fontweight='bold', loc='left')

# Try to load the PyMOL figure
try:
    img_path = Path('results/figures/pymol/Fig1B_p53_hotspots.png')
    if img_path.exists():
        img = Image.open(img_path)
        ax_b.imshow(img)
        ax_b.axis('off')
    else:
        # Create placeholder
        ax_b.text(0.5, 0.5, 'p53 DBD Structure\nwith Hotspots\n\n(See Fig1B_p53_hotspots.png)', 
                 ha='center', va='center', fontsize=12, transform=ax_b.transAxes,
                 bbox=dict(boxstyle='round', facecolor='lightgray'))
        ax_b.set_xlim(0, 1)
        ax_b.set_ylim(0, 1)
        ax_b.axis('off')
except Exception as e:
    ax_b.text(0.5, 0.5, f'Structure image\n{e}', ha='center', va='center', transform=ax_b.transAxes)
    ax_b.axis('off')

# ================================================================
# Panel C: Cancer Mutation Frequency
# ================================================================
ax_c = fig.add_axes([0.55, 0.42, 0.4, 0.32])
ax_c.set_title('C) Cancer Mutation Hotspots', fontsize=14, fontweight='bold', loc='left')

# Hotspot data
hotspots = {
    'R175H': {'freq': 4.8, 'type': 'Structural', 'color': '#E74C3C'},
    'G245S': {'freq': 2.7, 'type': 'Structural', 'color': '#E74C3C'},
    'R248Q/W': {'freq': 6.4, 'type': 'Contact', 'color': '#3498DB'},
    'R249S': {'freq': 2.5, 'type': 'Structural', 'color': '#E74C3C'},
    'R273H/C': {'freq': 4.5, 'type': 'Contact', 'color': '#3498DB'},
    'R282W': {'freq': 1.8, 'type': 'Structural', 'color': '#E74C3C'},
}

mutations = list(hotspots.keys())
frequencies = [hotspots[m]['freq'] for m in mutations]
colors_bar = [hotspots[m]['color'] for m in mutations]

bars = ax_c.barh(mutations, frequencies, color=colors_bar, edgecolor='black', linewidth=1.2)

ax_c.set_xlabel('Mutation Frequency (%)', fontsize=11)
ax_c.set_xlim(0, 8)
ax_c.invert_yaxis()

# Add frequency labels
for bar, freq in zip(bars, frequencies):
    ax_c.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
             f'{freq}%', va='center', fontsize=10)

# Add legend
structural_patch = mpatches.Patch(color='#E74C3C', label='Structural mutation')
contact_patch = mpatches.Patch(color='#3498DB', label='Contact mutation')
ax_c.legend(handles=[structural_patch, contact_patch], loc='lower right', frameon=True)

ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)

# ================================================================
# Panel D: Study Workflow
# ================================================================
ax_d = fig.add_axes([0.05, 0.02, 0.9, 0.35])
ax_d.set_xlim(0, 100)
ax_d.set_ylim(0, 40)
ax_d.axis('off')
ax_d.set_title('D) Study Workflow', fontsize=14, fontweight='bold', loc='left', x=-0.02)

# Workflow boxes
workflow_steps = [
    {'x': 5, 'y': 25, 'w': 18, 'h': 12, 'color': '#3498DB', 
     'title': 'MD Simulations', 'text': '• 6 temperatures\n• 100 ns each\n• 60,006 frames'},
    {'x': 28, 'y': 25, 'w': 18, 'h': 12, 'color': '#27AE60',
     'title': 'Feature Extraction', 'text': '• RMSD, Rg\n• Native contacts\n• H-bonds'},
    {'x': 51, 'y': 25, 'w': 18, 'h': 12, 'color': '#9B59B6',
     'title': 'Machine Learning', 'text': '• Random Forest\n• LOTO validation\n• R² = 0.994'},
    {'x': 74, 'y': 25, 'w': 18, 'h': 12, 'color': '#E74C3C',
     'title': 'SHAP Analysis', 'text': '• Feature importance\n• Native contacts\n  = 79.6%'},
]

for step in workflow_steps:
    # Main box
    rect = FancyBboxPatch((step['x'], step['y']), step['w'], step['h'],
                          boxstyle="round,pad=0.02,rounding_size=0.5",
                          facecolor=step['color'], edgecolor='black', linewidth=2,
                          alpha=0.9)
    ax_d.add_patch(rect)
    
    # Title
    ax_d.text(step['x'] + step['w']/2, step['y'] + step['h'] - 1.5, step['title'],
             ha='center', va='top', fontsize=11, fontweight='bold', color='white')
    
    # Text
    ax_d.text(step['x'] + step['w']/2, step['y'] + step['h']/2 - 1, step['text'],
             ha='center', va='center', fontsize=9, color='white')

# Arrows between boxes
arrow_style = dict(arrowstyle='->', color='black', lw=2)
for i in range(len(workflow_steps) - 1):
    x1 = workflow_steps[i]['x'] + workflow_steps[i]['w']
    x2 = workflow_steps[i+1]['x']
    y = workflow_steps[i]['y'] + workflow_steps[i]['h']/2
    ax_d.annotate('', xy=(x2, y), xytext=(x1, y),
                 arrowprops=dict(arrowstyle='->', color='black', lw=2.5))

# Temperature scale below first box
temps = ['248K', '273K', '310K', '350K', '400K', '450K']
temp_colors = ['#0000FF', '#00BFFF', '#00FF00', '#FFFF00', '#FF8C00', '#FF0000']

for i, (temp, color) in enumerate(zip(temps, temp_colors)):
    x = 5 + i * 3
    circle = Circle((x + 1.5, 20), 1.2, facecolor=color, edgecolor='black', linewidth=0.5)
    ax_d.add_patch(circle)
    ax_d.text(x + 1.5, 17.5, temp, ha='center', va='top', fontsize=7, rotation=45)

ax_d.text(14, 15, 'Cold → Hot', ha='center', va='top', fontsize=9, style='italic')

# Key finding box at bottom
finding_box = FancyBboxPatch((20, 2), 60, 10,
                             boxstyle="round,pad=0.02,rounding_size=0.8",
                             facecolor='#FFF9C4', edgecolor='#F57F17', linewidth=3)
ax_d.add_patch(finding_box)

ax_d.text(50, 7, 'KEY FINDING', ha='center', va='center', fontsize=12, fontweight='bold',
         color='#E65100')
ax_d.text(50, 4, 'Native contacts (79.6%) dominate p53 thermal response\n'
         'r = −0.887 | Cohen\'s d = 8.16 | 47% contact loss at 450K',
         ha='center', va='center', fontsize=10)

# Arrow from SHAP to finding
ax_d.annotate('', xy=(50, 12), xytext=(83, 25),
             arrowprops=dict(arrowstyle='->', color='#E74C3C', lw=2.5,
                           connectionstyle='arc3,rad=-0.2'))

# Save figure
plt.savefig(output_dir / 'Figure1_Introduction.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure1_Introduction.pdf', bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

print(f"Saved: {output_dir}/Figure1_Introduction.png")
print(f"Saved: {output_dir}/Figure1_Introduction.pdf")

# ================================================================
# Create Figure 2: Results Summary Figure
# ================================================================
print("\nCreating Results Summary Figure...")

fig2, axes = plt.subplots(2, 3, figsize=(16, 10))

# Panel A: Native contacts vs Temperature
ax_a2 = axes[0, 0]
temps = [248, 273, 310, 350, 400, 450]
nc_means = [0.962, 0.962, 0.846, 0.642, 0.615, 0.509]
nc_stds = [0.012, 0.014, 0.096, 0.099, 0.106, 0.078]

ax_a2.errorbar(temps, nc_means, yerr=nc_stds, fmt='o-', markersize=10, 
              color='#2ECC71', linewidth=2, capsize=5, capthick=2,
              markeredgecolor='black', markeredgewidth=1.5)
ax_a2.fill_between(temps, np.array(nc_means)-np.array(nc_stds), 
                   np.array(nc_means)+np.array(nc_stds), alpha=0.3, color='#2ECC71')
ax_a2.set_xlabel('Temperature (K)', fontsize=11)
ax_a2.set_ylabel('Native Contact Fraction', fontsize=11)
ax_a2.set_title('A) Native Contact Loss', fontsize=12, fontweight='bold')
ax_a2.set_ylim(0.4, 1.05)

# Add correlation text
ax_a2.text(0.95, 0.95, f'r = −0.887\np < 0.001', transform=ax_a2.transAxes,
          ha='right', va='top', fontsize=10, 
          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

ax_a2.spines['top'].set_visible(False)
ax_a2.spines['right'].set_visible(False)

# Panel B: Feature Correlations
ax_b2 = axes[0, 1]
features = ['Native\ncontacts', 'RMSD', 'H-bonds', 'Asphericity', 'Rg', 'End-to-end']
correlations = [-0.887, 0.608, -0.556, 0.509, 0.490, -0.155]
colors_corr = ['#E74C3C' if c < 0 else '#3498DB' for c in correlations]

bars_b = ax_b2.barh(features, correlations, color=colors_corr, edgecolor='black', linewidth=1.2)
ax_b2.axvline(x=0, color='black', linewidth=1)
ax_b2.set_xlabel('Pearson Correlation (r)', fontsize=11)
ax_b2.set_title('B) Feature-Temperature Correlations', fontsize=12, fontweight='bold')
ax_b2.set_xlim(-1, 0.8)

for bar, corr in zip(bars_b, correlations):
    x_pos = corr + 0.03 if corr > 0 else corr - 0.03
    ha = 'left' if corr > 0 else 'right'
    ax_b2.text(x_pos, bar.get_y() + bar.get_height()/2, f'{corr:.3f}', 
              ha=ha, va='center', fontsize=9, fontweight='bold')

ax_b2.spines['top'].set_visible(False)
ax_b2.spines['right'].set_visible(False)

# Panel C: SHAP Importance
ax_c2 = axes[0, 2]
shap_features = ['Native contacts', 'End-to-end', 'Rg', 'RMSD', 'Asphericity', 'H-bonds']
shap_values = [79.6, 7.2, 6.1, 5.3, 1.1, 0.8]
colors_shap = plt.cm.Greens(np.linspace(0.3, 0.9, len(shap_features)))[::-1]

bars_c = ax_c2.barh(shap_features, shap_values, color=colors_shap, edgecolor='black', linewidth=1.2)
ax_c2.set_xlabel('SHAP Contribution (%)', fontsize=11)
ax_c2.set_title('C) SHAP Feature Importance', fontsize=12, fontweight='bold')

for bar, val in zip(bars_c, shap_values):
    ax_c2.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2, f'{val}%', 
              ha='left', va='center', fontsize=10, fontweight='bold')

ax_c2.set_xlim(0, 95)
ax_c2.spines['top'].set_visible(False)
ax_c2.spines['right'].set_visible(False)

# Panel D: LOTO Validation
ax_d2 = axes[1, 0]
loto_temps = [248, 273, 310, 350, 400, 450]
loto_predicted = [289.5, 326.6, 332.5, 349.6, 398.7, 398.9]
loto_mae = [41.5, 57.7, 38.8, 36.6, 40.0, 51.1]

x = np.arange(len(loto_temps))
width = 0.35

bars1 = ax_d2.bar(x - width/2, loto_temps, width, label='True', color='#3498DB', edgecolor='black')
bars2 = ax_d2.bar(x + width/2, loto_predicted, width, label='Predicted', color='#E74C3C', edgecolor='black')

ax_d2.set_xlabel('Test Temperature', fontsize=11)
ax_d2.set_ylabel('Temperature (K)', fontsize=11)
ax_d2.set_title('D) LOTO Validation', fontsize=12, fontweight='bold')
ax_d2.set_xticks(x)
ax_d2.set_xticklabels([f'{t}K' for t in loto_temps])
ax_d2.legend(loc='upper left')

ax_d2.text(0.95, 0.95, f'MAE = 44.3K', transform=ax_d2.transAxes,
          ha='right', va='top', fontsize=10, fontweight='bold',
          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

ax_d2.spines['top'].set_visible(False)
ax_d2.spines['right'].set_visible(False)

# Panel E: Effect Size
ax_e2 = axes[1, 1]
effect_features = ['Native\ncontacts', 'H-bonds', 'RMSD', 'Asphericity', 'Rg', 'End-to-end']
cohens_d = [8.16, 3.2, 2.8, 2.1, 1.9, 0.8]

colors_effect = ['#E74C3C' if d > 2 else '#F39C12' if d > 0.8 else '#3498DB' for d in cohens_d]
bars_e = ax_e2.barh(effect_features, cohens_d, color=colors_effect, edgecolor='black', linewidth=1.2)

ax_e2.axvline(x=0.8, color='gray', linestyle='--', linewidth=1.5, label='Large effect (0.8)')
ax_e2.set_xlabel("Cohen's d (248K vs 450K)", fontsize=11)
ax_e2.set_title('E) Effect Sizes', fontsize=12, fontweight='bold')
ax_e2.legend(loc='lower right')

for bar, d in zip(bars_e, cohens_d):
    ax_e2.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2, f'{d:.2f}', 
              ha='left', va='center', fontsize=10, fontweight='bold')

ax_e2.spines['top'].set_visible(False)
ax_e2.spines['right'].set_visible(False)

# Panel F: Summary Statistics
ax_f2 = axes[1, 2]
ax_f2.axis('off')
ax_f2.set_title('F) Summary Statistics', fontsize=12, fontweight='bold', loc='left')

summary_text = """
┌─────────────────────────────────────┐
│  SIMULATION PARAMETERS              │
├─────────────────────────────────────┤
│  Total simulation time:    600 ns   │
│  Temperatures:             6        │
│  Frames analyzed:          60,006   │
├─────────────────────────────────────┤
│  MODEL PERFORMANCE                  │
├─────────────────────────────────────┤
│  Training R²:              0.994    │
│  LOTO MAE:                 44.3K    │
│  Variance captured:        78%      │
├─────────────────────────────────────┤
│  KEY FINDING                        │
├─────────────────────────────────────┤
│  Native contacts:          79.6%    │
│  Correlation (r):          −0.887   │
│  Cohen's d:                8.16     │
│  Contact loss (450K):      47.1%    │
└─────────────────────────────────────┘
"""

ax_f2.text(0.1, 0.95, summary_text, transform=ax_f2.transAxes, fontsize=10,
          fontfamily='monospace', verticalalignment='top',
          bbox=dict(boxstyle='round', facecolor='#F8F9FA', edgecolor='gray'))

plt.tight_layout()
plt.savefig(output_dir / 'Figure2_Results.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure2_Results.pdf', bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

print(f"Saved: {output_dir}/Figure2_Results.png")
print(f"Saved: {output_dir}/Figure2_Results.pdf")

# ================================================================
# Create Figure 3: Structural Analysis Figure
# ================================================================
print("\nCreating Structural Analysis Figure...")

fig3 = plt.figure(figsize=(16, 12))

# Load PyMOL images if available
pymol_dir = Path('results/figures/pymol')

# Panel positions
panels = [
    ('A) p53 Structure', 'Fig1A_p53_cartoon.png', [0.02, 0.52, 0.3, 0.42]),
    ('B) Cancer Hotspots', 'Fig1B_p53_hotspots.png', [0.35, 0.52, 0.3, 0.42]),
    ('C) Zinc Coordination', 'Fig1D_p53_zinc_site.png', [0.68, 0.52, 0.3, 0.42]),
    ('D) Native Contact Network', 'Fig6_p53_native_contacts.png', [0.02, 0.02, 0.3, 0.42]),
    ('E) Thermal Sensitivity', 'Fig5_p53_thermal_sensitivity.png', [0.35, 0.02, 0.3, 0.42]),
    ('F) p53 Family Superposition', 'Fig4_p53_p63_p73_superposition.png', [0.68, 0.02, 0.3, 0.42]),
]

for title, filename, pos in panels:
    ax = fig3.add_axes(pos)
    ax.set_title(title, fontsize=12, fontweight='bold', loc='left')
    
    img_path = pymol_dir / filename
    if img_path.exists():
        try:
            img = Image.open(img_path)
            ax.imshow(img)
        except:
            ax.text(0.5, 0.5, f'Error loading\n{filename}', ha='center', va='center',
                   transform=ax.transAxes)
    else:
        ax.text(0.5, 0.5, f'Image not found:\n{filename}', ha='center', va='center',
               transform=ax.transAxes, fontsize=10)
    ax.axis('off')

plt.savefig(output_dir / 'Figure3_Structures.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure3_Structures.pdf', bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

print(f"Saved: {output_dir}/Figure3_Structures.png")
print(f"Saved: {output_dir}/Figure3_Structures.pdf")

# ================================================================
# Create Figure 4: SHAP Detailed Analysis
# ================================================================
print("\nCreating SHAP Analysis Figure...")

fig4, axes = plt.subplots(2, 2, figsize=(14, 12))

# Try to load existing SHAP figures
shap_dir = Path('results/figures/xai/classification')

# Panel A: SHAP Summary
ax_a4 = axes[0, 0]
shap_summary_path = shap_dir / 'shap_summary_regression.png'
if shap_summary_path.exists():
    img = Image.open(shap_summary_path)
    ax_a4.imshow(img)
    ax_a4.set_title('A) SHAP Summary Plot', fontsize=12, fontweight='bold', loc='left')
else:
    ax_a4.text(0.5, 0.5, 'SHAP Summary\n(shap_summary_regression.png)', 
              ha='center', va='center', transform=ax_a4.transAxes)
    ax_a4.set_title('A) SHAP Summary Plot', fontsize=12, fontweight='bold', loc='left')
ax_a4.axis('off')

# Panel B: SHAP Bar
ax_b4 = axes[0, 1]
shap_bar_path = shap_dir / 'shap_importance_regression.png'
if shap_bar_path.exists():
    img = Image.open(shap_bar_path)
    ax_b4.imshow(img)
    ax_b4.set_title('B) Feature Importance', fontsize=12, fontweight='bold', loc='left')
else:
    ax_b4.text(0.5, 0.5, 'SHAP Importance\n(shap_importance_regression.png)', 
              ha='center', va='center', transform=ax_b4.transAxes)
    ax_b4.set_title('B) Feature Importance', fontsize=12, fontweight='bold', loc='left')
ax_b4.axis('off')

# Panel C: Prediction scatter
ax_c4 = axes[1, 0]
scatter_path = shap_dir / 'prediction_scatter.png'
if scatter_path.exists():
    img = Image.open(scatter_path)
    ax_c4.imshow(img)
    ax_c4.set_title('C) Training Performance', fontsize=12, fontweight='bold', loc='left')
else:
    ax_c4.text(0.5, 0.5, 'Prediction Scatter\n(prediction_scatter.png)', 
              ha='center', va='center', transform=ax_c4.transAxes)
    ax_c4.set_title('C) Training Performance', fontsize=12, fontweight='bold', loc='left')
ax_c4.axis('off')

# Panel D: LOTO
ax_d4 = axes[1, 1]
loto_path = shap_dir / 'loto_regression.png'
if loto_path.exists():
    img = Image.open(loto_path)
    ax_d4.imshow(img)
    ax_d4.set_title('D) LOTO Validation', fontsize=12, fontweight='bold', loc='left')
else:
    ax_d4.text(0.5, 0.5, 'LOTO Validation\n(loto_regression.png)', 
              ha='center', va='center', transform=ax_d4.transAxes)
    ax_d4.set_title('D) LOTO Validation', fontsize=12, fontweight='bold', loc='left')
ax_d4.axis('off')

plt.tight_layout()
plt.savefig(output_dir / 'Figure4_SHAP_Analysis.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.savefig(output_dir / 'Figure4_SHAP_Analysis.pdf', bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()

print(f"Saved: {output_dir}/Figure4_SHAP_Analysis.png")
print(f"Saved: {output_dir}/Figure4_SHAP_Analysis.pdf")

# ================================================================
# Summary
# ================================================================
print("\n" + "="*70)
print("PUBLICATION FIGURES COMPLETE")
print("="*70)
print(f"\nGenerated figures in {output_dir}:")
print("  1. Figure1_Introduction.png/pdf - Overview and study design")
print("  2. Figure2_Results.png/pdf - Key results summary")
print("  3. Figure3_Structures.png/pdf - 3D structural analysis")
print("  4. Figure4_SHAP_Analysis.png/pdf - SHAP explainability")
print("\nAll figures are publication-ready (300 DPI)")
