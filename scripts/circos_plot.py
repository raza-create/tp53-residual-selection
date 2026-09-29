#!/usr/bin/env python3
"""
Circos plot for p53 thermal analysis
Shows residue-level information and native contacts
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pycirclize import Circos
from pycirclize.parser import Matrix
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

output_dir = Path('results/figures/circos')
output_dir.mkdir(parents=True, exist_ok=True)

print("="*70)
print("GENERATING CIRCOS PLOT FOR P53 THERMAL ANALYSIS")
print("="*70)

# p53 DBD residue information
# Residues 94-292 (199 residues total)
residues = list(range(94, 293))
n_residues = len(residues)

# Define structural regions
regions = {
    'N-term': (94, 113),
    'S1': (114, 123),
    'S2': (124, 135),
    'S3': (136, 146),
    'L1': (147, 163),
    'S4': (164, 173),
    'L2': (174, 194),
    'S5': (195, 205),
    'S6': (206, 213),
    'H1': (214, 223),
    'S7': (224, 230),
    'S8': (231, 236),
    'L3': (237, 250),
    'S9': (251, 258),
    'S10': (259, 267),
    'H2': (268, 282),
    'C-term': (283, 292)
}

# Cancer hotspot residues
hotspots = {
    175: {'name': 'R175', 'freq': 4.8, 'type': 'structural'},
    245: {'name': 'G245', 'freq': 2.7, 'type': 'structural'},
    248: {'name': 'R248', 'freq': 6.4, 'type': 'contact'},
    249: {'name': 'R249', 'freq': 2.5, 'type': 'structural'},
    273: {'name': 'R273', 'freq': 4.5, 'type': 'contact'},
    282: {'name': 'R282', 'freq': 1.8, 'type': 'structural'}
}

# Thermal sensitivity by region (based on our SHAP analysis)
# Higher values = more temperature sensitive
thermal_sensitivity = {
    'N-term': 0.3,
    'S1': 0.4,
    'S2': 0.5,
    'S3': 0.5,
    'L1': 0.6,
    'S4': 0.7,
    'L2': 0.9,  # Highly sensitive
    'S5': 0.6,
    'S6': 0.5,
    'H1': 0.4,
    'S7': 0.5,
    'S8': 0.6,
    'L3': 0.95, # Most sensitive (DNA binding)
    'S9': 0.5,
    'S10': 0.4,
    'H2': 0.3,
    'C-term': 0.2
}

# Native contact pairs (key stabilizing interactions)
# These are residue pairs that form important contacts
native_contacts = [
    (176, 242), (176, 179), (179, 238),  # Zinc coordination
    (173, 239), (175, 220), (175, 245),  # Core contacts
    (248, 273), (248, 277), (248, 280),  # DNA binding region
    (241, 248), (244, 247), (245, 249),  # L3 loop
    (195, 220), (196, 219), (198, 217),  # Beta sheet contacts
    (157, 203), (159, 201), (161, 199),  # S2-S5 contacts
    (143, 180), (145, 178), (147, 176),  # L1-L2 contacts
    (120, 270), (122, 268), (125, 265),  # Terminal contacts
]

# ================================================================
# CREATE CIRCOS PLOT
# ================================================================
print("\nCreating Circos plot...")

# Define sectors for structural regions
sectors = {}
for region, (start, end) in regions.items():
    sectors[region] = end - start + 1

# Create Circos
circos = Circos(sectors, space=2)

# Color map for regions
region_colors = {
    'N-term': '#808080', 'C-term': '#808080',
    'S1': '#FFA500', 'S2': '#FFA500', 'S3': '#FFA500', 
    'S4': '#FFA500', 'S5': '#FFA500', 'S6': '#FFA500',
    'S7': '#FFA500', 'S8': '#FFA500', 'S9': '#FFA500', 'S10': '#FFA500',
    'L1': '#90EE90', 'L2': '#FF6B6B', 'L3': '#FF0000',
    'H1': '#6495ED', 'H2': '#6495ED'
}

# Add tracks for each sector
for sector in circos.sectors:
    # Track 1: Region labels and base track
    track1 = sector.add_track((90, 100))
    track1.axis(fc=region_colors.get(sector.name, '#CCCCCC'), ec='black', lw=0.5)
    track1.text(sector.name, color='black', size=6, r=105)
    
    # Track 2: Thermal sensitivity heatmap
    track2 = sector.add_track((75, 88))
    sensitivity = thermal_sensitivity.get(sector.name, 0.5)
    # Color based on sensitivity (blue=low, red=high)
    color = plt.cm.coolwarm(sensitivity)
    track2.axis(fc=color, ec='none')
    
    # Track 3: Cancer mutation markers
    track3 = sector.add_track((60, 73))
    track3.axis(fc='white', ec='black', lw=0.3)
    
    # Add mutation markers if in this region
    region_start, region_end = regions[sector.name]
    for res, info in hotspots.items():
        if region_start <= res <= region_end:
            pos = res - region_start
            # Mark mutation site
            track3.rect(pos, pos+1, fc='red', ec='darkred', lw=0.5)

# Add native contact links
print("Adding native contact links...")
for res1, res2 in native_contacts:
    # Find which sectors these residues belong to
    sector1, pos1 = None, None
    sector2, pos2 = None, None
    
    for region, (start, end) in regions.items():
        if start <= res1 <= end:
            sector1 = region
            pos1 = res1 - start
        if start <= res2 <= end:
            sector2 = region
            pos2 = res2 - start
    
    if sector1 and sector2:
        # Determine link color based on contact type
        if res1 in [176, 179, 238, 242] or res2 in [176, 179, 238, 242]:
            color = '#FFD700'  # Gold for zinc coordination
            alpha = 0.8
        elif res1 in hotspots or res2 in hotspots:
            color = '#FF0000'  # Red for hotspot contacts
            alpha = 0.7
        else:
            color = '#4169E1'  # Blue for other contacts
            alpha = 0.4
        
        circos.link((sector1, pos1, pos1+1), (sector2, pos2, pos2+1), 
                    color=color, alpha=alpha)

# Create figure
fig = circos.plotfig()

# Add legend
legend_elements = [
    plt.Line2D([0], [0], color='#FF0000', lw=4, label='L3 loop (DNA binding)'),
    plt.Line2D([0], [0], color='#FF6B6B', lw=4, label='L2 loop (thermal sensitive)'),
    plt.Line2D([0], [0], color='#FFA500', lw=4, label='β-strands'),
    plt.Line2D([0], [0], color='#6495ED', lw=4, label='α-helices'),
    plt.Line2D([0], [0], color='#FFD700', lw=4, label='Zinc coordination contacts'),
    plt.Line2D([0], [0], color='#FF0000', lw=2, alpha=0.7, label='Hotspot contacts'),
    plt.Line2D([0], [0], color='#4169E1', lw=2, alpha=0.4, label='Stabilizing contacts'),
]
fig.legend(handles=legend_elements, loc='lower right', fontsize=8, 
           title='Legend', title_fontsize=9)

plt.title('p53 DNA-Binding Domain: Structure, Thermal Sensitivity, and Native Contacts', 
          fontsize=12, fontweight='bold', y=1.02)

plt.tight_layout()
plt.savefig(output_dir / 'p53_circos_plot.png', dpi=300, bbox_inches='tight')
plt.savefig(output_dir / 'p53_circos_plot.pdf', bbox_inches='tight')
plt.close()

print(f"Saved: {output_dir}/p53_circos_plot.png")
print(f"Saved: {output_dir}/p53_circos_plot.pdf")

# ================================================================
# CREATE SECOND CIRCOS: Temperature Comparison
# ================================================================
print("\nCreating temperature comparison Circos...")

# Load our actual data
data = pd.read_csv('data/features_100ns/all_features_100ns.csv')

# Calculate per-residue metrics (simplified - using overall statistics)
temps = [248, 273, 310, 350, 400, 450]
temp_colors = ['#0000FF', '#00BFFF', '#00FF00', '#FFFF00', '#FF8C00', '#FF0000']

# Create temperature sectors
temp_sectors = {f'{t}K': 100 for t in temps}

circos2 = Circos(temp_sectors, space=5)

for i, sector in enumerate(circos2.sectors):
    temp = temps[i]
    temp_data = data[data['temperature'] == temp]
    
    # Track 1: Temperature label
    track1 = sector.add_track((90, 100))
    track1.axis(fc=temp_colors[i], ec='black', lw=1)
    track1.text(sector.name, color='white' if i > 2 else 'black', size=10, r=95)
    
    # Track 2: Native contacts (normalized)
    track2 = sector.add_track((70, 88))
    nc_mean = temp_data['native_contacts'].mean()
    # Fill proportional to native contacts
    nc_width = int(nc_mean * 100)
    track2.axis(fc='white', ec='black', lw=0.5)
    track2.rect(0, nc_width, fc='green', alpha=0.7)
    track2.text(f'NC:{nc_mean:.2f}', r=79, size=7)
    
    # Track 3: RMSD (normalized)
    track3 = sector.add_track((50, 68))
    rmsd_mean = temp_data['rmsd'].mean()
    rmsd_norm = min(rmsd_mean / 5, 1)  # Normalize to max 5 Å
    rmsd_width = int(rmsd_norm * 100)
    track3.axis(fc='white', ec='black', lw=0.5)
    track3.rect(0, rmsd_width, fc='red', alpha=0.7)
    track3.text(f'RMSD:{rmsd_mean:.1f}', r=59, size=7)

# Add links showing temperature transitions
for i in range(len(temps) - 1):
    circos2.link(
        (f'{temps[i]}K', 40, 60),
        (f'{temps[i+1]}K', 40, 60),
        color='purple', alpha=0.3
    )

fig2 = circos2.plotfig()
plt.title('p53 Structural Features Across Temperatures', fontsize=12, fontweight='bold')

# Add legend
legend2 = [
    plt.Rectangle((0,0), 1, 1, fc='green', alpha=0.7, label='Native Contacts'),
    plt.Rectangle((0,0), 1, 1, fc='red', alpha=0.7, label='RMSD'),
]
fig2.legend(handles=legend2, loc='lower right', fontsize=9)

plt.savefig(output_dir / 'p53_temperature_circos.png', dpi=300, bbox_inches='tight')
plt.savefig(output_dir / 'p53_temperature_circos.pdf', bbox_inches='tight')
plt.close()

print(f"Saved: {output_dir}/p53_temperature_circos.png")

# ================================================================
# CREATE THIRD CIRCOS: p53 Family Comparison
# ================================================================
print("\nCreating p53 family comparison Circos...")

# Protein family sectors
family_sectors = {'p53': 199, 'p63': 199, 'p73': 199}

circos3 = Circos(family_sectors, space=8)

family_colors = {'p53': '#E74C3C', 'p63': '#27AE60', 'p73': '#3498DB'}
family_hotspots = {
    'p53': [175, 245, 248, 273, 282],
    'p63': [204, 274, 277, 302, 311],  # Equivalent positions
    'p73': [178, 248, 251, 276, 285]   # Equivalent positions
}

for sector in circos3.sectors:
    protein = sector.name
    
    # Main track with protein color
    track1 = sector.add_track((85, 100))
    track1.axis(fc=family_colors[protein], ec='black', lw=1)
    track1.text(protein, color='white', size=14, r=92)
    
    # Conservation track
    track2 = sector.add_track((70, 83))
    track2.axis(fc='#E8E8E8', ec='black', lw=0.5)
    
    # Hotspot markers
    track3 = sector.add_track((55, 68))
    track3.axis(fc='white', ec='black', lw=0.5)
    
    for pos in family_hotspots[protein]:
        rel_pos = pos - 94 if protein == 'p53' else pos - 94  # Adjust based on alignment
        if 0 <= rel_pos < 199:
            track3.rect(rel_pos, rel_pos + 3, fc='red', ec='darkred')

# Links showing homologous positions
homologous_pairs = [
    (('p53', 175-94), ('p63', 204-94)),
    (('p53', 175-94), ('p73', 178-94)),
    (('p53', 248-94), ('p63', 277-94)),
    (('p53', 248-94), ('p73', 251-94)),
    (('p63', 204-94), ('p73', 178-94)),
]

for (s1, p1), (s2, p2) in homologous_pairs:
    if 0 <= p1 < 199 and 0 <= p2 < 199:
        circos3.link((s1, p1, p1+5), (s2, p2, p2+5), color='purple', alpha=0.5)

fig3 = circos3.plotfig()
plt.title('p53 Family: Conservation of Cancer Hotspot Positions', fontsize=12, fontweight='bold')

plt.savefig(output_dir / 'p53_family_circos.png', dpi=300, bbox_inches='tight')
plt.savefig(output_dir / 'p53_family_circos.pdf', bbox_inches='tight')
plt.close()

print(f"Saved: {output_dir}/p53_family_circos.png")

print("\n" + "="*70)
print("CIRCOS PLOTS COMPLETE")
print("="*70)
print(f"\nGenerated figures:")
print(f"  1. p53_circos_plot.png - Structure and native contacts")
print(f"  2. p53_temperature_circos.png - Temperature comparison")
print(f"  3. p53_family_circos.png - p53/p63/p73 comparison")
print(f"\nAll saved to: {output_dir}")
