#!/usr/bin/env python3
"""
PyMOL script to create publication-quality protein structure images
"""

import subprocess
from pathlib import Path

# Paths
pdb_dir = Path("/home/raza/p53_thermal_xai/results/figures/pdb_files")
output_dir = Path("/home/raza/p53_thermal_xai/results/figures/pymol_images")
output_dir.mkdir(parents=True, exist_ok=True)

# Create PyMOL script for individual structures
pymol_script = '''
# PyMOL Script for p53 Visualization

# Settings for high quality
set ray_trace_mode, 1
set ray_shadows, 1
set antialias, 2
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_flat_sheets, 1
set cartoon_side_chain_helper, 1
set depth_cue, 1
set spec_reflect, 0.5
set spec_power, 200
bg_color white

# Function to visualize a structure
'''

temperatures = [248, 273, 310, 350, 400, 450, 500]
temp_labels = ['m25C', '0C', '37C', '77C', '127C', '177C', '227C']

# Create individual structure scripts
for i, temp in enumerate(temperatures):
    pdb_file = pdb_dir / f"p53_{temp}K_last.pdb"
    
    if not pdb_file.exists():
        print(f"Warning: {pdb_file} not found")
        continue
    
    script_content = f'''
# PyMOL Script for p53 at {temp}K
reinitialize

# High quality settings
set ray_trace_mode, 1
set ray_shadows, 1
set antialias, 2
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_flat_sheets, 1
set cartoon_side_chain_helper, 1
set depth_cue, 1
set spec_reflect, 0.5
set spec_power, 200
set cartoon_transparency, 0
bg_color white

# Load structure
load {pdb_file}, p53_{temp}K

# Show as cartoon
hide everything
show cartoon, p53_{temp}K

# Color by secondary structure
color marine, ss h  # helices - blue
color magenta, ss s  # sheets - magenta
color cyan, ss l+''  # loops - cyan

# Highlight cancer hotspot residues as spheres
select hotspots, resi 175+245+248+249+273+282
show spheres, hotspots and name CA
color red, hotspots
set sphere_scale, 0.8, hotspots

# Add labels to hotspots
label hotspots and name CA, "%s%s" % (resn, resi)
set label_color, black
set label_size, 14
set label_font_id, 7

# Center and orient
center p53_{temp}K
orient p53_{temp}K
zoom p53_{temp}K, 5

# Save front view
ray 1200, 1000
png {output_dir}/p53_{temp}K_front.png, dpi=300

# Rotate for side view
rotate y, 90
ray 1200, 1000
png {output_dir}/p53_{temp}K_side.png, dpi=300

# Top view
rotate x, 90
ray 1200, 1000
png {output_dir}/p53_{temp}K_top.png, dpi=300

quit
'''
    
    script_file = output_dir / f"script_{temp}K.pml"
    with open(script_file, 'w') as f:
        f.write(script_content)
    
    print(f"Created: script_{temp}K.pml")

# Create comparison script (cold vs hot overlay)
comparison_script = f'''
# PyMOL Script for Cold vs Hot Comparison
reinitialize

# High quality settings
set ray_trace_mode, 1
set ray_shadows, 1
set antialias, 2
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set depth_cue, 1
bg_color white

# Load cold structure (248K)
load {pdb_dir}/p53_248K_last.pdb, cold_248K

# Load hot structure (500K)
load {pdb_dir}/p53_500K_last.pdb, hot_500K

# Align structures
align hot_500K, cold_248K

# Show as cartoon
hide everything
show cartoon, cold_248K
show cartoon, hot_500K

# Color: Cold = Blue, Hot = Red
color blue, cold_248K
color red, hot_500K

# Set transparency
set cartoon_transparency, 0.3, hot_500K

# Hotspots on cold structure
select hotspots_cold, cold_248K and resi 175+245+248+249+273+282
show spheres, hotspots_cold and name CA
color limegreen, hotspots_cold
set sphere_scale, 1.0, hotspots_cold

# Center and orient
center cold_248K
orient cold_248K
zoom all, 5

# Save
ray 1400, 1000
png {output_dir}/p53_cold_vs_hot_overlay.png, dpi=300

# Side view
rotate y, 90
ray 1400, 1000
png {output_dir}/p53_cold_vs_hot_side.png, dpi=300

quit
'''

with open(output_dir / "script_comparison.pml", 'w') as f:
    f.write(comparison_script)
print("Created: script_comparison.pml")

# Create all temperatures overlay script
all_temps_script = f'''
# PyMOL Script for All Temperatures Overlay
reinitialize

# High quality settings
set ray_trace_mode, 1
set ray_shadows, 1
set antialias, 2
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set depth_cue, 1
bg_color white

# Load all structures
load {pdb_dir}/p53_248K_last.pdb, temp_248K
load {pdb_dir}/p53_273K_last.pdb, temp_273K
load {pdb_dir}/p53_310K_last.pdb, temp_310K
load {pdb_dir}/p53_350K_last.pdb, temp_350K
load {pdb_dir}/p53_400K_last.pdb, temp_400K
load {pdb_dir}/p53_450K_last.pdb, temp_450K
load {pdb_dir}/p53_500K_last.pdb, temp_500K

# Align all to reference (248K)
align temp_273K, temp_248K
align temp_310K, temp_248K
align temp_350K, temp_248K
align temp_400K, temp_248K
align temp_450K, temp_248K
align temp_500K, temp_248K

# Show as cartoon
hide everything
show cartoon, all

# Color gradient: cold (blue) to hot (red)
color blue, temp_248K
color deepteal, temp_273K
color green, temp_310K
color yellow, temp_350K
color orange, temp_400K
color firebrick, temp_450K
color red, temp_500K

# Set transparency for better visibility
set cartoon_transparency, 0.2, all

# Center and orient
center temp_248K
orient temp_248K
zoom all, 5

# Save
ray 1400, 1000
png {output_dir}/p53_all_temps_overlay.png, dpi=300

# Side view
rotate y, 90
ray 1400, 1000
png {output_dir}/p53_all_temps_overlay_side.png, dpi=300

quit
'''

with open(output_dir / "script_all_temps.pml", 'w') as f:
    f.write(all_temps_script)
print("Created: script_all_temps.pml")

# Create surface representation script (like image 2 panel C)
surface_script = f'''
# PyMOL Script for Surface Representation
reinitialize

# High quality settings
set ray_trace_mode, 1
set ray_shadows, 1
set antialias, 2
bg_color white

# Load structure
load {pdb_dir}/p53_310K_last.pdb, p53

# Show surface
hide everything
show surface, p53

# Color by conservation/property - use B-factor or residue type
# Color hydrophobic residues
select hydrophobic, resn ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO
select polar, resn SER+THR+ASN+GLN+TYR+CYS
select charged_pos, resn LYS+ARG+HIS
select charged_neg, resn ASP+GLU

color white, p53
color cyan, hydrophobic
color magenta, polar
color blue, charged_pos
color red, charged_neg

# Highlight hotspots
select hotspots, resi 175+245+248+249+273+282
color yellow, hotspots

# Center and orient
center p53
orient p53
zoom p53, 5

# Save front view
ray 1200, 1000
png {output_dir}/p53_surface_front.png, dpi=300

# Rotate 180 for back view
rotate y, 180
ray 1200, 1000
png {output_dir}/p53_surface_back.png, dpi=300

quit
'''

with open(output_dir / "script_surface.pml", 'w') as f:
    f.write(surface_script)
print("Created: script_surface.pml")

print(f"\nAll scripts created in: {output_dir}")
print("\nNow run PyMOL with these scripts...")
