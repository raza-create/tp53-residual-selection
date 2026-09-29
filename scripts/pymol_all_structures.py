#!/usr/bin/env python3
"""
Generate 3D structure figures for p53, p63, p73 at all temperatures
"""
import subprocess
from pathlib import Path
import os

output_dir = Path('results/figures/pymol')
output_dir.mkdir(parents=True, exist_ok=True)

# Check for trajectory files
sim_dir = Path('data/structures/long_simulations')
temps = [248, 273, 310, 350, 400, 450]

print("="*70)
print("GENERATING 3D STRUCTURES FOR ALL TEMPERATURES")
print("="*70)

# First, extract representative frames from each temperature
extract_script = """
#!/bin/bash
cd /home/raza/p53_thermal_xai/data/structures/long_simulations

for TEMP in 248 273 310 350 400 450; do
    if [ -f "md_100ns_${TEMP}K.xtc" ] || [ -f "md_${TEMP}K.xtc" ]; then
        echo "Extracting frame from ${TEMP}K..."
        # Extract last frame (most equilibrated)
        if [ -f "md_100ns_${TEMP}K.xtc" ]; then
            echo 1 | gmx trjconv -f md_100ns_${TEMP}K.xtc -s md_100ns_${TEMP}K.tpr -o ../extracted/p53_${TEMP}K.pdb -dump -1 2>/dev/null
        elif [ -f "md_${TEMP}K.xtc" ]; then
            echo 1 | gmx trjconv -f md_${TEMP}K.xtc -s *.tpr -o ../extracted/p53_${TEMP}K.pdb -dump -1 2>/dev/null
        fi
    fi
done
"""

# Create extraction directory
os.makedirs('data/structures/extracted', exist_ok=True)

# PyMOL script for comprehensive visualization
pymol_script = """
# ================================================================
# P53 THERMAL ANALYSIS - 3D STRUCTURE FIGURES
# ================================================================

# Settings for publication quality
bg_color white
set ray_shadows, 0
set antialias, 2
set ray_trace_mode, 1
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set stick_radius, 0.2
set sphere_scale, 0.3

# ================================================================
# PART 1: REFERENCE STRUCTURE (2XWR - p53 DBD)
# ================================================================
fetch 2XWR, p53_ref
remove solvent
remove not chain A
hide all

# Figure 1A: Overall cartoon structure
show cartoon, p53_ref
color slate, ss h
color orange, ss s  
color gray80, ss l+''
orient p53_ref
zoom p53_ref, 5
ray 2400, 1800
png {output}/Fig1A_p53_cartoon.png, dpi=300

# Figure 1B: Cancer hotspot mutations
color gray80, p53_ref
select hotspots, p53_ref and resi 175+245+248+249+273+282
color red, hotspots
show spheres, hotspots and name CA
set sphere_scale, 1.2, hotspots
ray 2400, 1800
png {output}/Fig1B_p53_hotspots.png, dpi=300

# Figure 1C: DNA binding residues
hide spheres
color gray80, p53_ref
select dna_bind, p53_ref and resi 248+273+276+277+280+281+120+277
color blue, dna_bind
show sticks, dna_bind
ray 2400, 1800
png {output}/Fig1C_p53_dna_binding.png, dpi=300

# Figure 1D: Zinc coordination site
hide sticks
color gray80, p53_ref
select zinc_coord, p53_ref and resi 176+179+238+242
color orange, zinc_coord
show sticks, zinc_coord
select zn, p53_ref and name ZN
show spheres, zn
color gray50, zn
set sphere_scale, 0.7, zn
ray 2400, 1800
png {output}/Fig1D_p53_zinc_site.png, dpi=300

# Figure 1E: Surface with hotspots
hide all
show surface, p53_ref
color white, p53_ref
set transparency, 0.6
show spheres, hotspots and name CA
color red, hotspots
set sphere_scale, 1.0, hotspots
ray 2400, 1800
png {output}/Fig1E_p53_surface.png, dpi=300

# ================================================================
# PART 2: P63 STRUCTURE (3QYN)
# ================================================================
delete all
fetch 3QYN, p63_ref
remove solvent
remove not chain A
hide all

# Figure 2A: p63 cartoon
show cartoon, p63_ref
color forest, ss h
color tv_orange, ss s
color gray80, ss l+''
orient p63_ref
zoom p63_ref, 5
ray 2400, 1800
png {output}/Fig2A_p63_cartoon.png, dpi=300

# Figure 2B: p63 equivalent hotspots
color gray80, p63_ref
# Equivalent positions to p53 hotspots
select p63_hotspots, p63_ref and resi 204+274+277+278+302+311
color magenta, p63_hotspots
show spheres, p63_hotspots and name CA
set sphere_scale, 1.2, p63_hotspots
ray 2400, 1800
png {output}/Fig2B_p63_hotspots.png, dpi=300

# ================================================================
# PART 3: P73 STRUCTURE (2XWC)
# ================================================================
delete all
fetch 2XWC, p73_ref
remove solvent
remove not chain A
hide all

# Figure 3A: p73 cartoon
show cartoon, p73_ref
color tv_blue, ss h
color salmon, ss s
color gray80, ss l+''
orient p73_ref
zoom p73_ref, 5
ray 2400, 1800
png {output}/Fig3A_p73_cartoon.png, dpi=300

# Figure 3B: p73 equivalent hotspots
color gray80, p73_ref
select p73_hotspots, p73_ref and resi 178+248+251+252+276+285
color cyan, p73_hotspots
show spheres, p73_hotspots and name CA
set sphere_scale, 1.2, p73_hotspots
ray 2400, 1800
png {output}/Fig3B_p73_hotspots.png, dpi=300

# ================================================================
# PART 4: SUPERPOSITION OF ALL THREE
# ================================================================
delete all
fetch 2XWR, p53
fetch 3QYN, p63
fetch 2XWC, p73
remove solvent
remove (not chain A)

# Align structures
align p63, p53
align p73, p53

hide all
show cartoon, all
color slate, p53
color forest, p63
color tv_blue, p73
orient all
zoom all, 5
ray 2400, 1800
png {output}/Fig4_p53_p63_p73_superposition.png, dpi=300

# ================================================================
# PART 5: TEMPERATURE COMPARISON (if structures available)
# ================================================================
delete all

# Try to load extracted structures from MD
python
import os
temps = [248, 273, 310, 350, 400, 450]
extracted_dir = "/home/raza/p53_thermal_xai/data/structures/extracted"
for t in temps:
    pdb_file = f"{extracted_dir}/p53_{t}K.pdb"
    if os.path.exists(pdb_file):
        cmd.load(pdb_file, f"p53_{t}K")
        print(f"Loaded: p53_{t}K")
python end

# If structures loaded, create temperature comparison
hide all
show cartoon, p53_*
# Color by temperature (blue=cold, red=hot)
color blue, p53_248K
color cyan, p53_273K
color green, p53_310K
color yellow, p53_350K
color orange, p53_400K
color red, p53_450K
orient
zoom all, 5
ray 2400, 1800
png {output}/Fig5_temperature_comparison.png, dpi=300

quit
""".format(output=output_dir)

# Save PyMOL script
script_path = Path('scripts/pymol_all_structures.pml')
with open(script_path, 'w') as f:
    f.write(pymol_script)

print(f"PyMOL script saved: {script_path}")
print("\nThis script will generate:")
print("  - Fig1A-E: p53 reference structure (cartoon, hotspots, DNA binding, zinc, surface)")
print("  - Fig2A-B: p63 structure and hotspots")
print("  - Fig3A-B: p73 structure and hotspots")
print("  - Fig4: Superposition of p53/p63/p73")
print("  - Fig5: Temperature comparison (if MD structures available)")
print("\nRun with: pymol -cq scripts/pymol_all_structures.pml")
