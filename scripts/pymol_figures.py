#!/usr/bin/env python3
"""
Generate PyMOL 3D structure figures for publication
"""
import subprocess
from pathlib import Path

output_dir = Path('results/figures/pymol')
output_dir.mkdir(parents=True, exist_ok=True)

# PyMOL script for p53 structure visualization
pymol_script = """
# Load p53 DNA-binding domain structure
fetch 2XWR, p53
remove solvent
remove not chain A

# Set background and display
bg_color white
set ray_shadows, 0
set antialias, 2
set ray_trace_mode, 1

# === FIGURE 1: Overall structure with secondary structure ===
hide all
show cartoon, p53
color slate, ss h
color orange, ss s
color white, ss l+''
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
orient p53
zoom p53, 5
ray 2400, 2400
png {output}/p53_structure_cartoon.png, dpi=300

# === FIGURE 2: Structure colored by B-factor (flexibility) ===
spectrum b, blue_white_red, p53
ray 2400, 2400
png {output}/p53_bfactor.png, dpi=300

# === FIGURE 3: Cancer mutation hotspots ===
color gray80, p53
# Hotspot residues (using 2XWR numbering)
select hotspots, resi 175+245+248+249+273+282
color red, hotspots
show spheres, hotspots and name CA
set sphere_scale, 1.5, hotspots
ray 2400, 2400
png {output}/p53_hotspots.png, dpi=300

# === FIGURE 4: DNA binding interface ===
color gray80, p53
# DNA contact residues
select dna_contact, resi 248+273+276+277+280+281
color blue, dna_contact
show sticks, dna_contact
ray 2400, 2400
png {output}/p53_dna_contact.png, dpi=300

# === FIGURE 5: Zinc binding site ===
color gray80, p53
select zinc_site, resi 176+179+238+242
color orange, zinc_site
show sticks, zinc_site
# Show zinc if present
select zinc, name ZN
show spheres, zinc
color gray50, zinc
set sphere_scale, 0.8, zinc
ray 2400, 2400
png {output}/p53_zinc_site.png, dpi=300

# === FIGURE 6: Surface representation with hotspots ===
hide all
show surface, p53
color white, p53
set transparency, 0.7
show spheres, hotspots and name CA
color red, hotspots
set sphere_scale, 1.2, hotspots
ray 2400, 2400
png {output}/p53_surface_hotspots.png, dpi=300

# === FIGURE 7: Cartoon with native contacts network (conceptual) ===
hide all
show cartoon, p53
color slate, p53
# Show key stabilizing residues
select core_residues, resi 157+158+159+160+161+203+204+205+206+207
color green, core_residues
show sticks, core_residues
ray 2400, 2400
png {output}/p53_core_contacts.png, dpi=300

quit
""".format(output=output_dir)

# Save PyMOL script
script_path = Path('scripts/pymol_visualization.pml')
with open(script_path, 'w') as f:
    f.write(pymol_script)

print("="*70)
print("PYMOL SCRIPT CREATED")
print("="*70)
print(f"Script saved: {script_path}")
print("\nTo generate figures, run:")
print(f"  pymol -cq {script_path}")
print("\nOr if PyMOL is not installed:")
print("  sudo apt install pymol")
print("  pymol -cq scripts/pymol_visualization.pml")
