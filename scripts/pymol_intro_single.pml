# ================================================================
# SINGLE COMPREHENSIVE INTRODUCTION FIGURE
# All features in ONE panel for high-impact journal
# ================================================================

# Fetch and prepare structure
fetch 2XWR, p53
remove solvent
remove not chain A

# ================================================================
# HIGHEST QUALITY SETTINGS
# ================================================================
bg_color white
set ray_shadows, 1
set ray_shadow_decay_factor, 0.1
set antialias, 3
set ray_trace_mode, 1
set ambient, 0.35
set spec_reflect, 0.6
set shininess, 20
set specular, 1
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_cylindrical_helices, 1
set cartoon_tube_radius, 0.35
set cartoon_loop_radius, 0.25
set stick_radius, 0.2
set sphere_scale, 0.25
set depth_cue, 1
set fog_start, 0.55
set label_font_id, 7
set label_size, 28
set label_color, black
set label_outline_color, white
set label_shadow_mode, 3

# ================================================================
# BASE STRUCTURE - Secondary Structure Coloring
# ================================================================
hide all
show cartoon, p53

# Color secondary structure elements
# Alpha helices - blue
color tv_blue, ss h

# Beta sheets - orange/gold
color tv_orange, ss s

# Loops - light gray
color gray70, ss l+''

# ================================================================
# THERMAL SENSITIVITY COLORING ON LOOPS
# ================================================================
# L2 loop (highly thermally sensitive) - salmon/red
select L2_loop, resi 164-194
color salmon, L2_loop

# L3 loop (most sensitive - DNA binding) - red
select L3_loop, resi 237-250
color firebrick, L3_loop

# ================================================================
# CANCER HOTSPOTS - Spheres with Labels
# ================================================================
# Define all hotspots
select R175, resi 175
select G245, resi 245  
select R248, resi 248
select R249, resi 249
select R273, resi 273
select R282, resi 282
select all_hotspots, R175 or G245 or R248 or R249 or R273 or R282

# Show as spheres (CA only for clean look)
show spheres, all_hotspots and name CA
set sphere_scale, 1.4, all_hotspots

# Color by mutation type
# Structural mutations - RED
color red, R175 and name CA
color red, G245 and name CA
color red, R249 and name CA
color red, R282 and name CA

# DNA Contact mutations - BLUE
color blue, R248 and name CA
color blue, R273 and name CA

# Add labels
label R175 and name CA, "R175H"
label G245 and name CA, "G245S"
label R248 and name CA, "R248Q"
label R273 and name CA, "R273H"
label R282 and name CA, "R282W"

# ================================================================
# ZINC COORDINATION SITE
# ================================================================
# Show zinc as large sphere
select zinc, name ZN
show spheres, zinc
color gray40, zinc
set sphere_scale, 1.2, zinc

# Show coordinating residues as sticks
select zinc_coord, resi 176+179+238+242
show sticks, zinc_coord and (name CA+CB+CG+CD+CE+NE2+SG+ND1 or name C+N+O)
color orange, zinc_coord and not name C+N+O
color gray50, zinc_coord and name C+N+O

# ================================================================
# DNA BINDING RESIDUES - Show key contacts
# ================================================================
select dna_bind, resi 248+273+280+281+277
show sticks, dna_bind and sidechain and not all_hotspots
color marine, dna_bind and sidechain and not all_hotspots

# ================================================================
# OPTIMAL ORIENTATION - Show all features clearly
# ================================================================
reset
orient p53

# Rotate to best view showing:
# - Hotspots visible
# - Zinc site visible
# - L2/L3 loops visible
# - DNA binding surface partially visible
turn y, 35
turn x, -15
turn z, 5

zoom p53, -2

# ================================================================
# RENDER HIGH QUALITY
# ================================================================
ray 3600, 3600
png /home/raza/p53_thermal_xai/results/figures/publication/Figure1_Introduction_3D.png, dpi=300

# Also save in different format
ray 3600, 3600
png /home/raza/p53_thermal_xai/results/figures/publication/Figure1_Introduction_3D_hires.png, dpi=600

quit
