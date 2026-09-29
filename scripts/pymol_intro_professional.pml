# ================================================================
# HIGH-IMPACT JOURNAL INTRODUCTION FIGURE
# Professional quality with no label overlap
# ================================================================

# Fetch and prepare
fetch 2XWR, p53
remove solvent
remove not chain A

# ================================================================
# PROFESSIONAL RENDERING SETTINGS
# ================================================================
bg_color white
set ray_shadows, 1
set ray_shadow_decay_factor, 0.15
set ray_shadow_decay_range, 5
set antialias, 3
set ray_trace_mode, 1
set ray_trace_gain, 0.015
set ambient, 0.45
set direct, 0.4
set spec_reflect, 0.5
set spec_power, 200
set shininess, 25
set specular, 1
set reflect, 0.3
set depth_cue, 1
set fog_start, 0.6

# Cartoon settings
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_cylindrical_helices, 1
set cartoon_tube_radius, 0.4
set cartoon_loop_radius, 0.3
set cartoon_rect_length, 1.4
set cartoon_oval_length, 1.3

# ================================================================
# STRUCTURE VISUALIZATION
# ================================================================
hide all
show cartoon, p53

# Elegant color scheme - Nature style
# Beta sheets - warm gold
select beta_sheets, ss s
color 0xD4A84B, beta_sheets

# Alpha helices - deep teal  
select alpha_helices, ss h
color 0x2C7873, alpha_helices

# Core loops - light warm gray
color 0xB8B8B8, ss l+''

# ================================================================
# THERMALLY SENSITIVE REGIONS - Gradient coloring
# ================================================================
# L2 loop (thermal sensitive) - coral
select L2, resi 164-194
color 0xE07B54, L2

# L3 loop (most sensitive/DNA binding) - deep red
select L3, resi 237-250
color 0xC0392B, L3

# ================================================================
# ZINC COORDINATION - Elegant display
# ================================================================
select zinc_ion, name ZN
show spheres, zinc_ion
color 0x7D7D7D, zinc_ion
set sphere_scale, 1.0, zinc_ion

# Zinc coordinating residues - orange sticks
select zinc_ligands, resi 176+179+238+242 and (name SG or name NE2 or name ND1)
show sticks, zinc_ligands
set stick_radius, 0.18, zinc_ligands
color 0xE67E22, zinc_ligands

# Show coordination bonds
distance zinc_bonds, zinc_ion, zinc_ligands, cutoff=3.0, mode=0
hide labels, zinc_bonds
color 0xF39C12, zinc_bonds
set dash_width, 3, zinc_bonds
set dash_gap, 0.15, zinc_bonds

# ================================================================
# CANCER HOTSPOTS - Professional spheres
# ================================================================
# Create selections
select mut_R175, resi 175 and name CA
select mut_G245, resi 245 and name CA
select mut_R248, resi 248 and name CA
select mut_R249, resi 249 and name CA
select mut_R273, resi 273 and name CA
select mut_R282, resi 282 and name CA

# Structural mutations - Crimson red gradient
show spheres, mut_R175
show spheres, mut_G245
show spheres, mut_R249
show spheres, mut_R282
color 0xC0392B, mut_R175
color 0xE74C3C, mut_G245
color 0xC0392B, mut_R249
color 0xE74C3C, mut_R282
set sphere_scale, 1.3, mut_R175 or mut_G245 or mut_R249 or mut_R282

# Contact mutations - Deep blue
show spheres, mut_R248
show spheres, mut_R273
color 0x2980B9, mut_R248
color 0x3498DB, mut_R273
set sphere_scale, 1.3, mut_R248 or mut_R273

# ================================================================
# OPTIMAL VIEWING ANGLE
# ================================================================
reset
orient p53
turn y, 145
turn x, -20
turn z, -5
zoom p53, -5

# ================================================================
# CREATE LABEL OBJECTS (pseudo-atoms for clean labels)
# ================================================================
# Create pseudoatoms positioned away from structure for labels

# R175H label position
pseudoatom label_R175, pos=[-8, 25, 15]
label label_R175, "R175H"

# G245S label position  
pseudoatom label_G245, pos=[18, -5, -20]
label label_G245, "G245S"

# R248Q label position
pseudoatom label_R248, pos=[25, 10, -5]
label label_R248, "R248Q"

# R249S label position
pseudoatom label_R249, pos=[20, -15, -15]
label label_R249, "R249S"

# R273H label position
pseudoatom label_R273, pos=[10, 20, 10]
label label_R273, "R273H"

# R282W label position
pseudoatom label_R282, pos=[-5, 18, 20]
label label_R282, "R282W"

# Zinc label
pseudoatom label_zinc, pos=[-15, 5, -25]
label label_zinc, "Zn2+"

# Loop labels
pseudoatom label_L2, pos=[-25, 0, 10]
label label_L2, "L2 Loop"

pseudoatom label_L3, pos=[5, -25, -10]
label label_L3, "L3 Loop"

# Hide pseudoatom spheres
hide spheres, label_*

# Label styling
set label_font_id, 7
set label_size, 24
set label_color, black
set label_bg_color, white
set label_bg_transparency, 0.3
set label_connector, 1
set label_connector_color, gray50
set label_connector_width, 2

# Draw connector lines from labels to residues
distance line_R175, label_R175, mut_R175, mode=0
distance line_G245, label_G245, mut_G245, mode=0
distance line_R248, label_R248, mut_R248, mode=0
distance line_R249, label_R249, mut_R249, mode=0
distance line_R273, label_R273, mut_R273, mode=0
distance line_R282, label_R282, mut_R282, mode=0

# Style connector lines
hide labels, line_*
set dash_color, gray40, line_*
set dash_width, 1.5, line_*
set dash_gap, 0.0, line_*
set dash_length, 0.3, line_*

# ================================================================
# FINAL ADJUSTMENTS
# ================================================================
# Fine-tune view
turn y, -5
zoom p53, -3

# ================================================================
# HIGH-RESOLUTION RENDER
# ================================================================
# First render - standard
ray 4000, 4000
png /home/raza/p53_thermal_xai/results/figures/publication/Figure1_Intro_Professional.png, dpi=300

quit
