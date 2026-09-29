# ================================================================
# HIGH-IMPACT JOURNAL FIGURE - FULLY ANNOTATED
# Clear domain, region, and structural element labels
# ================================================================

fetch 2XWR, p53
remove solvent
remove not chain A

# ================================================================
# PROFESSIONAL RENDERING SETTINGS
# ================================================================
bg_color white
set ray_shadows, 1
set ray_shadow_decay_factor, 0.15
set antialias, 3
set ray_trace_mode, 1
set ambient, 0.45
set direct, 0.4
set spec_reflect, 0.5
set shininess, 25
set specular, 1
set depth_cue, 1
set fog_start, 0.6

set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_cylindrical_helices, 1
set cartoon_tube_radius, 0.4
set cartoon_loop_radius, 0.3

# ================================================================
# COLOR SCHEME BY STRUCTURAL REGIONS
# ================================================================
hide all
show cartoon, p53

# Define all structural regions of p53 DBD (94-292)
# Beta sandwich core
select S1, resi 112-117
select S2, resi 128-133
select S3, resi 141-143
select S4, resi 156-158
select S5, resi 199-204
select S6, resi 207-212
select S7, resi 229-232
select S8, resi 253-258
select S9, resi 264-269
select S10, resi 272-274

# Alpha helices
select H1, resi 177-181
select H2, resi 278-286

# Loops
select L1, resi 118-127
select L2, resi 164-176
select L3, resi 237-250

# N and C termini of DBD
select N_term, resi 94-111
select C_term, resi 287-292

# ================================================================
# COLOR BY REGION - Distinct colors
# ================================================================
# Beta sheets - Gold/Orange tones
color 0xF39C12, S1 or S2 or S3 or S4
color 0xE67E22, S5 or S6 or S7
color 0xD68910, S8 or S9 or S10

# Alpha helices - Teal/Blue
color 0x1ABC9C, H1
color 0x16A085, H2

# L1 loop - Light purple
color 0xA569BD, L1

# L2 loop (Zn binding, thermal sensitive) - Salmon
color 0xE74C3C, L2

# L3 loop (DNA binding, most thermal sensitive) - Deep red
color 0xC0392B, L3

# Termini - Gray
color 0x95A5A6, N_term
color 0x7F8C8D, C_term

# ================================================================
# ZINC ION AND COORDINATION
# ================================================================
select zinc_ion, name ZN
show spheres, zinc_ion
color 0x5D6D7E, zinc_ion
set sphere_scale, 1.2, zinc_ion

# Zinc coordinating residues
select zinc_coord, resi 176+179+238+242
show sticks, zinc_coord and (name SG or name NE2 or name ND1 or name CA or name CB or name CG)
color 0xF4D03F, zinc_coord

# ================================================================
# CANCER HOTSPOTS
# ================================================================
select mut_R175, resi 175 and name CA
select mut_G245, resi 245 and name CA
select mut_R248, resi 248 and name CA
select mut_R273, resi 273 and name CA
select mut_R282, resi 282 and name CA

# Show as spheres
show spheres, mut_R175 or mut_G245 or mut_R248 or mut_R273 or mut_R282

# Structural mutations - Magenta
color 0xFF006E, mut_R175
color 0xFF006E, mut_G245
color 0xFF006E, mut_R282

# Contact mutations - Cyan
color 0x00D4FF, mut_R248
color 0x00D4FF, mut_R273

set sphere_scale, 1.4, mut_R175 or mut_G245 or mut_R248 or mut_R273 or mut_R282

# ================================================================
# SET OPTIMAL VIEW
# ================================================================
reset
orient p53
turn y, 140
turn x, -25
turn z, 0
zoom p53, -8

# ================================================================
# LABELS - Positioned carefully to avoid overlap
# ================================================================
set label_font_id, 7
set label_size, 20
set label_color, black
set label_bg_color, white
set label_bg_transparency, 0.2

# --- STRUCTURAL REGION LABELS ---

# N-terminus
pseudoatom lab_Nterm, pos=[-30, 20, 30]
label lab_Nterm, "N-term (94)"

# C-terminus  
pseudoatom lab_Cterm, pos=[25, 30, 15]
label lab_Cterm, "C-term (292)"

# Beta sandwich label
pseudoatom lab_beta, pos=[0, -35, 0]
label lab_beta, "Beta-Sandwich Core"

# H1 helix
pseudoatom lab_H1, pos=[-35, -10, -15]
label lab_H1, "H1"

# H2 helix (DNA binding)
pseudoatom lab_H2, pos=[35, 15, 5]
label lab_H2, "H2 (DNA binding)"

# L1 loop
pseudoatom lab_L1, pos=[-20, 30, -5]
label lab_L1, "L1 Loop"

# L2 loop
pseudoatom lab_L2, pos=[-40, 0, 5]
label lab_L2, "L2 Loop"

# L3 loop
pseudoatom lab_L3, pos=[10, -20, -30]
label lab_L3, "L3 Loop (DNA contact)"

# Zinc site
pseudoatom lab_Zn, pos=[-25, -20, -25]
label lab_Zn, "Zn2+ site"

# --- HOTSPOT MUTATION LABELS ---
set label_size, 18

pseudoatom lab_R175, pos=[-15, 25, -5]
label lab_R175, "R175H*"

pseudoatom lab_G245, pos=[15, -10, -35]
label lab_G245, "G245S*"

pseudoatom lab_R248, pos=[30, -5, -20]
label lab_R248, "R248Q"

pseudoatom lab_R273, pos=[20, 25, -10]
label lab_R273, "R273H"

pseudoatom lab_R282, pos=[40, 5, 15]
label lab_R282, "R282W*"

# Hide pseudoatom spheres
hide spheres, lab_*

# ================================================================
# CONNECTOR LINES (optional - for clarity)
# ================================================================
# Draw thin lines from labels to structures
set dash_gap, 0.0
set dash_length, 0.5
set dash_width, 1.0
set dash_color, gray60

distance line_L2, lab_L2, resi 170 and name CA, mode=0
distance line_L3, lab_L3, resi 245 and name CA, mode=0
distance line_Zn, lab_Zn, zinc_ion, mode=0
distance line_H2, lab_H2, resi 280 and name CA, mode=0

distance line_m175, lab_R175, mut_R175, mode=0
distance line_m245, lab_G245, mut_G245, mode=0
distance line_m248, lab_R248, mut_R248, mode=0
distance line_m273, lab_R273, mut_R273, mode=0
distance line_m282, lab_R282, mut_R282, mode=0

hide labels, line_*

# ================================================================
# FINAL RENDER
# ================================================================
ray 4000, 4000
png /home/raza/p53_thermal_xai/results/figures/publication/Figure1_Annotated.png, dpi=300

quit
