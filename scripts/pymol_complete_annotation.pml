# ================================================================
# P53 DBD COMPLETE STRUCTURAL ANNOTATION
# ALL residues (94-292) with exact ranges
# For high-impact journal
# ================================================================

fetch 2XWR, p53
remove solvent
remove not chain A

# ================================================================
# RENDERING SETTINGS
# ================================================================
bg_color white
set ray_shadows, 1
set antialias, 3
set ray_trace_mode, 1
set ambient, 0.45
set spec_reflect, 0.5
set shininess, 25
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_cylindrical_helices, 1
set cartoon_tube_radius, 0.4
set cartoon_loop_radius, 0.3
set depth_cue, 1
set fog_start, 0.6

hide all
show cartoon, p53

# ================================================================
# COMPLETE P53 DBD SECONDARY STRUCTURE (94-292)
# Total: 199 residues
# ================================================================

# --- N-TERMINAL REGION (94-102) ---
select Nterm, resi 94-102
color 0x95A5A6, Nterm

# --- STRAND S1 (103-109) ---
select S1, resi 103-109
color 0xF5B041, S1

# --- LOOP L0 (110-111) ---
select L0, resi 110-111
color 0xD5D8DC, L0

# --- STRAND S2 (112-117) ---
select S2, resi 112-117
color 0xF5B041, S2

# --- LOOP BETWEEN S2-S3 (118-123) ---
select loop_S2S3, resi 118-123
color 0xD5D8DC, loop_S2S3

# --- STRAND S3 (124-127) ---
select S3, resi 124-127
color 0xF5B041, S3

# --- STRAND S4 (128-132) ---
select S4, resi 128-132
color 0xF5B041, S4

# --- LOOP L1 (133-141) ---
select L1_loop, resi 133-141
color 0xA569BD, L1_loop

# --- STRAND S5 (142-146) ---
select S5, resi 142-146
color 0xF5B041, S5

# --- LOOP BETWEEN S5-S6 (147-150) ---
select loop_S5S6, resi 147-150
color 0xD5D8DC, loop_S5S6

# --- STRAND S6 (151-155) ---
select S6, resi 151-155
color 0xF5B041, S6

# --- LOOP BETWEEN S6-S7 (156-162) ---
select loop_S6S7, resi 156-162
color 0xD5D8DC, loop_S6S7

# --- STRAND S7 (163-166) ---
select S7, resi 163-166
color 0xF5B041, S7

# --- L2 LOOP - ZINC BINDING REGION (167-195) ---
# This is the large loop containing zinc coordination
# Subdivided:
# L2a (167-176): Before zinc
select L2a, resi 167-176
color 0xE74C3C, L2a

# HELIX H1 within L2 (177-181)
select H1, resi 177-181
color 0x1ABC9C, H1

# L2b (182-195): After H1, containing C176, H179 zinc ligands
select L2b, resi 182-195
color 0xE74C3C, L2b

# Full L2 selection
select L2_full, resi 167-195

# --- STRAND S8 (196-198) ---
select S8, resi 196-198
color 0xF5B041, S8

# --- STRAND S9 (199-205) ---
select S9, resi 199-205
color 0xF5B041, S9

# --- LOOP BETWEEN S9-S10 (206-212) ---
select loop_S9S10, resi 206-212
color 0xD5D8DC, loop_S9S10

# --- STRAND S10 (213-220) ---
select S10, resi 213-220
color 0xF5B041, S10

# --- LOOP BETWEEN S10-S11 (221-228) ---
select loop_S10S11, resi 221-228
color 0xD5D8DC, loop_S10S11

# --- STRAND S11 (229-232) ---
select S11, resi 229-232
color 0xF5B041, S11

# --- LOOP BETWEEN S11-S12 (233-236) ---
select loop_S11S12, resi 233-236
color 0xD5D8DC, loop_S11S12

# --- L3 LOOP - DNA BINDING (237-250) ---
# Most critical for DNA contact and thermal sensitivity
select L3_loop, resi 237-250
color 0xC0392B, L3_loop

# --- STRAND S12 (251-254) ---
select S12, resi 251-254
color 0xF5B041, S12

# --- LOOP BETWEEN S12-S13 (255-260) ---
select loop_S12S13, resi 255-260
color 0xD5D8DC, loop_S12S13

# --- STRAND S13 (261-265) ---
select S13, resi 261-265
color 0xF5B041, S13

# --- LOOP BEFORE H2 (266-277) ---
select loop_preH2, resi 266-277
color 0xD5D8DC, loop_preH2

# --- HELIX H2 - DNA BINDING HELIX (278-286) ---
select H2, resi 278-286
color 0x16A085, H2

# --- C-TERMINAL (287-292) ---
select Cterm, resi 287-292
color 0x7F8C8D, Cterm

# ================================================================
# ZINC ION AND COORDINATION
# ================================================================
select zinc_ion, name ZN
show spheres, zinc_ion
color 0x5D6D7E, zinc_ion
set sphere_scale, 1.3, zinc_ion

# Zinc coordinating residues: C176, H179, C238, C242
select zinc_C176, resi 176
select zinc_H179, resi 179
select zinc_C238, resi 238
select zinc_C242, resi 242
select zinc_ligands, zinc_C176 or zinc_H179 or zinc_C238 or zinc_C242

show sticks, zinc_ligands and (name SG or name NE2 or name ND1 or name CA or name CB or name CG or name CD2 or name CE1)
color 0xF4D03F, zinc_ligands

# ================================================================
# CANCER HOTSPOT MUTATIONS
# ================================================================
# R175 - in L2 loop (structural mutation)
select mut_R175, resi 175 and name CA
show spheres, mut_R175
color 0xFF006E, mut_R175
set sphere_scale, 1.5, mut_R175

# G245 - in L3 loop (structural mutation)
select mut_G245, resi 245 and name CA
show spheres, mut_G245
color 0xFF006E, mut_G245
set sphere_scale, 1.5, mut_G245

# R248 - in L3 loop (DNA contact mutation)
select mut_R248, resi 248 and name CA
show spheres, mut_R248
color 0x00D4FF, mut_R248
set sphere_scale, 1.5, mut_R248

# R249 - in L3 loop (structural mutation)
select mut_R249, resi 249 and name CA
show spheres, mut_R249
color 0xFF006E, mut_R249
set sphere_scale, 1.5, mut_R249

# R273 - in loop before H2 (DNA contact mutation)
select mut_R273, resi 273 and name CA
show spheres, mut_R273
color 0x00D4FF, mut_R273
set sphere_scale, 1.5, mut_R273

# R282 - in H2 helix (structural mutation)
select mut_R282, resi 282 and name CA
show spheres, mut_R282
color 0xFF006E, mut_R282
set sphere_scale, 1.5, mut_R282

# ================================================================
# SET OPTIMAL VIEWING ANGLE
# ================================================================
reset
orient p53
turn y, 145
turn x, -20
zoom p53, -10

# ================================================================
# LABELS - All structural elements
# ================================================================
set label_font_id, 7
set label_size, 16
set label_color, black
set label_bg_color, white
set label_bg_transparency, 0.15

# Hide all pseudoatom spheres after creation
# N and C termini
pseudoatom pN, pos=[-35, 25, 25]
label pN, "N-term (94-102)"
pseudoatom pC, pos=[30, 30, 10]
label pC, "C-term (287-292)"

# Beta strands (grouped)
pseudoatom pBeta, pos=[5, -40, 0]
label pBeta, "Beta-sandwich (S1-S13)"

# L1 loop
pseudoatom pL1, pos=[-25, 25, -15]
label pL1, "L1 (133-141)"

# L2 loop with details
pseudoatom pL2, pos=[-45, -5, 0]
label pL2, "L2 (167-195)"

# H1 helix
pseudoatom pH1, pos=[-40, -15, -20]
label pH1, "H1 (177-181)"

# L3 loop
pseudoatom pL3, pos=[15, -25, -30]
label pL3, "L3 (237-250)"

# H2 helix
pseudoatom pH2, pos=[40, 15, 0]
label pH2, "H2 (278-286)"

# Zinc site
pseudoatom pZn, pos=[-30, -25, -30]
label pZn, "Zn2+ (C176,H179,C238,C242)"

# Cancer hotspots with positions
set label_size, 14
pseudoatom pm175, pos=[-20, 20, 0]
label pm175, "R175H (4.8%)"

pseudoatom pm245, pos=[10, -15, -40]
label pm245, "G245S (2.7%)"

pseudoatom pm248, pos=[30, -10, -25]
label pm248, "R248Q/W (6.4%)"

pseudoatom pm249, pos=[25, -20, -35]
label pm249, "R249S (2.5%)"

pseudoatom pm273, pos=[25, 22, -15]
label pm273, "R273H/C (4.5%)"

pseudoatom pm282, pos=[45, 5, 10]
label pm282, "R282W (1.8%)"

# Hide pseudoatom spheres
hide spheres, p*

# ================================================================
# CONNECTOR LINES
# ================================================================
set dash_gap, 0.0
set dash_length, 0.4
set dash_width, 1.0
set dash_color, gray60

distance ln_L2, pL2, resi 180 and name CA, mode=0
distance ln_L3, pL3, resi 243 and name CA, mode=0
distance ln_H1, pH1, resi 179 and name CA, mode=0
distance ln_H2, pH2, resi 282 and name CA, mode=0
distance ln_Zn, pZn, zinc_ion, mode=0

distance ln_175, pm175, mut_R175, mode=0
distance ln_245, pm245, mut_G245, mode=0
distance ln_248, pm248, mut_R248, mode=0
distance ln_249, pm249, mut_R249, mode=0
distance ln_273, pm273, mut_R273, mode=0
distance ln_282, pm282, mut_R282, mode=0

hide labels, ln_*

# ================================================================
# LEGEND TEXT (bottom of figure)
# ================================================================
pseudoatom leg1, pos=[-50, -45, 0]
label leg1, "Colors: Beta-strands (gold), H1/H2 helices (teal), L2 loop (red), L3 loop (dark red), Other loops (gray)"

pseudoatom leg2, pos=[-50, -50, 0]
label leg2, "Hotspots: Structural mutations (magenta), DNA contact mutations (cyan)"

hide spheres, leg*
set label_size, 12, leg*

# ================================================================
# RENDER
# ================================================================
ray 4500, 4500
png /home/raza/p53_thermal_xai/results/figures/publication/Figure1_Complete_Annotated.png, dpi=300

quit
