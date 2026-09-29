# ================================================================
# HIGH-IMPACT JOURNAL INTRODUCTION FIGURE
# p53 DNA-Binding Domain - 3D Structural Overview
# ================================================================

# Fetch structure
fetch 2XWR, p53
remove solvent
remove not chain A

# High-quality settings
bg_color white
set ray_shadows, 1
set ray_shadow_decay_factor, 0.1
set ray_shadow_decay_range, 2
set antialias, 2
set ray_trace_mode, 1
set ray_trace_gain, 0.1
set ray_trace_color, black
set ambient, 0.4
set spec_reflect, 0.5
set shininess, 10
set specular, 1
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_cylindrical_helices, 1
set cartoon_tube_radius, 0.4
set stick_radius, 0.25
set sphere_scale, 0.3
set depth_cue, 1
set fog_start, 0.45

# ================================================================
# FIGURE 1A: Overall Structure - Secondary Structure Colored
# ================================================================
hide all
show cartoon, p53

# Color by secondary structure
color tv_blue, ss h
color tv_orange, ss s
color white, ss l+''

# Show zinc as sphere
select zinc, name ZN
show spheres, zinc
color gray50, zinc
set sphere_scale, 0.8, zinc

# Optimal orientation
orient p53
turn y, 20
turn x, -10
zoom p53, -5

# Render
ray 2400, 2400
png /home/raza/p53_thermal_xai/results/figures/publication/Fig1_p53_structure_3D.png, dpi=300

# ================================================================
# FIGURE 1B: Cancer Hotspots Highlighted
# ================================================================
# Reset coloring
color gray80, p53

# Define hotspot residues
select hotspot_R175, resi 175
select hotspot_G245, resi 245
select hotspot_R248, resi 248
select hotspot_R249, resi 249
select hotspot_R273, resi 273
select hotspot_R282, resi 282
select all_hotspots, hotspot_R175 or hotspot_G245 or hotspot_R248 or hotspot_R249 or hotspot_R273 or hotspot_R282

# Color hotspots - Structural (red) vs Contact (blue)
color firebrick, hotspot_R175  # Structural
color firebrick, hotspot_G245  # Structural
color marine, hotspot_R248     # Contact
color firebrick, hotspot_R249  # Structural
color marine, hotspot_R273     # Contact
color firebrick, hotspot_R282  # Structural

# Show hotspots as spheres
show spheres, all_hotspots and name CA
set sphere_scale, 1.8, all_hotspots

# Show sidechains
show sticks, all_hotspots and not name C+N+O

# Add labels
set label_color, black
set label_size, 22
set label_font_id, 7
set label_shadow_mode, 2

label hotspot_R175 and name CA, "R175"
label hotspot_G245 and name CA, "G245"
label hotspot_R248 and name CA, "R248"
label hotspot_R249 and name CA, "R249"
label hotspot_R273 and name CA, "R273"
label hotspot_R282 and name CA, "R282"

# Optimal view
orient p53
turn y, 30
turn x, -15
zoom p53, -3

ray 2400, 2400
png /home/raza/p53_thermal_xai/results/figures/publication/Fig1_p53_hotspots_3D.png, dpi=300

# ================================================================
# FIGURE 1C: DNA Binding Surface
# ================================================================
hide labels
hide spheres
hide sticks
color gray80, p53

# L2 loop (164-194) - important for DNA binding
select L2_loop, resi 164-194
color salmon, L2_loop

# L3 loop (237-250) - critical for DNA contact
select L3_loop, resi 237-250
color tv_red, L3_loop

# Helix H2 (278-286) - DNA binding helix
select H2_helix, resi 278-286
color deepblue, H2_helix

# DNA contact residues
select dna_contact, resi 248+273+276+277+280+281
show sticks, dna_contact
color marine, dna_contact

# Zinc site
select zinc_site, resi 176+179+238+242
show sticks, zinc_site
color orange, zinc_site
show spheres, zinc
color gray60, zinc

# View focusing on DNA binding surface
orient p53
turn y, -60
turn x, 10
zoom p53, -5

ray 2400, 2400
png /home/raza/p53_thermal_xai/results/figures/publication/Fig1_p53_dna_surface_3D.png, dpi=300

# ================================================================
# FIGURE 1D: Thermal Sensitivity Regions
# ================================================================
hide sticks
hide spheres
color white, p53

# Color by thermal sensitivity (based on SHAP analysis)
# Most stable (blue) to most sensitive (red)

# Stable core - beta sandwich
select beta_core, resi 115-135 or resi 195-210 or resi 225-235 or resi 251-265
color lightblue, beta_core

# Moderately stable
select moderate, resi 94-114 or resi 268-292
color palegreen, moderate

# Thermally sensitive - loops
select sensitive_L1, resi 136-163
color yellow, sensitive_L1

# Highly sensitive - L2
select sensitive_L2, resi 164-194
color orange, sensitive_L2

# Most sensitive - L3 (DNA binding)
select sensitive_L3, resi 237-250
color firebrick, sensitive_L3

# Zinc sphere
show spheres, zinc
color gray50, zinc
set sphere_scale, 0.8, zinc

orient p53
turn y, 20
turn x, -10
zoom p53, -5

ray 2400, 2400
png /home/raza/p53_thermal_xai/results/figures/publication/Fig1_p53_thermal_sensitivity_3D.png, dpi=300

# ================================================================
# FIGURE 1E: Surface View with Hotspots
# ================================================================
hide all
show surface, p53
set transparency, 0.75, p53
color white, p53

# Show cartoon inside
show cartoon, p53
color gray60, p53
set cartoon_transparency, 0.0

# Hotspots as colored spheres on surface
show spheres, all_hotspots and name CA
set sphere_scale, 1.5, all_hotspots
color firebrick, hotspot_R175
color firebrick, hotspot_G245
color marine, hotspot_R248
color firebrick, hotspot_R249
color marine, hotspot_R273
color firebrick, hotspot_R282

orient p53
turn y, 25
turn x, -15
zoom p53, -3

ray 2400, 2400
png /home/raza/p53_thermal_xai/results/figures/publication/Fig1_p53_surface_hotspots_3D.png, dpi=300

# ================================================================
# FIGURE 1F: Complete Overview (Final Intro Figure)
# ================================================================
hide all
show cartoon, p53

# Rainbow coloring from N to C terminus
spectrum count, rainbow, p53

# Show key structural elements
# Zinc
show spheres, zinc
color gray50, zinc
set sphere_scale, 1.0, zinc

# Hotspots as red spheres
show spheres, all_hotspots and name CA
color red, all_hotspots
set sphere_scale, 1.2, all_hotspots

# Nice orientation
reset
orient p53
turn y, 30
turn x, -20
zoom p53, -5

ray 3000, 3000
png /home/raza/p53_thermal_xai/results/figures/publication/Fig1_p53_intro_final_3D.png, dpi=300

# ================================================================
# FIGURE 1G: Multiple Views Panel
# ================================================================
# Front view
reset
orient p53
zoom p53, -5
ray 1800, 1800
png /home/raza/p53_thermal_xai/results/figures/publication/Fig1_view_front.png, dpi=300

# Side view (90 degree rotation)
turn y, 90
ray 1800, 1800
png /home/raza/p53_thermal_xai/results/figures/publication/Fig1_view_side.png, dpi=300

# Top view
reset
orient p53
turn x, 90
zoom p53, -5
ray 1800, 1800
png /home/raza/p53_thermal_xai/results/figures/publication/Fig1_view_top.png, dpi=300

quit
