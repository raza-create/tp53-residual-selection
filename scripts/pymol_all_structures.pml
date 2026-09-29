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

# Output directory
output = "/home/raza/p53_thermal_xai/results/figures/pymol"

# ================================================================
# PART 1: P53 REFERENCE STRUCTURE (2XWR)
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
png /home/raza/p53_thermal_xai/results/figures/pymol/Fig1A_p53_cartoon.png, dpi=300

# Figure 1B: Cancer hotspot mutations
color gray80, p53_ref
select hotspots, p53_ref and resi 175+245+248+249+273+282
color red, hotspots
show spheres, hotspots and name CA
set sphere_scale, 1.2, hotspots
ray 2400, 1800
png /home/raza/p53_thermal_xai/results/figures/pymol/Fig1B_p53_hotspots.png, dpi=300

# Figure 1C: DNA binding residues
hide spheres
color gray80, p53_ref
select dna_bind, p53_ref and resi 248+273+276+277+280+281+120
color blue, dna_bind
show sticks, dna_bind
ray 2400, 1800
png /home/raza/p53_thermal_xai/results/figures/pymol/Fig1C_p53_dna_binding.png, dpi=300

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
png /home/raza/p53_thermal_xai/results/figures/pymol/Fig1D_p53_zinc_site.png, dpi=300

# Figure 1E: Surface with hotspots highlighted
hide all
show surface, p53_ref
color white, p53_ref
set transparency, 0.6
show spheres, hotspots and name CA
color red, hotspots
set sphere_scale, 1.0, hotspots
ray 2400, 1800
png /home/raza/p53_thermal_xai/results/figures/pymol/Fig1E_p53_surface.png, dpi=300

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
png /home/raza/p53_thermal_xai/results/figures/pymol/Fig2A_p63_cartoon.png, dpi=300

# Figure 2B: p63 hotspot equivalent positions
color gray80, p63_ref
select p63_hotspots, p63_ref and resi 204+274+277+278+302+311
color magenta, p63_hotspots
show spheres, p63_hotspots and name CA
set sphere_scale, 1.2, p63_hotspots
ray 2400, 1800
png /home/raza/p53_thermal_xai/results/figures/pymol/Fig2B_p63_hotspots.png, dpi=300

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
png /home/raza/p53_thermal_xai/results/figures/pymol/Fig3A_p73_cartoon.png, dpi=300

# Figure 3B: p73 hotspot equivalent positions
color gray80, p73_ref
select p73_hotspots, p73_ref and resi 178+248+251+252+276+285
color cyan, p73_hotspots
show spheres, p73_hotspots and name CA
set sphere_scale, 1.2, p73_hotspots
ray 2400, 1800
png /home/raza/p53_thermal_xai/results/figures/pymol/Fig3B_p73_hotspots.png, dpi=300

# ================================================================
# PART 4: SUPERPOSITION OF ALL THREE PROTEINS
# ================================================================
delete all
fetch 2XWR, p53
fetch 3QYN, p63
fetch 2XWC, p73
remove solvent
remove not chain A

# Align structures
align p63, p53
align p73, p53

hide all
show cartoon, all
color red, p53
color green, p63
color blue, p73
orient all
zoom all, 5
ray 2400, 1800
png /home/raza/p53_thermal_xai/results/figures/pymol/Fig4_p53_p63_p73_superposition.png, dpi=300

# ================================================================
# PART 5: P53 COLORED BY THERMAL SENSITIVITY
# ================================================================
delete all
fetch 2XWR, p53_thermal
remove solvent
remove not chain A
hide all

# Color by regions that show temperature sensitivity
# Core (most affected) - red
# Intermediate - orange/yellow
# Stable - blue
show cartoon, p53_thermal
color white, p53_thermal

# L2 loop (highly flexible, temperature sensitive)
select L2_loop, p53_thermal and resi 164-194
color red, L2_loop

# L3 loop (DNA binding, temperature sensitive)
select L3_loop, p53_thermal and resi 237-250
color orange, L3_loop

# Beta sandwich core (more stable)
select beta_core, p53_thermal and resi 195-236
color yellow, beta_core

# Stable regions
select stable, p53_thermal and resi 94-163
color lightblue, stable

orient p53_thermal
zoom p53_thermal, 5
ray 2400, 1800
png /home/raza/p53_thermal_xai/results/figures/pymol/Fig5_p53_thermal_sensitivity.png, dpi=300

# ================================================================
# PART 6: NATIVE CONTACTS VISUALIZATION
# ================================================================
delete all
fetch 2XWR, p53_contacts
remove solvent
remove not chain A
hide all

show cartoon, p53_contacts
color gray80, p53_contacts

# Highlight residues involved in key native contacts
# These are core stabilizing interactions
select core_contacts, p53_contacts and resi 173+176+179+195+196+198+220+239+240+242+243+245+247+248
color green, core_contacts
show sticks, core_contacts and sidechain

orient p53_contacts
zoom p53_contacts, 5
ray 2400, 1800
png /home/raza/p53_thermal_xai/results/figures/pymol/Fig6_p53_native_contacts.png, dpi=300

quit
