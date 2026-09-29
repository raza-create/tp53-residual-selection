#!/usr/bin/env python3
"""
Professional Illustrative 3D p53 Figure
Similar to high-impact journal style
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle, FancyBboxPatch, Wedge, Arc, Polygon, Ellipse
from matplotlib.collections import PatchCollection
import matplotlib.patheffects as path_effects
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from pathlib import Path as FilePath

output_dir = FilePath('results/figures/publication')
output_dir.mkdir(parents=True, exist_ok=True)

print("="*70)
print("CREATING ILLUSTRATIVE 3D p53 FIGURE")
print("="*70)

# Create figure with light background
fig, ax = plt.subplots(figsize=(14, 16))
ax.set_xlim(0, 100)
ax.set_ylim(0, 120)
ax.axis('off')
ax.set_facecolor('#E8EDF2')
fig.patch.set_facecolor('#E8EDF2')

# ================================================================
# HELPER FUNCTIONS FOR 3D EFFECTS
# ================================================================

def draw_3d_sphere(ax, cx, cy, radius, color, label='', label_color='white', 
                   highlight=True, glow=False, glow_color=None):
    """Draw a 3D-looking sphere with highlight and optional glow"""
    from matplotlib.colors import to_rgb
    
    if glow:
        glow_c = glow_color if glow_color else color
        for i in range(5, 0, -1):
            glow_circle = Circle((cx, cy), radius + i*1.5, 
                                facecolor=glow_c, alpha=0.1, edgecolor='none')
            ax.add_patch(glow_circle)
    
    # Main sphere
    sphere = Circle((cx, cy), radius, facecolor=color, edgecolor='#333333', 
                   linewidth=2, zorder=10)
    ax.add_patch(sphere)
    
    # Highlight (top-left)
    if highlight:
        highlight_circle = Circle((cx - radius*0.3, cy + radius*0.3), radius*0.35,
                                  facecolor='white', alpha=0.5, edgecolor='none', zorder=11)
        ax.add_patch(highlight_circle)
    
    # Label
    if label:
        ax.text(cx, cy, label, ha='center', va='center', fontsize=14, 
               fontweight='bold', color=label_color, zorder=12)

def draw_gradient_bar(ax, x1, y1, x2, y2, width, colors, zorder=5):
    """Draw a gradient-colored bar between two points"""
    from matplotlib.colors import to_rgb
    
    # Calculate angle
    angle = np.arctan2(y2 - y1, x2 - x1)
    length = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    
    # Create segments
    n_segments = 20
    for i in range(n_segments):
        t1 = i / n_segments
        t2 = (i + 1) / n_segments
        
        # Interpolate position
        px1 = x1 + t1 * (x2 - x1)
        py1 = y1 + t1 * (y2 - y1)
        px2 = x1 + t2 * (x2 - x1)
        py2 = y1 + t2 * (y2 - y1)
        
        # Interpolate color
        if len(colors) == 2:
            r = colors[0][0] + t1 * (colors[1][0] - colors[0][0])
            g = colors[0][1] + t1 * (colors[1][1] - colors[0][1])
            b = colors[0][2] + t1 * (colors[1][2] - colors[0][2])
        else:
            r, g, b = colors[0]
        
        # Draw segment as thick line
        ax.plot([px1, px2], [py1, py2], color=(r, g, b), linewidth=width, 
               solid_capstyle='round', zorder=zorder)

def draw_helix(ax, x_start, y_start, x_end, y_end, n_turns=4, width=3, color='#9B59B6'):
    """Draw a 3D-looking helix"""
    n_points = n_turns * 20
    t = np.linspace(0, n_turns * 2 * np.pi, n_points)
    
    # Direction vector
    dx = x_end - x_start
    dy = y_end - y_start
    length = np.sqrt(dx**2 + dy**2)
    
    # Perpendicular vector
    px, py = -dy/length, dx/length
    
    for i in range(len(t) - 1):
        # Position along helix
        frac = i / (len(t) - 1)
        x = x_start + frac * dx
        y = y_start + frac * dy
        
        # Helix displacement
        offset = np.sin(t[i]) * width
        x += px * offset
        y += py * offset
        
        # Z-depth for color variation
        z = np.cos(t[i])
        
        # Color based on depth
        if z > 0:
            c = color
            alpha = 0.9
            lw = 6
        else:
            from matplotlib.colors import to_rgb
            rgb = to_rgb(color)
            c = tuple(max(0, v - 0.2) for v in rgb)
            alpha = 0.7
            lw = 4
        
        # Draw ellipse for each turn segment
        if i % 5 == 0:
            ellipse = Ellipse((x, y), 3, 5, angle=np.degrees(np.arctan2(dy, dx)),
                            facecolor=c, edgecolor='#333333', linewidth=0.5, 
                            alpha=alpha, zorder=5)
            ax.add_patch(ellipse)

def draw_dna_helix(ax, x_center, y_center, length=25, width=15):
    """Draw a double helix DNA"""
    n_points = 50
    t = np.linspace(0, 4 * np.pi, n_points)
    
    # Strand 1
    x1 = x_center + width/2 * np.sin(t)
    y1 = y_center - length/2 + (length/n_points) * np.arange(n_points)
    
    # Strand 2 (180 degrees out of phase)
    x2 = x_center + width/2 * np.sin(t + np.pi)
    y2 = y1
    
    # Draw backbone strands
    ax.plot(x1, y1, color='#3498DB', linewidth=4, zorder=3)
    ax.plot(x2, y2, color='#3498DB', linewidth=4, zorder=3)
    
    # Draw base pairs
    for i in range(0, n_points, 3):
        # Only draw when strands are at similar depth
        if abs(np.sin(t[i])) < 0.7:
            # Alternate colors for base pairs
            bp_color = '#2ECC71' if i % 6 == 0 else '#27AE60'
            ax.plot([x1[i], x2[i]], [y1[i], y2[i]], color=bp_color, linewidth=2, zorder=2)
            
            # Nucleotide circles
            circle1 = Circle((x1[i], y1[i]), 1, facecolor='#3498DB', edgecolor='#2980B9', 
                            linewidth=1, zorder=4)
            circle2 = Circle((x2[i], y2[i]), 1, facecolor='#2ECC71', edgecolor='#27AE60',
                            linewidth=1, zorder=4)
            ax.add_patch(circle1)
            ax.add_patch(circle2)

# ================================================================
# MAIN FIGURE ELEMENTS
# ================================================================

# --- N-terminus (green sphere) ---
draw_3d_sphere(ax, 15, 105, 5, '#2ECC71', 'N', 'white', glow=True, glow_color='#2ECC71')
ax.text(22, 107, '(1-40)', fontsize=10, color='#5D6D7E')

# --- TAD1 label ---
bbox_tad1 = dict(boxstyle='round,pad=0.3', facecolor='#E74C3C', edgecolor='none')
ax.text(18, 100, 'TAD1', fontsize=11, fontweight='bold', color='white', bbox=bbox_tad1)

# --- Gradient bar from N-term to TAD2 ---
draw_gradient_bar(ax, 18, 102, 30, 85, 8, [(0.91, 0.3, 0.24), (0.95, 0.5, 0.3)])

# --- TAD2 label ---
ax.text(33, 83, '(41-61)', fontsize=10, color='#5D6D7E')
bbox_tad2 = dict(boxstyle='round,pad=0.3', facecolor='#F39C12', edgecolor='none')
ax.text(28, 78, 'TAD2', fontsize=11, fontweight='bold', color='white', bbox=bbox_tad2)

# --- PRD (Proline-rich domain) ---
ax.text(38, 72, '(62-93)', fontsize=10, color='#5D6D7E')
bbox_prd = dict(boxstyle='round,pad=0.3', facecolor='#F5B041', edgecolor='none')
ax.text(33, 68, 'PRD', fontsize=11, fontweight='bold', color='white', bbox=bbox_prd)

# --- Gradient bar from TAD2 to DBD ---
draw_gradient_bar(ax, 32, 80, 40, 62, 7, [(0.95, 0.5, 0.3), (0.96, 0.7, 0.3)])

# ================================================================
# DBD - Large 3D Sphere (Central element)
# ================================================================
dbd_x, dbd_y, dbd_r = 50, 55, 18

# Shadow
shadow = Circle((dbd_x + 2, dbd_y - 2), dbd_r, facecolor='#2C3E50', alpha=0.3, 
               edgecolor='none', zorder=1)
ax.add_patch(shadow)

# Main DBD sphere with gradient effect
dbd_main = Circle((dbd_x, dbd_y), dbd_r, facecolor='#3498DB', edgecolor='#2C3E50',
                 linewidth=3, zorder=5)
ax.add_patch(dbd_main)

# Highlight arc (top)
highlight_arc = Wedge((dbd_x - 5, dbd_y + 8), dbd_r * 0.6, 30, 150, 
                     facecolor='white', alpha=0.3, edgecolor='none', zorder=6)
ax.add_patch(highlight_arc)

# DBD label
ax.text(dbd_x, dbd_y + 2, 'DBD', fontsize=22, fontweight='bold', color='white',
       ha='center', va='center', zorder=7)
ax.text(dbd_x, dbd_y - 6, '(94-292)', fontsize=12, color='white',
       ha='center', va='center', zorder=7)

# --- Cancer hotspot mutations (red spheres on DBD) ---
hotspots = [
    (dbd_x - 12, dbd_y + 12, 'R175H'),
    (dbd_x + 12, dbd_y + 12, 'R248'),
    (dbd_x - 14, dbd_y - 2, 'G245'),
    (dbd_x + 14, dbd_y - 2, 'R273'),
    (dbd_x - 8, dbd_y - 12, 'R282'),
    (dbd_x + 8, dbd_y + 5, 'R249'),
]

for hx, hy, name in hotspots:
    # Red sphere with glow
    glow = Circle((hx, hy), 3.5, facecolor='#E74C3C', alpha=0.3, edgecolor='none', zorder=6)
    ax.add_patch(glow)
    hotspot_circle = Circle((hx, hy), 2.2, facecolor='#E74C3C', edgecolor='#922B21',
                           linewidth=1, zorder=8)
    ax.add_patch(hotspot_circle)

# --- Zinc binding sites (yellow spheres) ---
zn_sites = [
    (dbd_x - 6, dbd_y + 6),
    (dbd_x + 4, dbd_y - 8),
]

for zx, zy in zn_sites:
    # Yellow glow
    glow = Circle((zx, zy), 3, facecolor='#F4D03F', alpha=0.4, edgecolor='none', zorder=6)
    ax.add_patch(glow)
    zn_circle = Circle((zx, zy), 1.8, facecolor='#F4D03F', edgecolor='#B7950B',
                       linewidth=1, zorder=8)
    ax.add_patch(zn_circle)

# Zn2+ label
ax.text(dbd_x - 10, dbd_y + 10, 'Zn²⁺', fontsize=9, fontweight='bold', color='#B7950B', zorder=9)

# ================================================================
# OD - Oligomerization Domain (Helix)
# ================================================================
# Draw helix
helix_x_start, helix_y_start = 70, 58
helix_x_end, helix_y_end = 82, 45

for i in range(6):
    ellipse = Ellipse((helix_x_start + i*2.2, helix_y_start - i*2.5), 4, 7,
                     angle=-50, facecolor='#9B59B6', edgecolor='#6C3483',
                     linewidth=1, alpha=0.9, zorder=5)
    ax.add_patch(ellipse)
    # Shadow
    ellipse_shadow = Ellipse((helix_x_start + i*2.2 + 0.5, helix_y_start - i*2.5 - 0.5), 4, 7,
                            angle=-50, facecolor='#4A235A', alpha=0.3, edgecolor='none', zorder=4)
    ax.add_patch(ellipse_shadow)

# OD label
bbox_od = dict(boxstyle='round,pad=0.3', facecolor='#8E44AD', edgecolor='none')
ax.text(88, 40, 'OD', fontsize=11, fontweight='bold', color='white', bbox=bbox_od)
ax.text(88, 35, '(324-355)', fontsize=10, color='#5D6D7E')

# --- Linker to CTD ---
ax.plot([85, 90], [42, 50], color='#7F8C8D', linewidth=5, solid_capstyle='round', zorder=3)

# --- CTD label ---
bbox_ctd = dict(boxstyle='round,pad=0.3', facecolor='#7F8C8D', edgecolor='none')
ax.text(88, 58, 'CTD', fontsize=10, fontweight='bold', color='white', bbox=bbox_ctd)
ax.text(88, 63, '(356-393)', fontsize=10, color='#5D6D7E')

# --- C-terminus (red sphere) ---
draw_3d_sphere(ax, 92, 68, 4, '#E74C3C', 'C', 'white', glow=True, glow_color='#E74C3C')

# ================================================================
# DNA Double Helix (Bottom)
# ================================================================
ax.annotate('', xy=(50, 30), xytext=(50, 37),
           arrowprops=dict(arrowstyle='->', color='#2C3E50', lw=3))

draw_dna_helix(ax, 50, 15, length=20, width=12)

ax.text(50, 2, 'Target DNA', fontsize=14, fontweight='bold', color='#2C3E50',
       ha='center', va='center')

# ================================================================
# LEGEND BOXES
# ================================================================

# Cancer Hotspots legend
legend1_x, legend1_y = 5, 65
rect1 = FancyBboxPatch((legend1_x, legend1_y), 22, 25,
                       boxstyle="round,pad=0.02,rounding_size=1",
                       facecolor='white', edgecolor='#E74C3C', linewidth=2)
ax.add_patch(rect1)

ax.text(legend1_x + 11, legend1_y + 22, 'Cancer Hotspots', fontsize=11, 
       fontweight='bold', color='#E74C3C', ha='center')

hotspot_list = ['R175H', 'R248', 'R273', 'R282', 'G245', 'R249']
for i, hs in enumerate(hotspot_list):
    y_pos = legend1_y + 17 - i * 3
    circle = Circle((legend1_x + 4, y_pos), 1.2, facecolor='#E74C3C', edgecolor='#922B21')
    ax.add_patch(circle)
    ax.text(legend1_x + 7, y_pos, hs, fontsize=9, color='#2C3E50', va='center')

# Zinc Binding Sites legend
legend2_x, legend2_y = 5, 45
rect2 = FancyBboxPatch((legend2_x, legend2_y), 22, 12,
                       boxstyle="round,pad=0.02,rounding_size=1",
                       facecolor='white', edgecolor='#F4D03F', linewidth=2)
ax.add_patch(rect2)

ax.text(legend2_x + 11, legend2_y + 9, 'Zinc Binding Sites', fontsize=11,
       fontweight='bold', color='#B7950B', ha='center')
circle_zn = Circle((legend2_x + 4, legend2_y + 4), 1.5, facecolor='#F4D03F', 
                   edgecolor='#B7950B')
ax.add_patch(circle_zn)
ax.text(legend2_x + 7, legend2_y + 4, 'Structural Zn²⁺', fontsize=9, color='#2C3E50', va='center')

# ================================================================
# TITLE
# ================================================================
ax.text(50, 115, 'Human p53 Tumor Suppressor Protein', fontsize=18, 
       fontweight='bold', color='#2C3E50', ha='center')
ax.text(50, 111, 'Domain Architecture (393 amino acids)', fontsize=12,
       color='#5D6D7E', ha='center', style='italic')

# ================================================================
# FRAME
# ================================================================
frame = FancyBboxPatch((1, 1), 98, 117,
                       boxstyle="round,pad=0.02,rounding_size=2",
                       facecolor='none', edgecolor='#3498DB', linewidth=4)
ax.add_patch(frame)

# ================================================================
# SAVE
# ================================================================
plt.savefig(output_dir / 'Figure1_p53_Illustrative_3D.png', dpi=300, bbox_inches='tight',
            facecolor='#E8EDF2', edgecolor='none')
plt.savefig(output_dir / 'Figure1_p53_Illustrative_3D.pdf', bbox_inches='tight',
            facecolor='#E8EDF2', edgecolor='none')
plt.savefig(output_dir / 'Figure1_p53_Illustrative_3D.tiff', dpi=300, bbox_inches='tight',
            facecolor='#E8EDF2', edgecolor='none')
plt.close()

print("\n" + "="*70)
print("ILLUSTRATIVE 3D FIGURE CREATED SUCCESSFULLY")
print("="*70)
print(f"\nSaved:")
print(f"  - {output_dir}/Figure1_p53_Illustrative_3D.png")
print(f"  - {output_dir}/Figure1_p53_Illustrative_3D.pdf")
print(f"  - {output_dir}/Figure1_p53_Illustrative_3D.tiff")
print("\nFigure includes:")
print("  ✓ N-terminus (green 3D sphere)")
print("  ✓ TAD1, TAD2 (gradient bars)")
print("  ✓ PRD (proline-rich domain)")
print("  ✓ DBD (large blue 3D sphere with hotspots & Zn)")
print("  ✓ OD (purple helix)")
print("  ✓ CTD and C-terminus")
print("  ✓ Target DNA double helix")
print("  ✓ Legend boxes")
print("  ✓ Professional styling")
