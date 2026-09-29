import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# Use 'Agg' backend for Linux environments without a display server
import matplotlib
matplotlib.use('Agg')

def generate_p53_figure():
    fig, ax = plt.subplots(figsize=(14, 7), dpi=300)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6)
    ax.axis('off')

    # Background styling
    fig.patch.set_facecolor('white')

    # Titles
    plt.text(6, 5.6, 'Human p53 Tumor Suppressor Protein Structure', 
             fontsize=18, fontweight='bold', ha='center', family='sans-serif')
    plt.text(6, 5.2, 'Total Length: 393 Amino Acids', fontsize=12, fontstyle='italic', ha='center')

    # --- N-TERMINUS & TADs ---
    # N-Term Circle
    ax.add_patch(patches.Circle((1, 3), 0.45, color='#77DD77', ec='black', lw=1.5, zorder=2))
    plt.text(1, 3, 'N', fontsize=14, fontweight='bold', ha='center', va='center')
    plt.text(1, 2.3, '1-40', fontsize=10, fontweight='bold', ha='center')

    # TAD1 & TAD2 Rectangles
    ax.add_patch(patches.Rectangle((1.8, 2.7), 0.9, 0.6, color='#FF6961', ec='black', lw=1.5, zorder=2))
    ax.add_patch(patches.Rectangle((2.7, 2.7), 0.9, 0.6, color='#FFB347', ec='black', lw=1.5, zorder=2))
    plt.text(2.25, 3, 'TAD1', fontsize=10, fontweight='bold', ha='center', va='center')
    plt.text(3.15, 3, 'TAD2', fontsize=10, fontweight='bold', ha='center', va='center')
    plt.text(2.7, 2.3, '41-61', fontsize=10, fontweight='bold', ha='center')

    # --- DNA BINDING DOMAIN (DBD) ---
    ax.add_patch(patches.Circle((5.8, 3), 1.3, color='#89CFF0', ec='black', lw=2, zorder=1))
    plt.text(5.8, 4.5, 'DNA Binding Domain (DBD)', fontsize=13, fontweight='bold', ha='center')
    plt.text(5.8, 1.5, 'Residues: 94-292', fontsize=10, fontweight='bold', ha='center')

    # Zinc Coordination
    for x_offset in [-0.3, 0.3]:
        ax.add_patch(patches.Circle((5.8 + x_offset, 3.2), 0.18, color='#FDFD96', ec='black', lw=1, zorder=3))
    plt.text(5.8, 3.2, 'Zn$^{2+}$', fontsize=10, fontweight='bold', ha='center', va='center')

    # Hotspot Mutations
    muts = [(4.8, 3.6, 'R175H'), (4.7, 2.8, 'G245'), (5.3, 2.1, 'R248'), 
            (6.3, 2.1, 'R249'), (6.9, 2.8, 'R273'), (6.8, 3.6, 'R282')]
    for x, y, txt in muts:
        ax.add_patch(patches.Circle((x, y), 0.09, color='red', ec='black', zorder=4))
        plt.text(x, y+0.18, txt, fontsize=9, fontweight='bold', ha='center', color='darkred')

    # --- OLIGOMERIZATION DOMAIN (OD) ---
    for i in range(6):
        ax.add_patch(patches.Ellipse((8.2 + i*0.12, 3), 0.25, 0.8, angle=15, color='#B19CD9', ec='black', alpha=0.8))
    plt.text(8.5, 4.0, 'OD', fontsize=12, fontweight='bold', ha='center')
    plt.text(8.5, 2.3, '324-355', fontsize=10, fontweight='bold', ha='center')

    # --- C-TERMINAL DOMAIN (CTD) ---
    ax.add_patch(patches.Rectangle((9.5, 2.7), 1.8, 0.6, color='#CFCFC4', ec='black', lw=1.5))
    plt.text(10.4, 3, 'CTD', fontsize=12, fontweight='bold', ha='center', va='center')
    plt.text(10.4, 2.3, '356-393', fontsize=10, fontweight='bold', ha='center')

    # --- TARGET DNA ---
    x_dna = np.linspace(4.5, 7.1, 100)
    y_dna1 = 0.6 + 0.15 * np.sin(2 * np.pi * (x_dna - 4.5) / 1.0)
    y_dna2 = 0.3 + 0.15 * np.sin(2 * np.pi * (x_dna - 4.5) / 1.0)
    plt.plot(x_dna, y_dna1, color='#008080', lw=2.5, zorder=5)
    plt.plot(x_dna, y_dna2, color='#008080', lw=2.5, zorder=5)
    plt.text(5.8, -0.1, 'Target DNA Interaction', color='#008080', fontweight='bold', ha='center')
    ax.annotate('', xy=(5.8, 0.8), xytext=(5.8, 1.7), arrowprops=dict(arrowstyle='->', lw=2))

    # Save the figure
    output_file = "p53_structure_publication.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    print(f"Successfully generated: {output_file}")

if __name__ == "__main__":
    generate_p53_figure()
