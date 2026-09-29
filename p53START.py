import matplotlib.pyplot as plt
import numpy as np

# Data derived from analysis
categories = ['248K', '273K', '310K', '350K', '400K', '450K']
true_temps = np.array([248, 273, 310, 350, 400, 450])
pred_temps = np.array([290, 331, 349, 350, 400, 399])
mae_labels = ["MAE=42K", "MAE=58K", "MAE=39K", "MAE=37K", "MAE=40K", "MAE=51K"]

x = np.arange(len(categories))
width = 0.35  # Width of the bars

# Setup figure with a clean, professional style
fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

# Colors: Professional Blue (True) and Red (Predicted)
color_true = '#377eb8' 
color_pred = '#e41a1c' 

# Plotting bars
rects1 = ax.bar(x - width/2, true_temps, width, label='True Temperature', 
                color=color_true, edgecolor='black', linewidth=0.8, alpha=0.9)
rects2 = ax.bar(x + width/2, pred_temps, width, label='Predicted Temperature', 
                color=color_pred, edgecolor='black', linewidth=0.8, alpha=0.9)

# Formatting Axes
ax.set_ylabel('Temperature (K)', fontsize=12, fontweight='bold')
ax.set_xlabel('Test Temperature', fontsize=12, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=11)
ax.set_ylim(0, 550)

# Legend
ax.legend(loc='upper left', frameon=True, fontsize=10, framealpha=0.95)

# Grid and Spines
ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=0.3)
ax.set_axisbelow(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Add MAE annotations above each bar group
for i in range(len(categories)):
    group_x = x[i]
    max_h = max(true_temps[i], pred_temps[i])
    ax.text(group_x, max_h + 10, mae_labels[i], 
            ha='center', va='bottom', fontsize=10, fontweight='bold', color='#333333')

plt.tight_layout()
plt.savefig('reconstructed_plot.png', dpi=300)
plt.show()
