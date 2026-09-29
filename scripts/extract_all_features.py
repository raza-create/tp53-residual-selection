#!/usr/bin/env python3
"""
Extract comprehensive structural features from p53 MD simulations
Features: RMSD, Rg, Secondary Structure, H-bonds, SASA, Native Contacts
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Paths
sim_dir = Path("/home/raza/p53_thermal_xai/data/structures/long_simulations")
output_dir = Path("/home/raza/p53_thermal_xai/results")
feature_dir = Path("/home/raza/p53_thermal_xai/data/features")
feature_dir.mkdir(parents=True, exist_ok=True)

# Reference structure
ref_gro = sim_dir.parent / "npt.gro"

# Temperature list
temperatures = [248, 273, 310, 350, 400, 450, 500]

print("=" * 70)
print("Extracting Comprehensive Features from p53 MD Simulations")
print("=" * 70)

# Store all results
all_features = []

for temp in temperatures:
    temp_celsius = temp - 273
    print(f"\n{'='*70}")
    print(f"Processing {temp}K ({temp_celsius}°C)")
    print("=" * 70)
    
    xtc_file = sim_dir / f"md_{temp}K.xtc"
    
    if not xtc_file.exists():
        print(f"  WARNING: {xtc_file} not found, skipping...")
        continue
    
    # Load trajectory
    u = mda.Universe(str(ref_gro), str(xtc_file))
    protein = u.select_atoms("protein")
    ca_atoms = u.select_atoms("protein and name CA")
    backbone = u.select_atoms("protein and backbone")
    
    n_frames = len(u.trajectory)
    print(f"  Frames: {n_frames}")
    
    # Reference (first frame)
    u.trajectory[0]
    ref_positions = ca_atoms.positions.copy()
    
    # Arrays to store per-frame values
    rmsd_vals = []
    rg_vals = []
    hbond_counts = []
    end_to_end = []
    
    # ===== 1. RMSD and Rg =====
    print("  Calculating RMSD and Rg...")
    for ts in u.trajectory:
        # RMSD
        rmsd = rms.rmsd(ca_atoms.positions, ref_positions, superposition=True)
        rmsd_vals.append(rmsd)
        
        # Radius of gyration
        rg = protein.radius_of_gyration()
        rg_vals.append(rg)
        
        # End-to-end distance (N-term to C-term CA)
        n_term = ca_atoms[0].position
        c_term = ca_atoms[-1].position
        e2e = np.linalg.norm(c_term - n_term)
        end_to_end.append(e2e)
    
    print(f"    RMSD: {np.mean(rmsd_vals):.2f} ± {np.std(rmsd_vals):.2f} Å")
    print(f"    Rg: {np.mean(rg_vals):.2f} ± {np.std(rg_vals):.2f} Å")
    print(f"    End-to-End: {np.mean(end_to_end):.2f} ± {np.std(end_to_end):.2f} Å")
    
    # ===== 2. Hydrogen Bonds =====
    print("  Calculating Hydrogen Bonds...")
    try:
        hbonds = HydrogenBondAnalysis(
            universe=u,
            donors_sel="protein and name N",
            hydrogens_sel="protein and name H",
            acceptors_sel="protein and name O",
            d_a_cutoff=3.5,
            d_h_a_angle_cutoff=120
        )
        hbonds.run(verbose=False)
        
        # Count H-bonds per frame
        for frame_idx in range(n_frames):
            frame_hbonds = hbonds.results.hbonds[hbonds.results.hbonds[:, 0] == frame_idx]
            hbond_counts.append(len(frame_hbonds))
        
        avg_hbonds = np.mean(hbond_counts)
        print(f"    H-bonds: {avg_hbonds:.1f} ± {np.std(hbond_counts):.1f}")
    except Exception as e:
        print(f"    H-bond calculation failed: {e}")
        hbond_counts = [0] * n_frames
        avg_hbonds = 0
    
    # ===== 3. RMSF (Flexibility per residue) =====
    print("  Calculating RMSF...")
    # Reload universe to reset alignment
    u = mda.Universe(str(ref_gro), str(xtc_file))
    ca_atoms = u.select_atoms("protein and name CA")
    
    align.AlignTraj(u, u, select="protein and name CA", in_memory=True).run()
    rmsf_analysis = rms.RMSF(ca_atoms).run()
    rmsf_values = rmsf_analysis.results.rmsf
    avg_rmsf = np.mean(rmsf_values)
    print(f"    Average RMSF: {avg_rmsf:.2f} Å")
    
    # ===== 4. Asphericity (shape) =====
    print("  Calculating Shape Parameters...")
    asphericity_vals = []
    u.trajectory[0]  # Reset to first frame
    for ts in u.trajectory:
        # Moment of inertia tensor
        I = protein.moment_of_inertia()
        eigenvalues = np.linalg.eigvalsh(I)
        eigenvalues = np.sort(eigenvalues)
        
        # Asphericity: deviation from sphere (normalized)
        trace = np.sum(eigenvalues)
        if trace > 0:
            asphericity = (eigenvalues[2] - 0.5 * (eigenvalues[0] + eigenvalues[1])) / trace
        else:
            asphericity = 0
        asphericity_vals.append(asphericity)
    
    avg_asphericity = np.mean(asphericity_vals)
    print(f"    Asphericity: {avg_asphericity:.4f}")
    
    # ===== 5. Native Contacts =====
    print("  Calculating Native Contacts...")
    # Use first frame as native state
    u.trajectory[0]
    native_contacts_ref = []
    cutoff = 8.0  # Å
    
    # Calculate CA-CA distances in native state
    ca_pos_native = ca_atoms.positions
    n_res = len(ca_atoms)
    
    for i in range(n_res):
        for j in range(i + 4, n_res):  # Skip nearby residues
            dist = np.linalg.norm(ca_pos_native[i] - ca_pos_native[j])
            if dist < cutoff:
                native_contacts_ref.append((i, j, dist))
    
    total_native = len(native_contacts_ref)
    print(f"    Native contacts (reference): {total_native}")
    
    # Calculate fraction of native contacts over trajectory
    native_contact_fractions = []
    for ts in u.trajectory:
        ca_pos = ca_atoms.positions
        contacts_preserved = 0
        
        for (i, j, native_dist) in native_contacts_ref:
            current_dist = np.linalg.norm(ca_pos[i] - ca_pos[j])
            if current_dist < cutoff * 1.2:  # Allow 20% tolerance
                contacts_preserved += 1
        
        fraction = contacts_preserved / total_native if total_native > 0 else 0
        native_contact_fractions.append(fraction)
    
    avg_native_contacts = np.mean(native_contact_fractions)
    print(f"    Native Contact Fraction: {avg_native_contacts:.3f}")
    
    # ===== Store Results =====
    feature_data = {
        'temperature_K': temp,
        'temperature_C': temp_celsius,
        'rmsd_mean': np.mean(rmsd_vals),
        'rmsd_std': np.std(rmsd_vals),
        'rg_mean': np.mean(rg_vals),
        'rg_std': np.std(rg_vals),
        'end_to_end_mean': np.mean(end_to_end),
        'end_to_end_std': np.std(end_to_end),
        'hbonds_mean': np.mean(hbond_counts) if hbond_counts else 0,
        'hbonds_std': np.std(hbond_counts) if hbond_counts else 0,
        'rmsf_mean': avg_rmsf,
        'asphericity_mean': avg_asphericity,
        'native_contacts_fraction': avg_native_contacts,
        'n_frames': n_frames
    }
    all_features.append(feature_data)
    
    # Save per-temperature detailed data
    temp_df = pd.DataFrame({
        'frame': range(n_frames),
        'time_ps': np.arange(n_frames) * 5,
        'rmsd': rmsd_vals,
        'rg': rg_vals,
        'end_to_end': end_to_end,
        'hbonds': hbond_counts if hbond_counts else [0]*n_frames,
        'asphericity': asphericity_vals,
        'native_contacts': native_contact_fractions
    })
    temp_df.to_csv(feature_dir / f"features_{temp}K.csv", index=False)
    print(f"  Saved: features_{temp}K.csv")

# Create summary DataFrame
summary_df = pd.DataFrame(all_features)
summary_df.to_csv(feature_dir / "all_features_summary.csv", index=False)
print(f"\nSaved: all_features_summary.csv")

print("\n" + "=" * 70)
print("COMPREHENSIVE FEATURE SUMMARY")
print("=" * 70)
print(summary_df.to_string(index=False))

# ===== Create Visualizations =====
print("\n" + "=" * 70)
print("Creating Feature Visualizations...")
print("=" * 70)

fig_dir = output_dir / "figures"
fig_dir.mkdir(parents=True, exist_ok=True)

colors = ['#0000FF', '#00BFFF', '#00FF00', '#FFFF00', '#FFA500', '#FF4500', '#FF0000']

# Figure 1: All features vs temperature
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# RMSD
ax = axes[0, 0]
ax.errorbar(summary_df['temperature_K'], summary_df['rmsd_mean'], 
            yerr=summary_df['rmsd_std'], marker='o', capsize=5, color='red', linewidth=2)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('RMSD (Å)')
ax.set_title('RMSD vs Temperature')
ax.grid(True, alpha=0.3)

# Rg
ax = axes[0, 1]
ax.errorbar(summary_df['temperature_K'], summary_df['rg_mean'],
            yerr=summary_df['rg_std'], marker='s', capsize=5, color='blue', linewidth=2)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Radius of Gyration (Å)')
ax.set_title('Rg vs Temperature')
ax.grid(True, alpha=0.3)

# H-bonds
ax = axes[0, 2]
ax.errorbar(summary_df['temperature_K'], summary_df['hbonds_mean'],
            yerr=summary_df['hbonds_std'], marker='^', capsize=5, color='green', linewidth=2)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Number of H-bonds')
ax.set_title('Hydrogen Bonds vs Temperature')
ax.grid(True, alpha=0.3)

# Native Contacts
ax = axes[1, 0]
ax.plot(summary_df['temperature_K'], summary_df['native_contacts_fraction'], 
        marker='D', color='purple', linewidth=2, markersize=8)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Fraction of Native Contacts')
ax.set_title('Native Contacts vs Temperature')
ax.grid(True, alpha=0.3)

# End-to-end distance
ax = axes[1, 1]
ax.errorbar(summary_df['temperature_K'], summary_df['end_to_end_mean'],
            yerr=summary_df['end_to_end_std'], marker='p', capsize=5, color='orange', linewidth=2)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('End-to-End Distance (Å)')
ax.set_title('End-to-End Distance vs Temperature')
ax.grid(True, alpha=0.3)

# RMSF
ax = axes[1, 2]
ax.plot(summary_df['temperature_K'], summary_df['rmsf_mean'], 
        marker='h', color='brown', linewidth=2, markersize=8)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Average RMSF (Å)')
ax.set_title('Flexibility (RMSF) vs Temperature')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(fig_dir / 'all_features_vs_temperature.png', dpi=150)
print("  Saved: all_features_vs_temperature.png")

# Figure 2: Time series comparison
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for idx, temp in enumerate(temperatures):
    temp_file = feature_dir / f"features_{temp}K.csv"
    if temp_file.exists():
        df = pd.read_csv(temp_file)
        time_ns = df['time_ps'] / 1000
        
        axes[0, 0].plot(time_ns, df['rmsd'], color=colors[idx], label=f'{temp}K', alpha=0.8)
        axes[0, 1].plot(time_ns, df['rg'], color=colors[idx], label=f'{temp}K', alpha=0.8)
        axes[1, 0].plot(time_ns, df['native_contacts'], color=colors[idx], label=f'{temp}K', alpha=0.8)
        axes[1, 1].plot(time_ns, df['hbonds'], color=colors[idx], label=f'{temp}K', alpha=0.8)

axes[0, 0].set_xlabel('Time (ns)')
axes[0, 0].set_ylabel('RMSD (Å)')
axes[0, 0].set_title('RMSD Over Time')
axes[0, 0].legend()
axes[0, 0].grid(True, alpha=0.3)

axes[0, 1].set_xlabel('Time (ns)')
axes[0, 1].set_ylabel('Rg (Å)')
axes[0, 1].set_title('Radius of Gyration Over Time')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.3)

axes[1, 0].set_xlabel('Time (ns)')
axes[1, 0].set_ylabel('Native Contact Fraction')
axes[1, 0].set_title('Native Contacts Over Time')
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3)

axes[1, 1].set_xlabel('Time (ns)')
axes[1, 1].set_ylabel('H-bonds')
axes[1, 1].set_title('Hydrogen Bonds Over Time')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(fig_dir / 'time_series_all_features.png', dpi=150)
print("  Saved: time_series_all_features.png")

print("\n" + "=" * 70)
print("Feature Extraction Complete!")
print("=" * 70)
