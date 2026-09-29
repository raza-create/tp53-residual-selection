#!/usr/bin/env python3
"""Extract features from 100ns MD simulations"""
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

sim_dir = Path("/home/raza/p53_thermal_xai/data/structures/long_simulations")
output_dir = Path("/home/raza/p53_thermal_xai/data/features_100ns")
output_dir.mkdir(parents=True, exist_ok=True)

ref_gro = sim_dir.parent / "npt.gro"
temperatures = [248, 273, 310, 350, 400, 450]

print("="*70)
print("EXTRACTING FEATURES FROM 100ns SIMULATIONS")
print("="*70)

all_data = []

for temp in temperatures:
    print(f"\nProcessing {temp}K...")
    xtc = sim_dir / f"md_100ns_{temp}K.xtc"
    
    if not xtc.exists():
        print(f"  WARNING: {xtc} not found!")
        continue
    
    u = mda.Universe(str(ref_gro), str(xtc))
    protein = u.select_atoms("protein")
    ca = u.select_atoms("protein and name CA")
    
    print(f"  Frames: {len(u.trajectory)}")
    
    # Align to first frame
    align.AlignTraj(u, u, select="protein and name CA", in_memory=True).run()
    
    # Reference for RMSD
    ref_pos = ca.positions.copy()
    u.trajectory[0]
    first_pos = ca.positions.copy()
    
    results = []
    for ts in u.trajectory:
        # RMSD
        rmsd_val = rms.rmsd(ca.positions, first_pos, superposition=False)
        
        # Radius of gyration
        rg = protein.radius_of_gyration()
        
        # End-to-end distance
        n_term = protein.select_atoms("name CA")[0].position
        c_term = protein.select_atoms("name CA")[-1].position
        end_to_end = np.linalg.norm(c_term - n_term)
        
        # H-bonds (simplified - count close N-O pairs)
        n_atoms = protein.select_atoms("name N")
        o_atoms = protein.select_atoms("name O")
        hbonds = 0
        for n in n_atoms.positions:
            dists = np.linalg.norm(o_atoms.positions - n, axis=1)
            hbonds += np.sum(dists < 3.5)
        
        # Native contacts
        if ts.frame == 0:
            ref_ca_pos = ca.positions.copy()
            ref_dists = []
            for i in range(len(ref_ca_pos)):
                for j in range(i+4, len(ref_ca_pos)):
                    d = np.linalg.norm(ref_ca_pos[i] - ref_ca_pos[j])
                    if d < 8.0:
                        ref_dists.append((i, j, d))
        
        native = 0
        for i, j, d0 in ref_dists:
            d = np.linalg.norm(ca.positions[i] - ca.positions[j])
            if d < d0 * 1.2:
                native += 1
        native_frac = native / len(ref_dists) if ref_dists else 1.0
        
        # Asphericity
        positions = protein.positions
        com = positions.mean(axis=0)
        centered = positions - com
        gyration_tensor = np.dot(centered.T, centered) / len(positions)
        eigenvalues = np.linalg.eigvalsh(gyration_tensor)
        eigenvalues = np.sort(eigenvalues)[::-1]
        asphericity = eigenvalues[0] - 0.5*(eigenvalues[1] + eigenvalues[2])
        
        results.append({
            'frame': ts.frame,
            'time_ns': ts.time / 1000,
            'temperature': temp,
            'rmsd': rmsd_val,
            'rg': rg,
            'end_to_end': end_to_end,
            'hbonds': hbonds,
            'native_contacts': native_frac,
            'asphericity': asphericity
        })
    
    df = pd.DataFrame(results)
    df.to_csv(output_dir / f"features_100ns_{temp}K.csv", index=False)
    print(f"  Saved: features_100ns_{temp}K.csv ({len(df)} frames)")
    all_data.append(df)

# Combine all
all_df = pd.concat(all_data, ignore_index=True)
all_df.to_csv(output_dir / "all_features_100ns.csv", index=False)
print(f"\n✅ Total: {len(all_df)} frames saved to all_features_100ns.csv")
print("="*70)
