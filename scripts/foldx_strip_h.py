import MDAnalysis as mda
import glob, os

SRC = "/home/raza/p53_thermal_xai/analysis/foldx/snapshots_all"
OUT = "/home/raza/p53_thermal_xai/analysis/foldx/snapshots_all_noH"
os.makedirs(OUT, exist_ok=True)

files = sorted(glob.glob(f"{SRC}/*.pdb"))
for f in files:
    u = mda.Universe(f)
    heavy = u.select_atoms('not name H*')
    out = f.replace(SRC, OUT)
    heavy.write(out)

print('done:', len(glob.glob(f"{OUT}/*.pdb")), 'files')
