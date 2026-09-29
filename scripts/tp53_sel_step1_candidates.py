#!/usr/bin/env python3
# Step 1: enumerate every possible TP53 DBD (94-292) missense SNV, annotate SBS96 channel, match CRISPR RFS.
# Validation only. No null model, no permutation.
import os, re, sys, urllib.request, numpy as np, pandas as pd

D = "/home/raza/p53_thermal_xai/data/external/tp53_selection"
os.makedirs(D, exist_ok=True)
CDS_START, CDS_END, CHR = 7571720, 7590868, "17"   # TP53 hg19 locus, minus strand
COMP = {"A":"T","C":"G","G":"C","T":"A"}
CODON = {}
_b = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
_m = "TTTTCCCCAAAAGGGG"*4
_l = "TCAG"*16
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for a, b, c, aa in zip(_b, _m, _l, _aa): CODON[a+b+c] = aa

# --- TP53 RefSeq NM_000546 CDS, hg19 exon coords (minus strand) ---
EXONS = [(7579912,7579839),(7579721,7579700),(7579590,7579312),(7578554,7578371),
         (7578289,7578177),(7577971,7577499),(7577155,7577019),(7576926,7576853),
         (7576853,7576853)]  # placeholder; replaced by Ensembl fetch below

def fetch(url, path, binary=False):
    if os.path.exists(path): return path
    urllib.request.urlretrieve(url, path); return path

print("# fetching TP53 CDS and genomic sequence from Ensembl GRCh37 REST")
try:
    import json
    seq_cds = json.load(urllib.request.urlopen(
        "https://grch37.rest.ensembl.org/sequence/id/ENST00000269305?type=cds;content-type=application/json", timeout=60))["seq"]
    exons = json.load(urllib.request.urlopen(
        "https://grch37.rest.ensembl.org/overlap/id/ENST00000269305?feature=cds;content-type=application/json", timeout=60))
except Exception as e:
    sys.exit("Ensembl GRCh37 REST failed: %s" % e)
ex = sorted([(e["start"], e["end"]) for e in exons if e.get("Parent") == "ENST00000269305"], reverse=True)
print("# CDS length %d nt (%d codons); %d CDS blocks" % (len(seq_cds), len(seq_cds)//3, len(ex)))
prot = "".join(CODON[seq_cds[i:i+3]] for i in range(0, len(seq_cds)-2, 3))
print("# translated protein length %d; starts %s ; residue 175=%s 248=%s 273=%s"
      % (len(prot), prot[:6], prot[174], prot[247], prot[272]))
if not (prot[174] == "R" and prot[247] == "R" and prot[272] == "R"):
    sys.exit("CDS translation does not match p53 — aborting")

# map CDS position -> genomic coordinate (minus strand, descending)
gpos = []
for s, e in ex: gpos.extend(range(e, s-1, -1))
if len(gpos) != len(seq_cds):
    sys.exit("CDS/genomic length mismatch: %d vs %d" % (len(gpos), len(seq_cds)))

# genomic +strand sequence around locus, for trinucleotide context
lo, hi = min(gpos)-2, max(gpos)+2
gseq = urllib.request.urlopen(
    "https://grch37.rest.ensembl.org/sequence/region/human/%s:%d..%d:1?content-type=text/plain" % (CHR, lo, hi),
    timeout=60).read().decode().strip().upper()
def gbase(p): return gseq[p-lo]
chk = "".join(COMP[gbase(p)] for p in gpos[:9])
print("# genomic->CDS check, first 9 nt: %s vs CDS %s" % (chk, seq_cds[:9]))
if chk != seq_cds[:9]: sys.exit("genomic/CDS orientation mismatch")

PY_ = {"C":"C","T":"T"}
rows = []
for ci in range(93, 292):                      # residues 94-292, 0-based
    codon = seq_cds[ci*3:ci*3+3]; wt_aa = CODON[codon]
    for k in range(3):
        cds_i = ci*3 + k; g = gpos[cds_i]
        ref_plus = gbase(g)
        for alt_plus in "ACGT":
            if alt_plus == ref_plus: continue
            alt_cds = COMP[alt_plus]
            new = codon[:k] + alt_cds + codon[k+1:]
            aa = CODON[new]
            if aa == wt_aa or aa == "*": continue
            # SBS96 channel: pyrimidine-referenced on the + strand
            if ref_plus in "CT":
                r, a = ref_plus, alt_plus
                l5, r3 = gbase(g-1), gbase(g+1)
            else:
                r, a = COMP[ref_plus], COMP[alt_plus]
                l5, r3 = COMP[gbase(g+1)], COMP[gbase(g-1)]
            ch = "%s[%s>%s]%s" % (l5, r, a, r3)
            rows.append(dict(residue=ci+1, wt_aa=wt_aa, alt_aa=aa, codon_pos=k+1,
                             gpos=g, ref_plus=ref_plus, alt_plus=alt_plus,
                             sbs96=ch, sub="%s>%s" % (r, a)))
c = pd.DataFrame(rows)
print("\n# candidate DBD missense SNVs: %d over %d residues" % (len(c), c.residue.nunique()))
print("# distinct protein changes: %d | distinct SBS96 channels used: %d/96"
      % (len(c.groupby(['residue','wt_aa','alt_aa'])), c.sbs96.nunique()))

sge = pd.read_excel("/home/raza/p53_thermal_xai/data/external/funk2025/41588_2024_2039_MOESM4_ESM.xlsx",
                    "Supp_Table_2", header=2)
sge = sge[sge.type_p == "mis"].dropna(subset=["codon","aa_ref","aa_alt","rfs_median"])
sge["codon"] = sge.codon.astype(int)
rfs = sge.groupby(["codon","aa_ref","aa_alt"]).rfs_median.median()
c["CRISPR_RFS"] = [rfs.get((r.residue, r.wt_aa, r.alt_aa), np.nan) for r in c.itertuples()]
ref_map = sge.groupby("codon").aa_ref.first()
c["sge_ref"] = c.residue.map(ref_map)
mism = ((c.sge_ref.notna()) & (c.sge_ref != c.wt_aa)).groupby(c.residue).any()
print("# residues where our WT aa disagrees with SGE aa_ref: %d" % int(mism.sum()))
ok = c.CRISPR_RFS.notna()
print("# candidates with CRISPR RFS: %d / %d (%.1f%%); residues covered %d"
      % (ok.sum(), len(c), 100*ok.mean(), c[ok].residue.nunique()))
c.to_csv(D + "/tp53_dbd_candidate_missense_snvs.tsv", sep="\t", index=False)

print("\n# per-substitution-type breakdown of candidates:")
print(c.groupby("sub").size().to_string())
print("\nSaved: " + D + "/tp53_dbd_candidate_missense_snvs.tsv")
