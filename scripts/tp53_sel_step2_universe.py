#!/usr/bin/env python3
# Step 2: restrict candidate universe and observed TCGA events to residues 126-292.
# No null model, no permutation.
import json, re, urllib.request, numpy as np, pandas as pd

D   = "/home/raza/p53_thermal_xai/data/external/tp53_selection"
SGE = "/home/raza/p53_thermal_xai/data/external/funk2025/41588_2024_2039_MOESM4_ESM.xlsx"
SBSF = "/home/raza/p53_thermal_xai/data/external/gse182927"   # not used here
LO, HI = 126, 292
CANCERS = {"lusc":"Lung-SCC","ov":"Ovary-AdenoCa","hnsc":"Head-SCC","brca":"Breast-cancer",
           "coadread":"ColoRect-AdenoCa","luad":"Lung-AdenoCa","stad":"Stomach-AdenoCa",
           "blca":"Transitional-cell-carcinoma","ucec":"Uterus-AdenoCa","gbm":"CNS-GBM",
           "lihc":"Liver-HCC","paad":"Panc-AdenoCa"}

c = pd.read_csv(D + "/tp53_dbd_candidate_missense_snvs.tsv", sep="\t")
u = c[(c.residue >= LO) & (c.residue <= HI)].copy()
has = u.CRISPR_RFS.notna()
print("Primary universe %d-%d\n" % (LO, HI))
print("Residues:                        %d" % u.residue.nunique())
print("Possible SNV missense mutations: %d" % len(u))
print("With CRISPR RFS:                 %d" % int(has.sum()))
print("Coverage %%:                      %.2f" % (100*has.mean()))
print("Distinct protein substitutions:  %d" % u[has].groupby(["residue","wt_aa","alt_aa"]).ngroups)
print("SBS96 channels represented:      %d/96" % u[has].sbs96.nunique())

miss = u[~has]
print("\nCandidate SNVs in %d-%d still lacking RFS: %d" % (LO, HI, len(miss)))
if len(miss):
    print(miss.groupby(["residue","wt_aa","alt_aa"]).size().reset_index(name="n")
            .assign(sub=lambda d: d.wt_aa+d.residue.astype(str)+d.alt_aa)[["sub","n"]].to_string(index=False))
    print("affected residues: %s" % sorted(miss.residue.unique()))
u.to_csv(D + "/tp53_universe_126_292.tsv", sep="\t", index=False)

sbs = pd.read_csv("/home/raza/p53_thermal_xai/data/external/pcawg_signatures/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
sbs["patient"] = sbs["Sample Names"].str.slice(0, 12)
sge = pd.read_excel(SGE, "Supp_Table_2", header=2)
sge = sge[sge.type_p == "mis"].dropna(subset=["codon","aa_ref","aa_alt","rfs_median"])
sge["codon"] = sge.codon.astype(int)
rfs = sge.groupby(["codon","aa_ref","aa_alt"]).rfs_median.median()

print("\nCancer      TP53_missense_SNV_126_292  with_matched_SBS  with_CRISPR_RFS  final_usable")
obs, tot = [], np.zeros(4, int)
for ct, sig in CANCERS.items():
    st = ct + "_tcga_pan_can_atlas_2018"
    body = json.dumps({"entrezGeneIds":[7157], "sampleListId": st+"_sequenced"}).encode()
    req = urllib.request.Request(
        "https://www.cbioportal.org/api/molecular-profiles/%s_mutations/mutations/fetch?projection=DETAILED" % st,
        data=body, headers={"Content-Type":"application/json"})
    m = json.load(urllib.request.urlopen(req, timeout=120))
    pats = set(sbs[sbs["Cancer Types"] == sig].patient)
    n_in = n_sbs = n_rfs = 0
    for x in m:
        if x.get("mutationType") != "Missense_Mutation": continue
        mm = re.match(r"^([A-Z])(\d+)([A-Z])$", x.get("proteinChange") or "")
        if not mm: continue
        wt, pos, alt = mm.group(1), int(mm.group(2)), mm.group(3)
        if not (LO <= pos <= HI): continue
        if not (len(x.get("referenceAllele","")) == 1 == len(x.get("variantAllele",""))): continue
        n_in += 1
        insbs = x["patientId"] in pats
        if insbs: n_sbs += 1
        v = rfs.get((pos, wt, alt), np.nan)
        if not np.isnan(v): n_rfs += 1
        if insbs and not np.isnan(v):
            obs.append(dict(cancer=ct.upper(), patient=x["patientId"], sample=x["sampleId"],
                            residue=pos, wt_aa=wt, alt_aa=alt, CRISPR_RFS=v))
    fin = sum(1 for o in obs if o["cancer"] == ct.upper())
    tot += np.array([n_in, n_sbs, n_rfs, fin])
    print("%-10s %8d %20d %17d %14d" % (ct.upper(), n_in, n_sbs, n_rfs, fin))
print("%-10s %8d %20d %17d %14d" % ("TOTAL", *tot))
pd.DataFrame(obs).to_csv(D + "/tcga_observed_events_126_292.tsv", sep="\t", index=False)
print("\nSaved: %s/tp53_universe_126_292.tsv\n       %s/tcga_observed_events_126_292.tsv" % (D, D))
