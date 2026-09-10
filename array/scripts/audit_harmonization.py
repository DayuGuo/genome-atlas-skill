#!/usr/bin/env python3
"""Audit: raw QC, allele harmonization vs 1000G (GRCh37) and HGDP (GRCh38 via liftover), ClinVar allele consistency."""
import pathlib, pandas as pd, numpy as np, gzip, re, json
from config import *
ROOT = pathlib.Path(__file__).resolve().parents[1]; AU = AUDIT; AU.mkdir(exist_ok=True)
raw = pd.read_csv(RAW, sep="\t", comment="#", header=None, names=["rsid","chrom","pos","genotype"], dtype={"chrom":str})
raw["genotype"] = raw.genotype.fillna("--").astype(str)
out = {}
out["n_rows"] = len(raw); out["n_rsid"] = int(raw.rsid.str.startswith("rs").sum()); out["n_wegene_id"] = int(raw.rsid.str.startswith("w").sum())
out["n_other_id"] = int((~raw.rsid.str.startswith(("rs","w"))).sum())
out["dup_rsid"] = int(raw.rsid.duplicated().sum()); out["dup_pos"] = int(raw.duplicated(["chrom","pos"]).sum())
out["genotype_alphabet"] = sorted(set("".join(raw.genotype)))
miss = raw.genotype.str.contains("-"); out["n_missing"] = int(miss.sum()); out["call_rate"] = round(1-miss.mean(), 5)
g = raw[~miss]
het = g.genotype.str[0] != g.genotype.str[1]
out["het_autosomal"] = round(het[g.chrom.isin([str(i) for i in range(1,23)])].mean(), 4)
for c in ["X","Y","MT"]:
    sub = g[g.chrom==c]; out[f"{c}_called"] = len(sub); out[f"{c}_het"] = int(het[g.chrom==c].sum())
# PAR check: X het within PAR1 (60001-2699520) / PAR2 (154931044-155260560) GRCh37
par = g[(g.chrom=="X") & ((g.pos<2699520)|(g.pos>154931044))]; out["X_PAR_sites"] = len(par); out["X_PAR_het"] = int((par.genotype.str[0]!=par.genotype.str[1]).sum())
out["Y_called_share"] = round(len(g[g.chrom=="Y"])/len(raw[raw.chrom=="Y"]), 4)
out["indel_sites"] = int(raw.genotype.str.contains("[ID]").sum())
out["per_chrom_missing"] = raw.assign(mis=raw.genotype.str.contains("-")).groupby("chrom").mis.sum().astype(int).to_dict()

# ---- harmonization vs 1000G (results/ancestry/kg.overlap.pvar = biallelic SNPs at the sample autosomal positions)
kg = pd.read_csv(RESULTS/"ancestry"/"kg.overlap.pvar", sep="\t", comment="#", header=None, usecols=[0,1,2,3,4], names=["chrom","pos","id","ref","alt"], dtype={"chrom":str})
d = g[g.chrom.isin([str(i) for i in range(1,23)]) & g.genotype.str.match(r"^[ACGT]{2}$")].copy()
d["a1"], d["a2"] = d.genotype.str[0], d.genotype.str[1]
m = d.merge(kg, on=["chrom","pos"], how="left", indicator="mrg")
COMP = str.maketrans("ACGT","TGCA")
def classify(r):
    if r.mrg != "both": return "not_in_1000G_biallelic_SNP"
    al = {r.ref, r.alt}; gset = {r.a1, r.a2}
    amb = al in ({"A","T"},{"C","G"})
    if gset <= al: return "ambiguous_AT_CG_dropped" if amb else "exact_match"
    if {x.translate(COMP) for x in gset} <= al: return "ambiguous_AT_CG_dropped" if amb else "complement_flip_needed"
    return "allele_mismatch"
m["class"] = [classify(r) for r in m.itertuples()]
cls = m["class"].value_counts().to_dict(); out["harmonization_1000G"] = cls
# among matched non-ambiguous: is observed allele set consistent with REF/ALT?  (REF/ALT swap is not applicable: array has no REF)
mm = m[m["class"]=="allele_mismatch"]; mm.head(50).to_csv(AU/"allele_mismatch_examples_1000G.tsv", sep="\t", index=False)
m[["rsid","chrom","pos","genotype","ref","alt","class"]].to_csv(AU/"allele_harmonization.tsv", sep="\t", index=False)
# how many 1000G-position matches had a different rsID in 1000G (rsID consistency)
both = m[m.mrg=="both"]; has_rs = both[both.id.str.startswith("rs")]
out["rsid_mismatch_vs_1000G"] = int((has_rs.id != has_rs.rsid).sum()); out["rsid_compared"] = len(has_rs)

# ---- ClinVar allele consistency: for sites used in clinvar screen, do genotype alleles ⊆ {REF, ALT}?
cv = pd.read_csv(RESULTS/"clinvar_screen.tsv", sep="\t", dtype={"chrom":str})
ok = [set(gt) <= {r, a} for gt, r, a in zip(cv.genotype, cv.ref, cv.alt)]
out["clinvar_rows"] = len(cv); out["clinvar_geno_subset_of_refalt"] = int(sum(ok)); out["clinvar_geno_has_other_allele"] = int(len(cv)-sum(ok))
cv[~np.array(ok)].to_csv(AU/"clinvar_allele_inconsistent.tsv", sep="\t", index=False)
# ClinVar file date
with gzip.open(REF/"clinvar_grch37.vcf.gz", "rt") as fh:
    for line in fh:
        if line.startswith("##fileDate"): out["clinvar_fileDate"] = line.strip().split("=")[1]
        if not line.startswith("#"): break

# ---- HGDP liftover audit
try:
    hg = pd.read_csv(RESULTS/"hgdp"/"hgdp.eas.pvar", sep="\t", comment="#", header=None, usecols=[0,1,2,3,4], names=["chrom","pos","id","ref","alt"], dtype={"chrom":str})
    out["hgdp_overlap_sites"] = len(hg)
    out["hgdp_vcf_sites_used"] = sum(1 for l in open(RESULTS/"hgdp"/"sample.vcf") if not l.startswith("#"))
except Exception as e: out["hgdp"] = str(e)
json.dump(out, open(AU/"qc_and_harmonization.json","w"), indent=1, ensure_ascii=False)
print(json.dumps(out, indent=1, ensure_ascii=False))
