#!/usr/bin/env python3
"""Look up the curated trait/pharmacogenomic panel in the user's genotypes.

Allele orientation is verified against Ensembl GRCh37 (plus strand): if the panel's
effect allele is not among the Ensembl alleles at that rsid, the row is flagged.
Output: results/trait_panel_report.tsv
"""
import json, pathlib, requests, pandas as pd, duckdb
from config import *
ROOT = pathlib.Path(__file__).resolve().parents[1]
panel = pd.read_csv(PANEL/"trait_panel.tsv", sep="\t")
con = duckdb.connect(str(DATA/"genome.duckdb"), read_only=True)
g = con.execute("select rsid, chrom, pos, genotype from snp where rsid in (select unnest(?))",
                [panel.rsid.tolist()]).df()
m = panel.merge(g, on="rsid", how="left")

# Ensembl GRCh37 batch lookup (cache to disk)
cache = RESULTS/"ensembl_grch37_variation_cache.json"
ens = json.load(open(cache)) if cache.exists() else {}
todo = [r for r in panel.rsid if r not in ens]
for i in range(0, len(todo), 100):
    r = requests.post("https://grch37.rest.ensembl.org/variation/homo_sapiens",
                      headers={"Content-Type": "application/json", "Accept": "application/json"},
                      json={"ids": todo[i:i+100]}, timeout=120)
    r.raise_for_status(); ens.update(r.json())
json.dump(ens, open(cache, "w"))

COMP = str.maketrans("ACGT", "TGCA")
rows = []
for _, r in m.iterrows():
    e = ens.get(r.rsid, {})
    allele_str = ""; ens_pos = None
    for mp in e.get("mappings", []):
        if mp.get("assembly_name") == "GRCh37" and mp.get("seq_region_name") in [str(i) for i in range(1,23)]+["X","Y","MT"]:
            allele_str = mp.get("allele_string", ""); ens_pos = mp.get("start")
    alleles = allele_str.split("/") if allele_str else []
    eff = r.effect_allele
    if eff in ("I", "D"):
        strand_ok = "indel"
    elif not alleles:
        strand_ok = "no-ensembl"
    elif eff in alleles:
        strand_ok = "ok"
    elif eff.translate(COMP) in alleles:
        strand_ok = "FLIP?"
    else:
        strand_ok = "mismatch"
    geno = r.genotype if isinstance(r.genotype, str) else None
    if geno is None:
        n_eff = None; interp = "芯片未覆盖"
    elif "-" in geno:
        n_eff = None; interp = "未检出 (no call)"
    else:
        n_eff = geno.count(eff)
        # X-linked male: hemizygous; count effect alleles as 0 or 2 for interpretation
        if r.chrom == "X" and geno[0] == geno[1]:
            n_eff = 2 if geno[0] == eff else 0
        interp = [r.interp_0, r.interp_1, r.interp_2][n_eff]
    rows.append(dict(rsid=r.rsid, gene=r.gene, category=r.category, trait=r.trait,
                     chrom=r.chrom, pos=r.pos, genotype=geno, effect_allele=eff, n_effect=n_eff,
                     ensembl_alleles=allele_str, strand_check=strand_ok,
                     ensembl_maf=e.get("MAF"), ancestral=e.get("ancestral_allele"),
                     clinical=";".join(e.get("clinical_significance", [])), interpretation=interp))
out = pd.DataFrame(rows)
out.to_csv(RESULTS/"trait_panel_report.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(out[["rsid","gene","genotype","effect_allele","n_effect","ensembl_alleles","strand_check","ancestral","interpretation"]].to_string())
print("\nstrand check summary:", out.strand_check.value_counts().to_dict())
