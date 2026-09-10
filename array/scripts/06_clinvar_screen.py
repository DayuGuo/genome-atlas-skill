#!/usr/bin/env python3
"""Screen the user's genotypes against ClinVar (GRCh37).

For every array site where the sample carries at least one allele that matches a ClinVar ALT
allele, report the ClinVar significance, review status, condition, gene. Sorted so that
Pathogenic / Likely pathogenic come first. Output: results/clinvar_screen.tsv
NOTE: array genotyping has a non-trivial error rate at rare sites; any P/LP hit needs
confirmation by sequencing before it means anything.
"""
import gzip, pathlib, re, pandas as pd
from config import *
ROOT = pathlib.Path(__file__).resolve().parents[1]
df = pd.read_parquet(DATA/"genome.parquet")
df = df[~df.is_missing]
key = {(c, p): (r, g) for c, p, r, g in zip(df.chrom.astype(str), df.pos, df.rsid, df.genotype)}

def info_field(info, k):
    m = re.search(rf"(?:^|;){k}=([^;]*)", info); return m.group(1) if m else ""

rows = []
with gzip.open(REF/"clinvar_grch37.vcf.gz", "rt") as fh:
    for line in fh:
        if line[0] == "#": continue
        chrom, pos, vid, ref, alt, _, _, info = line.rstrip("\n").split("\t")[:8]
        k = (chrom, int(pos))
        if k not in key: continue
        rsid, geno = key[k]
        if len(ref) != 1 or len(alt) != 1:  # SNV only; array indel encoding (I/D) is ambiguous
            continue
        n_alt = geno.count(alt)
        if n_alt == 0: continue
        sig = info_field(info, "CLNSIG"); rev = info_field(info, "CLNREVSTAT")
        rows.append(dict(rsid=rsid, chrom=chrom, pos=int(pos), ref=ref, alt=alt, genotype=geno, n_alt=n_alt,
                         clinvar_id=vid, significance=sig, review_status=rev,
                         gene=info_field(info, "GENEINFO").split(":")[0],
                         condition=info_field(info, "CLNDN").replace("_", " ")[:150],
                         hgvs=info_field(info, "CLNHGVS")))
out = pd.DataFrame(rows)
rank = {"Pathogenic": 0, "Pathogenic/Likely_pathogenic": 0, "Likely_pathogenic": 1, "Pathogenic|risk_factor": 0,
        "risk_factor": 2, "drug_response": 2, "Affects": 2, "association": 3, "protective": 3,
        "Conflicting_classifications_of_pathogenicity": 4, "Uncertain_significance": 5,
        "Likely_benign": 8, "Benign/Likely_benign": 8, "Benign": 9}
out["rank"] = out.significance.map(lambda s: min([rank.get(x, 6) for x in re.split(r"[|,/]", s)] + [6]))
out = out.sort_values(["rank", "chrom", "pos"]).reset_index(drop=True)
out.to_csv(RESULTS/"clinvar_screen.tsv", sep="\t", index=False)
print("array sites with a ClinVar ALT allele carried:", len(out))
print(out.significance.str.split(r"[|/]").str[0].value_counts().head(12).to_string())
top = out[out["rank"] <= 2]
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
print("\n=== Pathogenic / Likely pathogenic / risk factor / drug response hits ===")
print(top[["rsid","gene","genotype","ref","alt","n_alt","significance","review_status","condition"]].to_string())
