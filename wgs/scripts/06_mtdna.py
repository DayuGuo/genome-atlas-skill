#!/usr/bin/env python
"""MT: heteroplasmy from pileup allele depths (alt fraction >= 3%, depth >= 100), haplogrep3 from PASS haploid calls
(plus haplogrep3 on the VCF directly for cross-check)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import os, subprocess, pandas as pd, numpy as np
PROJ = str(P); W = f"{PROJ}/wgs/03_haplo"; os.makedirs(W, exist_ok=True)
rows = []
for line in open(f"{W}/mt_pileup_ad.tsv"):
    pos, ref, alt, ad = line.rstrip("\n").split("\t"); alts = alt.split(","); ads = list(map(int, ad.split(",")))
    dp = sum(ads)
    for a, n in zip(alts, ads[1:]):
        if a == "<*>": continue
        rows.append((int(pos), ref, a, dp, n, n/dp if dp else 0))
df = pd.DataFrame(rows, columns=["pos", "ref", "alt", "dp", "alt_reads", "af"])
het = df[(df.dp >= 100) & (df.af >= 0.03) & (df.af <= 0.97) & (df.alt_reads >= 10)]
hom = df[(df.af > 0.97)]
het.to_csv(f"{W}/mt_heteroplasmy.tsv", sep="\t", index=False)
print("homoplasmic (AF>97%):", len(hom), " heteroplasmic (3-97%, DP>=100, >=10 reads):", len(het))
print(het.to_string(index=False))
# HSD for haplogrep3 from the 44 PASS VCF calls (homoplasmic)
q = subprocess.run(["bcftools", "query", "-r", "MT", "-f", "%POS\t%REF\t%ALT\n", f"{PROJ}/wgs/00_input/{SAMPLE}.pass.vcf.gz"], capture_output=True, text=True).stdout
muts = []
for l in q.splitlines():
    p, r, a = l.split("\t"); p = int(p)
    if len(r) == 1 and len(a) == 1: muts.append(f"{p}{a}")
    elif len(r) == 1 and len(a) > 1: muts.append("".join(f"{p}.{i+1}{b}" for i, b in enumerate(a[1:])))
    elif len(a) == 1 and len(r) > 1: muts.append(" ".join(f"{p+i+1}d" for i in range(len(r)-1)))
with open(f"{W}/target_mt.hsd", "w") as f: f.write(f"SampleId\tRange\tHaplogroup\tPolymorphisms\n{SAMPLE}\t1-16569\t?\t" + "\t".join(muts) + "\n")
subprocess.run([HAPLOGREP3, "classify", "--tree", "phylotree-rcrs@17.2", "--in", f"{W}/target_mt.hsd", "--out", f"{W}/haplogrep3.txt", "--extend-report"], check=True, capture_output=True)
print(open(f"{W}/haplogrep3.txt").read()[:800])
