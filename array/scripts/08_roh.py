#!/usr/bin/env python3
"""Runs of homozygosity (ROH) on autosomes: consecutive homozygous called SNPs spanning
>= 1.5 Mb with >= 100 SNPs, allowing 1 heterozygous call per run (genotyping error) and breaking runs at
any >500 kb gap between array sites (centromeres, assembly gaps).
Long ROH (> 5-10 Mb) indicate recent parental relatedness; a handful of 1.5-3 Mb runs is
normal in any East Asian genome."""
import pathlib, pandas as pd
from config import *
ROOT = pathlib.Path(__file__).resolve().parents[1]
df = pd.read_parquet(DATA/"genome.parquet")
df = df[(~df.is_missing) & df.chrom.isin([str(i) for i in range(1,23)]) & df.a1.isin(list("ACGT"))]
runs = []
for c, g in df.groupby("chrom", observed=True):
    pos = g.pos.to_numpy(); het = g.is_het.to_numpy()
    i = 0; n = len(pos)
    while i < n:
        if het[i]: i += 1; continue
        j = i; hets = 0
        while j + 1 < n:
            if pos[j+1] - pos[j] > 500_000: break  # centromere / assembly gap
            if het[j+1]:
                if hets >= 1: break
                hets += 1
            j += 1
        while het[j]: j -= 1
        length = pos[j] - pos[i]; nsnp = j - i + 1
        if length >= 1.5e6 and nsnp >= 100:
            runs.append(dict(chrom=c, start=int(pos[i]), end=int(pos[j]), length_mb=round(length/1e6, 2), n_snps=nsnp))
        i = j + 1
r = pd.DataFrame(runs).sort_values("length_mb", ascending=False)
r.to_csv(RESULTS/"roh.tsv", sep="\t", index=False)
tot = r.length_mb.sum() if len(r) else 0
print(f"ROH segments >=1.5Mb: {len(r)}   total: {tot:.1f} Mb   longest: {r.length_mb.max() if len(r) else 0} Mb   >5Mb: {(r.length_mb>5).sum()}")
print(r.head(15).to_string(index=False))
