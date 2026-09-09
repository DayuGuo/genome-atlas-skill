#!/usr/bin/env python3
"""Convert the WeGene genotype export into analysis-friendly formats.

Outputs (all under results/ or data/):
  data/genome.parquet        columnar table: rsid, chrom, pos, genotype, a1, a2, is_het, is_missing
  data/genome.duckdb         DuckDB database with table `snp` (same columns) for SQL queries
  data/genome.23andme.txt    23andMe-style text (rsid, chrom, pos, genotype) for third-party tools
  results/summary_basic.tsv  per-chromosome counts, missingness, heterozygosity
"""
import sys, pathlib
import pandas as pd, duckdb

from config import *
ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = RAW

df = pd.read_csv(RAW, sep="\t", comment="#", header=None,
                 names=["rsid", "chrom", "pos", "genotype"], dtype={"chrom": str})
df["genotype"] = df["genotype"].fillna("--").astype(str)
df["a1"] = df["genotype"].str[0]
df["a2"] = df["genotype"].str[1]
df["is_missing"] = df["genotype"].str.contains("-")
df["is_het"] = (~df["is_missing"]) & (df["a1"] != df["a2"])

CHR_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
df["chrom"] = pd.Categorical(df["chrom"], categories=CHR_ORDER, ordered=True)
df = df.sort_values(["chrom", "pos"]).reset_index(drop=True)

df.to_parquet(DATA / "genome.parquet", index=False)

con = duckdb.connect(str(DATA / "genome.duckdb"))
con.execute("DROP TABLE IF EXISTS snp")
con.execute("CREATE TABLE snp AS SELECT rsid, CAST(chrom AS VARCHAR) AS chrom, pos, genotype, a1, a2, is_het, is_missing FROM df")
con.execute("CREATE INDEX IF NOT EXISTS idx_rsid ON snp(rsid)")
con.close()

# 23andMe-style file (widely accepted by third-party tools)
out = DATA / "genome.23andme.txt"
with open(out, "w") as fh:
    fh.write("# rsid\tchromosome\tposition\tgenotype\n")
    df[["rsid", "chrom", "pos", "genotype"]].to_csv(fh, sep="\t", header=False, index=False)

# summary
called = df[~df.is_missing]
summ = df.groupby("chrom", observed=True).agg(
    n_sites=("rsid", "size"),
    n_missing=("is_missing", "sum"),
    n_het=("is_het", "sum"),
).reset_index()
summ["call_rate"] = 1 - summ.n_missing / summ.n_sites
summ["het_rate"] = summ.n_het / (summ.n_sites - summ.n_missing)
summ.to_csv(RESULTS / "summary_basic.tsv", sep="\t", index=False, float_format="%.4f")
print(summ.to_string(index=False))
print(f"\nTotal sites: {len(df):,}  called: {len(called):,}  autosomal het rate: "
      f"{called[called.chrom.isin([str(i) for i in range(1,23)])].is_het.mean():.4f}")
