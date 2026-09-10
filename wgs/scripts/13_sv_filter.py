#!/usr/bin/env python
"""Filter Delly SVs (PASS, non-ref, 50bp-5Mb, PRECISE or strong PE support), validate DEL/DUP >=2kb by read depth,
annotate overlapping genes/exons (Ensembl 87 GFF3)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import subprocess, gzip, io, collections, bisect, os
import pandas as pd, numpy as np
PROJ = str(P); W = f"{PROJ}/wgs/08_sv"; MEAN_DP = 28.88
q = subprocess.run(["bcftools", "query", "-i", 'FILTER="PASS" && GT="alt"', "-f",
    "%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/PRECISE\t%INFO/PE\t%INFO/SR\t%INFO/MAPQ\t%INFO/CHR2\t%INFO/POS2\t[%GT\t%GQ\t%RC\t%RCL\t%RCR\t%DR\t%DV\t%RR\t%RV]\n",
    f"{W}/delly/target.sv.bcf"], capture_output=True, text=True, check=True).stdout
sv = pd.read_csv(io.StringIO(q), sep="\t", header=None, names="chrom pos end svtype svlen precise pe sr mapq chr2 pos2 geno gq rc rcl rcr dr dv rr rv".split(), dtype={"chrom": str, "chr2": str}, na_values=".")
sv["precise"] = sv.precise.fillna(0).astype(int)
sv["size"] = np.where(sv.svtype == "BND", np.nan, (sv.end - sv.pos).abs())
main = [str(i) for i in range(1, 23)] + ["X", "Y"]
sv = sv[sv.chrom.isin(main) & (sv.gq >= 20)]
nonbnd = sv[(sv.svtype != "BND") & (sv["size"] >= 50) & (sv["size"] <= 5e6) & ((sv.precise == 1) | ((sv.pe >= 5) & (sv.mapq >= 40)))].copy()
# read-depth validation for DEL/DUP >= 2kb via mosdepth 1kb bins
def region_depth(c, s, e):
    out = subprocess.run(["tabix", f"{PROJ}/wgs/01_qc/depth.regions.bed.gz", f"{c}:{max(1,s)}-{e}"], capture_output=True, text=True).stdout
    v = [float(l.split("\t")[3]) for l in out.splitlines()]
    base = MEAN_DP * (0.5 if c in ("X", "Y") else 1.0)
    return np.median(v) / base if v else np.nan
big = nonbnd[(nonbnd.svtype.isin(["DEL", "DUP"])) & (nonbnd["size"] >= 2000)]
nonbnd["dp_ratio"] = np.nan
nonbnd.loc[big.index, "dp_ratio"] = [region_depth(r.chrom, r.pos, r.end) for r in big.itertuples()]
def rd_ok(r):
    if np.isnan(r.dp_ratio): return "n/a"
    if r.svtype == "DEL": return "ok" if r.dp_ratio < 0.75 else "mismatch"
    return "ok" if r.dp_ratio > 1.25 else "mismatch"
def geno_dp(r):
    if np.isnan(r.dp_ratio): return r.geno
    if r.svtype == "DEL": return "1/1" if r.dp_ratio < 0.2 else "0/1"
    return "1/1" if r.dp_ratio > 1.8 else "0/1"
nonbnd["depth_check"] = nonbnd.apply(rd_ok, axis=1)
nonbnd["geno"] = nonbnd.apply(geno_dp, axis=1)
# drop depth-inconsistent large CNVs and large imprecise inversions
nonbnd = nonbnd[(nonbnd.depth_check != "mismatch") & ~((nonbnd.svtype == "INV") & (nonbnd["size"] > 100000) & ~((nonbnd.precise == 1) & (nonbnd.sr >= 5)))]
# gene annotation
genes = collections.defaultdict(list); exons = collections.defaultdict(list)
with gzip.open(f"{PROJ}/data/ref/annot/Homo_sapiens.GRCh37.87.gff3.gz", "rt") as f:
    for l in f:
        if l.startswith("#"): continue
        p = l.split("\t")
        if p[2] == "gene" and "biotype=protein_coding" in p[8]:
            name = [x for x in p[8].split(";") if x.startswith("Name=")]; genes[p[0]].append((int(p[3]), int(p[4]), name[0][5:] if name else p[8][:20]))
        elif p[2] == "CDS":
            exons[p[0]].append((int(p[3]), int(p[4])))
for d in (genes, exons):
    for c in d: d[c].sort()
def overlap(d, c, s, e, with_name=True):
    hits = []
    for a, b, *n in d.get(c, []):
        if a > e: break
        if b >= s: hits.append(n[0] if n else True)
    return hits
nonbnd["genes"] = [",".join(sorted(set(overlap(genes, r.chrom, r.pos, r.end)))) for r in nonbnd.itertuples()]
nonbnd["cds_overlap"] = [len(overlap(exons, r.chrom, r.pos, r.end)) > 0 for r in nonbnd.itertuples()]
nonbnd["whole_gene_del"] = [",".join(g for a, b, g in genes.get(r.chrom, []) if r.svtype == "DEL" and r.pos <= a and r.end >= b) for r in nonbnd.itertuples()]
nonbnd = nonbnd.sort_values(["chrom", "pos"], key=lambda s: s.map(lambda x: main.index(x) if x in main else x) if s.name == "chrom" else s)
nonbnd.to_csv(f"{W}/sv_filtered.tsv", sep="\t", index=False)
print("filtered non-BND SVs:", len(nonbnd)); print(nonbnd.svtype.value_counts().to_dict())
print("size classes:", pd.cut(nonbnd["size"], [0, 100, 1000, 10000, 100000, 5e6]).value_counts().sort_index().to_dict())
print("depth check for DEL/DUP>=2kb:", nonbnd.depth_check.value_counts().to_dict())
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
cds = nonbnd[nonbnd.cds_overlap & (nonbnd.depth_check != "mismatch")]
print(f"\nSVs overlapping CDS (depth-consistent): {len(cds)}; hom: {(cds.geno=='1/1').sum()}")
print(cds[cds.geno == "1/1"][["chrom", "pos", "end", "svtype", "size", "geno", "dp_ratio", "genes", "whole_gene_del"]].to_string(index=False))
print("\nhet SVs deleting/duplicating whole genes or >=10kb with CDS:")
print(cds[(cds.geno != "1/1") & ((cds.whole_gene_del != "") | (cds["size"] >= 10000))][["chrom", "pos", "end", "svtype", "size", "geno", "dp_ratio", "genes", "whole_gene_del"]].to_string(index=False))
print("\ndp_ratio distribution by type/genotype (DEL/DUP >= 2kb):")
b = nonbnd[nonbnd.dp_ratio.notna()]
print(b.groupby(["svtype", "geno"]).dp_ratio.describe()[["count", "25%", "50%", "75%"]].round(2))
