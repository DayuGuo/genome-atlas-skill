#!/usr/bin/env python
"""Build {NAME_EN}'s genotype at every 1000G biallelic SNV (autosomes + X): PASS WGS call, hom-ref fill in callable
regions, otherwise omitted (= missing). Output VCF -> plink2 pgen with 1000G variant IDs."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import os, sys, subprocess, gzip, io
import numpy as np, pandas as pd

PROJ = os.environ.get("PROJ", str(P)); WGS = f"{PROJ}/wgs"
OUT = f"{WGS}/02_complete"; os.makedirs(OUT, exist_ok=True)
PVAR = f"{PROJ}/data/ref/all_phase3.pvar"
XVCF = f"{WGS}/00_input/X.recall.vcf.gz"
CHRS = [str(c) for c in range(1, 23)] + (["X"] if os.path.exists(XVCF) else [])

pv = pd.read_csv(PVAR, sep="\t", comment=None, skiprows=lambda i: False, header=None, dtype=str,
                 names=["chrom", "pos", "id", "ref", "alt"], engine="c", na_filter=False)
pv = pv[~pv.chrom.str.startswith("#")]
pv = pv[(pv.ref.str.len() == 1) & (pv.alt.str.len() == 1) & pv.ref.isin(list("ACGT")) & pv.alt.isin(list("ACGT"))]
pv["pos"] = pv.pos.astype(np.int64)
print("1000G biallelic SNVs:", len(pv), file=sys.stderr)

call = {}
for line in open(f"{WGS}/00_input/callable.bed"):
    c, s, e = line.split()[:3]; call.setdefault(c, []).append((int(s), int(e)))
def is_callable(c, pos):
    iv = np.array(call[c]); s, e = iv[:, 0], iv[:, 1]
    i = np.searchsorted(s, pos - 1, side="right") - 1
    ok = i >= 0; i = np.clip(i, 0, len(s) - 1)
    return ok & (s[i] <= pos - 1) & (pos - 1 < e[i])

def q(args): return subprocess.run(args, capture_output=True, text=True, check=True).stdout

out = gzip.open(f"{OUT}/{SAMPLE}.1kg_sites.vcf.gz", "wt", compresslevel=3)
out.write("##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
for c in CHRS: out.write(f"##contig=<ID={c}>\n")
out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}\n")
stats = []
for c in CHRS:
    sites = pv[pv.chrom == c]
    src = XVCF if c == "X" else f"{WGS}/00_input/{SAMPLE}.norm.vcf.gz"
    w = pd.read_csv(io.StringIO(q(["bcftools", "query", "-r", c, "-f", "%POS\t%REF\t%ALT\t%FILTER\t[%GT]\n", src])),
                    sep="\t", header=None, names=["pos", "ref", "alt", "flt", "gt"], dtype={"pos": np.int64}, na_filter=False)
    ok = w[(w.flt == "PASS") | (c == "X")]  # X recall was already QUAL/DP filtered
    m = sites.merge(ok[["pos", "ref", "alt", "gt"]], on=["pos", "ref", "alt"], how="left")
    has_any = m.pos.isin(set(w.pos.values))
    fill = m["gt"].isna() & ~has_any & is_callable(c, m.pos.values)
    m.loc[fill, "gt"] = "0/0"
    if c == "X":  # haploid non-PAR calls -> keep as single allele (plink2 handles haploid X for male)
        pass
    keep = m[m["gt"].notna()].copy()
    keep["gt"] = keep["gt"].str.replace("|", "/", regex=False)
    stats.append((c, len(sites), int((m["gt"].notna() & ~fill).sum()), int(fill.sum()), int((m["gt"].isna()).sum())))
    print(c, stats[-1], file=sys.stderr)
    body = keep.pos.astype(str) + "\t" + keep.id + "\t" + keep.ref + "\t" + keep.alt + "\t.\tPASS\t.\tGT\t" + keep["gt"]
    out.write("\n".join((c + "\t" + body).tolist()) + "\n")
out.close()
pd.DataFrame(stats, columns=["chrom", "sites_1kg", "wgs_called", "homref_filled", "missing"]).to_csv(f"{OUT}/complete_set_stats.tsv", sep="\t", index=False)
