#!/usr/bin/env python3
"""Build an HSD file of mtDNA variants vs rCRS and run haplogrep3."""
import pathlib, subprocess, pandas as pd
from config import *
ROOT = pathlib.Path(__file__).resolve().parents[1]
rcrs = "".join(l.strip() for l in open(REF/"rCRS.fasta") if not l.startswith(">"))
df = pd.read_parquet(DATA/"genome.parquet")
mt = df[(df.chrom=="MT") & (~df.is_missing)].copy()
mt["ref"] = [rcrs[p-1] for p in mt.pos]
mt["diff"] = mt.a1 != mt.ref
# indel-style ids (a1 may be I/D) -> skip those for haplogrep
snv = mt[mt.a1.isin(list("ACGT"))]
variants = [f"{p}{a}" for p, a in zip(snv[snv["diff"]].pos, snv[snv["diff"]].a1)]
covered = sorted(set(snv.pos))
# range: list covered positions as ranges (haplogrep uses this to know what was typed)
ranges, start, prev = [], covered[0], covered[0]
for p in covered[1:]:
    if p != prev+1: ranges.append(f"{start}-{prev}" if start!=prev else f"{start}"); start=p
    prev=p
ranges.append(f"{start}-{prev}" if start!=prev else f"{start}")
hsd = DATA/"sample.MT.hsd"
with open(hsd,"w") as f:
    f.write("SampleId\tRange\tHaplogroup\tPolymorphisms\n")
    f.write(f"{SAMPLE}\t" + ";".join(ranges) + "\t?\t" + "\t".join(variants) + "\n")
print(f"MT sites typed: {len(snv)}, differing from rCRS: {len(variants)}")
print("Variants:", " ".join(variants))
mt.to_csv(RESULTS/"mt_sites_vs_rCRS.tsv", sep="\t", index=False)
out = RESULTS/"mt_haplogrep3"
out.mkdir(exist_ok=True)
subprocess.run([str(pathlib.Path(HAPLOGREP3)), "classify", "--tree", "phylotree-rcrs@17.2",
                "--in", str(hsd), "--out", str(out/"haplogroups.txt"), "--extend-report"], check=True)
print(open(out/"haplogroups.txt").read())
