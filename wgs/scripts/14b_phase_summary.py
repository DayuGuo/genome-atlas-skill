#!/usr/bin/env python
"""Summarise read-backed phasing and resolve compound-heterozygote / cis-trans questions at clinically relevant loci."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import subprocess, io, pandas as pd, numpy as np
P=str(P); W=f"{P}/wgs/10_phase"
q=subprocess.run(["bcftools","query","-f","%CHROM\t%POS\t[%GT\t%PS]\n",f"{W}/{SAMPLE}.pass.readphased.vcf.gz"],capture_output=True,text=True).stdout
d=pd.read_csv(io.StringIO(q),sep="\t",header=None,names=["chrom","pos","gt","ps"],dtype={"chrom":str,"ps":str})
het=d[d["gt"].str.contains(r"^(0[|/]1|1[|/]0)$",regex=True)]
ph=het[het["gt"].str.contains(r"\|")]
blocks=ph[ph.ps!="."].groupby(["chrom","ps"]).pos.agg(["min","max","size"])
blocks["len"]=blocks["max"]-blocks["min"]
n50=None; L=np.sort(blocks.len.values)[::-1]; c=np.cumsum(L); n50=L[np.searchsorted(c,c[-1]/2)]
print(f"heterozygous sites {len(het):,}; read-phased {len(ph):,} ({len(ph)/len(het)*100:.1f}%); blocks {len(blocks):,}; block N50 {n50/1e3:.1f} kb; longest {blocks.len.max()/1e6:.2f} Mb")
open(f"{W}/summary.txt","w").write(f"het\t{len(het)}\nphased\t{len(ph)}\nblocks\t{len(blocks)}\nn50_kb\t{n50/1e3:.1f}\nmax_mb\t{blocks.len.max()/1e6:.2f}\n")
# --- cis/trans at loci of interest
LOCI={"UGT1A1 *28 + *80":("2",234668500,234669200),"CFTR":("7",117120017,117308718),"ATP7B":("13",52506805,52585630),
      "NUDT15":("13",48609376,48624627),"CYP2C19":("10",96522463,96612671),"HLA-B":("6",31321649,31324989),"ALDH2":("12",112204691,112247782)}
for name,(c,s,e) in LOCI.items():
    t=subprocess.run(["bcftools","query","-r",f"{c}:{s}-{e}","-i",'GT="het"',"-f","%POS\t%ID\t%REF>%ALT\t[%GT\t%PS]\n",f"{W}/{SAMPLE}.pass.readphased.vcf.gz"],capture_output=True,text=True).stdout.splitlines()
    ps=set(l.split("\t")[4] for l in t if len(l.split("\t"))>4 and l.split("\t")[4]!=".")
    print(f"\n{name}: {len(t)} het sites, {len(ps)} phase block(s)")
    for l in t[:8]: print("   ",l)
