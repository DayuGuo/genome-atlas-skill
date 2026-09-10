#!/usr/bin/env python
"""The classic 'personality gene' variants, read straight from the WGS, with what the evidence actually supports."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pathlib
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import subprocess, pandas as pd, os
P=str(P); V=f"{P}/wgs/00_input/{SAMPLE}.norm.vcf.gz"; PV=KG_PFILE+".pvar"; W=f"{P}/wgs/20_behaviour"
os.makedirs(W,exist_ok=True)
PANEL = pathlib.Path(__file__).resolve().parents[1] / "panel" / "candidate_genes.tsv"
_cg = pd.read_csv(PANEL, sep="\t").fillna("")
M = [(r.gene, r.rsid, r.variant_zh, r.claim_zh, r.evidence_zh) for r in _cg.itertuples()]
EN = {r.rsid: (r.variant_en, r.claim_en, r.evidence_en) for r in _cg.itertuples()}
pos={}
want=set(x[1] for x in M)
with open(PV) as f:
    for line in f:
        if line.startswith("#"): continue
        c,p,i,r,a=line.rstrip("\n").split("\t")[:5]
        for tok in i.split(";"):
            if tok in want and tok not in pos: pos[tok]=(c,int(p),r,a)
rows=[]
for gene,rs,var,claim,truth in M:
    if rs not in pos: rows.append((gene,rs,var,"-","未在 1000G 面板中",claim,truth)); continue
    c,p,r,a=pos[rs]
    t=subprocess.run(["bcftools","query","-r",f"{c}:{p}-{p}","-f","%REF\t%ALT\t%FILTER\t[%GT\t%DP]\n",V],capture_output=True,text=True).stdout.strip().splitlines()
    if t:
        rec=[x.split("\t") for x in t]; hit=[x for x in rec if x[0]==r and x[1]==a] or rec
        ref,alt,flt,gt,dp=hit[0]; n=gt.replace("|","/").count("1")
        g=(alt if n else ref) if c=="X" else {0:f"{ref}/{ref}",1:f"{ref}/{alt}",2:f"{alt}/{alt}"}.get(n,gt)
        rows.append((gene,rs,var,f"{c}:{p}",g,claim,truth))
    else:
        rows.append((gene,rs,var,f"{c}:{p}",r if c=="X" else f"{r}/{r}",claim,truth))
df=pd.DataFrame(rows,columns=["gene","rsid","variant","locus","genotype","popular_claim","what_evidence_supports"])
df["variant_en"]=[EN.get(r,("","",""))[0] for r in df.rsid]
df["popular_claim_en"]=[EN.get(r,("","",""))[1] for r in df.rsid]
df["evidence_en"]=[EN.get(r,("","",""))[2] for r in df.rsid]
df.to_csv(f"{W}/candidate_genes.tsv",sep="\t",index=False)
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",34)
print(df[["gene","rsid","variant","genotype","what_evidence_supports"]].to_string(index=False))
