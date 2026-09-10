#!/usr/bin/env python
"""Extended blood-group typing. Marker positions are resolved from the 1000G phase3 pvar by rsID, so no
coordinate is taken on trust; genotypes come from {NAME_EN}'s normalised VCF plus the callable mask."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pathlib
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import subprocess, io, os, bisect, collections, pandas as pd, pysam
P=str(P); V=f"{P}/wgs/00_input/{SAMPLE}.norm.vcf.gz"; W=f"{P}/wgs/16_panels"; os.makedirs(W,exist_ok=True)
PANEL = pathlib.Path(__file__).resolve().parents[1] / "panel" / "blood_groups.tsv"
MARK = [(r.system, r.rsid, r.variant, r.meaning) for r in pd.read_csv(PANEL, sep="\t").fillna("").itertuples()]
pv=KG_PFILE+".pvar"
want=set(m[1] for m in MARK)
pos={}
with open(pv) as f:
    for line in f:
        if line.startswith("#"): continue
        c,p,i,r,a=line.rstrip("\n").split("\t")[:5]
        for tok in i.split(";"):
            if tok in want and tok not in pos: pos[tok]=(c,int(p),r,a)
print("resolved via 1000G pvar:",len(pos),"of",len(want))
# callable mask
call=collections.defaultdict(list)
for l in open(f"{P}/wgs/00_input/callable.bed"):
    c,s,e=l.split()[:3]; call[c].append((int(s),int(e)))
starts={c:[s for s,_ in v] for c,v in call.items()}
def callable_(c,p):
    v=call.get(c)
    if not v: return False
    i=bisect.bisect_right(starts[c],p-1)-1
    return i>=0 and v[i][0]<=p-1<v[i][1]
fa=pysam.FastaFile(FASTA)
rows=[]
for sysname,rs,desc,note in MARK:
    if rs not in pos: rows.append((sysname,rs,desc,"-","not in 1000G panel","",note)); continue
    c,p,r,a=pos[rs]
    t=subprocess.run(["bcftools","query","-r",f"{c}:{p}-{p}","-f","%REF\t%ALT\t%FILTER\t[%GT\t%DP\t%AD]\n",V],capture_output=True,text=True).stdout.strip().splitlines()
    if t:
        rec=[x.split("\t") for x in t]
        hit=[x for x in rec if x[0]==r and x[1]==a] or rec
        ref,alt,flt,gt,dp,ad=hit[0]
        n=gt.replace("|","/").count("1"); geno={0:f"{r}/{r}",1:f"{r}/{a}",2:f"{a}/{a}"}.get(n,gt) if len(gt.replace("|","/").split("/"))==2 else (a if n else r)
        rows.append((sysname,rs,desc,f"{c}:{p}",geno,f"{flt} DP={dp} AD={ad}",note))
    else:
        ok=callable_(c,p); b=fa.fetch(c,p-1,p).upper()
        rows.append((sysname,rs,desc,f"{c}:{p}",f"{b}/{b}" if ok else "no call (low depth)","hom-ref (callable)" if ok else "uncallable",note))
df=pd.DataFrame(rows,columns=["system","rsid","variant","locus","genotype","evidence","meaning"])
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",46)
print(df.to_string(index=False)); df.to_csv(f"{W}/bloodgroups.tsv",sep="\t",index=False)
