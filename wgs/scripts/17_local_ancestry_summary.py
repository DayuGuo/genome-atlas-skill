#!/usr/bin/env python
"""Genome-wide and per-segment local ancestry from FLARE output."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import glob, gzip, io, os, subprocess, pandas as pd, numpy as np
P=str(P); W=f"{P}/wgs/12_localanc"
ANC=["European","SouthEA","SouthAsian","NorthEA"]
segs=[]; tot=np.zeros(len(ANC)); n_mark=0
for f in sorted(glob.glob(f"{W}/la.*.anc.vcf.gz"),key=lambda x:int(x.split("la.")[1].split(".")[0])):
    c=f.split("la.")[1].split(".")[0]
    q=subprocess.run(["bcftools","query","-f","%POS\t[%AN1\t%AN2]\n",f],capture_output=True,text=True).stdout
    d=pd.read_csv(io.StringIO(q),sep="\t",header=None,names=["pos","a1","a2"])
    if not len(d): continue
    n_mark+=len(d)
    for col in ("a1","a2"):
        v=d[col].values
        for i in range(len(ANC)): tot[i]+=(v==i).sum()
        # segments
        chg=np.flatnonzero(np.diff(v))+1; bounds=np.concatenate([[0],chg,[len(v)]])
        for s,e in zip(bounds[:-1],bounds[1:]):
            segs.append((c,col,int(d.pos.iloc[s]),int(d.pos.iloc[e-1]),ANC[v[s]],e-s))
    print(f"chr{c}: {len(d)} markers",flush=True)
sg=pd.DataFrame(segs,columns=["chrom","hap","start","end","anc","n_markers"])
sg["mb"]=(sg.end-sg.start)/1e6
sg.to_csv(f"{W}/segments.tsv",sep="\t",index=False)
frac=tot/tot.sum()
print("\ngenome-wide ancestry fractions (both haplotypes):")
for a,f_ in zip(ANC,frac): print(f"  {a:12s} {f_*100:6.2f}%")
big=sg[(sg.anc=="SouthEA")&(sg.mb>=2)].sort_values("mb",ascending=False)
print(f"\nSouthEA segments >=2 Mb: {len(big)}; longest {big.mb.max():.1f} Mb" if len(big) else "\nno SouthEA segment >=2 Mb")
print(big.head(10).to_string(index=False))
sw=sg.groupby("anc").agg(n=("mb","size"),mb=("mb","sum"),median_mb=("mb","median"))
print("\nsegment counts:"); print(sw.round(2).to_string())
pd.Series(dict(zip(ANC,frac))).to_csv(f"{W}/global.tsv",sep="\t",header=False)
