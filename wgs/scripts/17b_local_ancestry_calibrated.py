#!/usr/bin/env python
"""Genome-wide north/south ancestry proportions for {NAME_EN}, expressed against held-out CHB and CHS individuals
scored with the same panels (FLARE posterior global ancestry, weighted by chromosome length)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import glob, gzip, pandas as pd, numpy as np, os
P=str(P); W=f"{P}/wgs/12_localanc"
LEN={str(i):l for i,l in zip(range(1,23),[249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566])}
rows=[]
for f in glob.glob(f"{W}/la.*.global.anc.gz"):
    c=os.path.basename(f).split(".")[1]
    d=pd.read_csv(f,sep="\t"); d["chrom"]=c; rows.append(d)
dy=pd.concat(rows); dy["w"]=dy.chrom.map(LEN)
anc=["NorthEA","SouthEA","European","SouthAsian"]
gw={a:np.average(dy[a],weights=dy.w) for a in anc}
print(f"{NAME_EN}, length-weighted over 22 autosomes:")
for a in anc: print(f"  {a:12s} {gw[a]*100:5.2f}%")
ps=pd.read_csv(KG_PFILE+".psam",sep="\t").rename(columns={"#IID":"SAMPLE"})
cal=[]
for f in glob.glob(f"{W}/calib.*.global.anc.gz"):
    c=os.path.basename(f).split(".")[1]
    d=pd.read_csv(f,sep="\t"); d["chrom"]=c; cal.append(d)
if cal:
    cd=pd.concat(cal).merge(ps[["SAMPLE","Population"]],on="SAMPLE"); cd["w"]=cd.chrom.map(LEN)
    print(f"\ncalibration on chromosomes {sorted(set(cd.chrom))} (held-out reference individuals):")
    out=[]
    for pop,g in cd.groupby("Population"):
        per=g.groupby("SAMPLE").apply(lambda x: pd.Series({a:np.average(x[a],weights=x.w) for a in anc}),include_groups=False)
        out.append((pop,len(per),per.NorthEA.mean(),per.NorthEA.std(),per.SouthEA.mean(),per.SouthEA.std()))
        print(f"  {pop}: n={len(per)}  NorthEA {per.NorthEA.mean()*100:.1f}% (sd {per.NorthEA.std()*100:.1f})  SouthEA {per.SouthEA.mean()*100:.1f}% (sd {per.SouthEA.std()*100:.1f})")
    # {NAME_EN} restricted to the same chromosomes for a fair comparison
    same=dy[dy.chrom.isin(set(cd.chrom))]
    dn=np.average(same.NorthEA,weights=same.w); ds=np.average(same.SouthEA,weights=same.w)
    print(f"  {NAME_EN} (same chromosomes): NorthEA {dn*100:.1f}%  SouthEA {ds*100:.1f}%")
    for pop,n,m,s,ms,ss in out:
        print(f"    vs {pop}: {NAME_EN} is {(dn-m)/s:+.1f} sd on NorthEA")
    pd.DataFrame(out,columns=["pop","n","north_mean","north_sd","south_mean","south_sd"]).to_csv(f"{W}/calibration.tsv",sep="\t",index=False)
pd.Series(gw).to_csv(f"{W}/dayu_global.tsv",sep="\t",header=False)
dy.to_csv(f"{W}/per_chrom.tsv",sep="\t",index=False)
print(f"\nper-chromosome SouthEA share ({NAME_EN}):")
print(dy.sort_values("chrom",key=lambda s:s.astype(int))[["chrom","NorthEA","SouthEA"]].round(3).to_string(index=False))
