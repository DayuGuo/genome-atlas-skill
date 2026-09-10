#!/usr/bin/env python
"""Archaic (Neanderthal / Denisovan) introgressed segments carried by {NAME_EN}, using Browning 2018 Sprime segments called in CHB and CHS.
A segment counts as carried when {NAME_EN} carries the archaic allele at >=30% of its archaic-matching variants (>=10 variants typed)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import glob, io, subprocess, collections, pandas as pd, numpy as np
P=str(P); W=f"{P}/wgs/15_archaic"; A=SPRIME_DIR
sp=pd.concat([pd.read_csv(f,sep="\t",dtype={"CHROM":str}).assign(POP=f.split("/")[-1].split(".")[0]) for f in glob.glob(f"{A}/CH*.chr*.ND_match")],ignore_index=True)
sp=sp[(sp.REF.str.len()==1)&(sp.ALT.str.len()==1)]
sp["seg"]=sp.POP+":"+sp.CHROM+":"+sp.SEGMENT.astype(str)
print("Sprime variants:",len(sp)," segments:",sp.seg.nunique())
# {NAME_EN} genotypes at these positions from the complete set
reg=f"{W}/sites.tsv"; sp[["CHROM","POS"]].drop_duplicates().sort_values(["CHROM","POS"],key=lambda s: s.astype(int) if s.name=="CHROM" else s).to_csv(reg,sep="\t",header=False,index=False)
q=subprocess.run(["bcftools","query","-R",reg,"-f","%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n",f"{P}/wgs/02_complete/{SAMPLE}.1kg_sites.vcf.gz"],capture_output=True,text=True).stdout
g=pd.read_csv(io.StringIO(q),sep="\t",header=None,names=["CHROM","POS","REF","ALT","GT"],dtype={"CHROM":str})
m=sp.merge(g,on=["CHROM","POS","REF","ALT"],how="left")
m["alt_dose"]=m.GT.str.count("1")
m["arch_dose"]=np.where(m.ALLELE==1,m.alt_dose,2-m.alt_dose)
m=m[m.GT.notna()]
rows=[]
for seg,d in m.groupby("seg"):
    nm=d[d.NMATCH=="match"]; dm=d[d.DMATCH=="match"]
    typed=len(d); carried=(d.arch_dose>0).sum(); hom=(d.arch_dose==2).sum()
    n_match=len(nm); n_car=(nm.arch_dose>0).sum(); d_match=len(dm); d_car=(dm.arch_dose>0).sum()
    src="Neanderthal" if (d.NMATCH=="match").sum()>(d.DMATCH=="match").sum() else ("Denisovan" if (d.DMATCH=="match").sum()>0 else "ambiguous")
    rows.append(dict(seg=seg,pop=d.POP.iloc[0],chrom=d.CHROM.iloc[0],start=int(d.POS.min()),end=int(d.POS.max()),kb=round((d.POS.max()-d.POS.min())/1e3),n_typed=typed,
                     frac_carried=round(carried/typed,2),frac_hom=round(hom/typed,2),n_nmatch=n_match,frac_nmatch_carried=round(n_car/max(1,n_match),2),n_dmatch=d_match,frac_dmatch_carried=round(d_car/max(1,d_match),2),source=src))
r=pd.DataFrame(rows); r["carried"]=(r.n_typed>=10)&(r.frac_carried>=0.3); r["homozygous"]=r.carried&(r.frac_hom>=0.3)
r.to_csv(f"{W}/segments_all.tsv",sep="\t",index=False)
c=r[r.carried].copy()
# merge CHB/CHS duplicates of the same region
c=c.sort_values(["chrom","start"]); merged=[]
for row in c.itertuples():
    if merged and merged[-1]["chrom"]==row.chrom and row.start<=merged[-1]["end"]:
        merged[-1]["end"]=max(merged[-1]["end"],row.end); merged[-1]["source"]=merged[-1]["source"] if merged[-1]["source"]==row.source else "mixed"; merged[-1]["hom"]=merged[-1]["hom"] or row.homozygous
    else: merged.append(dict(chrom=row.chrom,start=row.start,end=row.end,source=row.source,hom=bool(row.homozygous)))
mg=pd.DataFrame(merged); mg["mb"]=(mg.end-mg.start)/1e6
mg.to_csv(f"{W}/segments_carried.bed",sep="\t",index=False,header=False)
tot=mg.mb.sum(); nea=mg[mg.source=="Neanderthal"].mb.sum(); den=mg[mg.source=="Denisovan"].mb.sum()
print(f"segments tested {len(r)}; carried {r.carried.sum()} ({len(mg)} after merging CHB/CHS); span carried: {tot:.1f} Mb (Neanderthal {nea:.1f}, Denisovan {den:.1f}, mixed/ambiguous {tot-nea-den:.1f}); homozygous {mg.hom.sum()}")
# population baseline: fraction of segments carried by an average individual ~ derive from Sprime? not available; report {NAME_EN}'s own numbers only.
print(mg.sort_values("mb",ascending=False).head(12).to_string(index=False))
open(f"{W}/summary.txt","w").write(f"segments_tested\t{len(r)}\ncarried\t{r.carried.sum()}\nmerged\t{len(mg)}\nspan_mb\t{tot:.1f}\nneanderthal_mb\t{nea:.1f}\ndenisovan_mb\t{den:.1f}\nhomozygous\t{int(mg.hom.sum())}\n")
