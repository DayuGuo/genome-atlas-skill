#!/usr/bin/env python
"""Which genes the carried archaic segments cover, and how that compares with a matched random expectation."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pathlib
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import subprocess, io, os, random, collections
import pandas as pd, numpy as np
P=str(P); W=f"{P}/wgs/15_archaic"
GENES=f"{W}/genes.bed"
if not os.path.exists(GENES):
    subprocess.run(f"zcat {P}/data/ref/annot/Homo_sapiens.GRCh37.87.gff3.gz | "
                   "awk -F'\\t' '$3==\"gene\" && $9~/biotype=protein_coding/ {match($9,/Name=[^;]+/); "
                   "n=substr($9,RSTART+5,RLENGTH-5); print $1\"\\t\"$4\"\\t\"$5\"\\t\"n}' | sort -k1,1 -k2,2n > "+GENES,
                   shell=True,check=True)
g=pd.read_csv(GENES,sep="\t",header=None,names=["chrom","start","end","gene"],dtype={"chrom":str})
g=g[g.chrom.isin([str(i) for i in range(1,23)])]
seg=pd.read_csv(f"{W}/segments_carried.bed",sep="\t",header=None,names=["chrom","start","end","src","hom","mb"],dtype={"chrom":str})
seg["id"]=seg.chrom+":"+seg.start.astype(str)+"-"+seg.end.astype(str)
# genes per segment
rows=[]
for s in seg.itertuples():
    sub=g[(g.chrom==s.chrom)&(g.end>=s.start)&(g.start<=s.end)]
    rows.append((s.id,s.chrom,s.start,s.end,round(s.mb,2),s.src,bool(s.hom),len(sub),",".join(sorted(sub.gene))))
d=pd.DataFrame(rows,columns=["seg","chrom","start","end","mb","src","hom","n_genes","genes"])
d.to_csv(f"{W}/segments_genes.tsv",sep="\t",index=False)
allg=sorted(set(x for gs in d.genes if gs for x in gs.split(",")))
print(f"片段 {len(d)} 段，覆盖蛋白编码基因 {len(allg)} 个，占全部 {len(g)} 个的 {len(allg)/len(g)*100:.1f}%")
# gene families that stand out
_fam = pd.read_csv(pathlib.Path(__file__).resolve().parents[1] / "panel" / "gene_families.tsv", sep="\t")
FAM = {r.family_zh: r.regex for r in _fam.itertuples()}
ser=pd.Series(allg)
fam=[]
for name,pat in FAM.items():
    hit=sorted(ser[ser.str.match(pat)]); tot=sorted(g.gene[g.gene.str.match(pat)].unique())
    if hit: fam.append((name,len(hit),len(tot),round(len(hit)/len(tot)*100,1),", ".join(hit[:12])+(" …" if len(hit)>12 else "")))
EN = {r.family_zh: r.family_en for r in _fam.itertuples()}
fdf=pd.DataFrame(fam,columns=["family","carried","total","pct","examples"]).sort_values("pct",ascending=False)
fdf["family_en"]=fdf.family.map(EN).fillna(fdf.family)
fdf.to_csv(f"{W}/gene_families.tsv",sep="\t",index=False)
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",70)
print("\n== 基因家族的覆盖比例（该家族有多少比例的基因落在渗入片段里）")
print(fdf.to_string(index=False))
# background expectation: shuffle segments within chromosomes 200x, count genes hit
CHRL={str(i):l for i,l in zip(range(1,23),[249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566])}
random.seed(7)
byc={c:sub.sort_values("start")[["start","end","gene"]].values for c,sub in g.groupby("chrom")}
def count_hits(segments):
    n=set()
    for c,s,e in segments:
        arr=byc.get(c)
        if arr is None: continue
        for gs,ge,gn in arr:
            if gs>e: break
            if ge>=s: n.add(gn)
    return len(n)
obs=len(allg)
null=[]
segs=[(r.chrom,r.start,r.end) for r in seg.itertuples()]
for _ in range(200):
    sh=[]
    for c,s,e in segs:
        L=e-s; st=random.randint(0,max(1,CHRL[c]-L)); sh.append((c,st,st+L))
    null.append(count_hits(sh))
null=np.array(null)
print(f"\n观察到覆盖 {obs} 个基因；随机打乱片段位置 200 次的期望 {null.mean():.0f} ± {null.std():.0f}，"
      f"z = {(obs-null.mean())/null.std():+.2f}")
open(f"{W}/gene_summary.txt","w").write(f"n_genes\t{obs}\nnull_mean\t{null.mean():.1f}\nnull_sd\t{null.std():.1f}\nz\t{(obs-null.mean())/null.std():.2f}\n")
print("\n== 覆盖基因最多的 10 段")
print(d.nlargest(10,"n_genes")[["seg","mb","src","hom","n_genes"]].to_string(index=False))
