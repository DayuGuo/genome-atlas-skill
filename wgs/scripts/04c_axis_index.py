#!/usr/bin/env python
"""Per-chromosome north/south position on the East Asian PC1 axis, WGS version (replaces the round-3 'north-south index')."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import subprocess, io, os, pandas as pd, numpy as np
P=str(P); W=f"{P}/wgs/04_ancestry"; O=f"{P}/wgs/18_ns"; os.makedirs(O,exist_ok=True)
wts=pd.read_csv(f"{W}/eas.pca.eigenvec.allele",sep="\t",dtype={"#CHROM":str})
wts.columns=["chrom","id","ref","alt","a1"]+[f"PC{i}" for i in range(1,11)]
af=pd.read_csv(f"{W}/eas.pca.afreq",sep="\t",dtype={"#CHROM":str}); af.columns=["chrom","id","ref","alt","alt_freq","n"]
kg=pd.read_csv(f"{W}/eas.proj.sscore",sep="\t").rename(columns={"#IID":"iid"})
me=pd.read_csv(f"{W}/eas.target.proj.sscore",sep="\t")
# per-chromosome scores: rerun --score restricted to each chromosome
rows=[]
for c in [str(i) for i in range(1,23)]:
    ids=wts[wts.chrom==c].id.drop_duplicates()
    idf=f"{O}/ids.{c}.txt"; ids.to_csv(idf,index=False,header=False)
    for who,pf,out in [("kg",f"{W}/kg.common",f"{O}/kg.{c}"),("target",f"{P}/wgs/02_complete/{SAMPLE}.1kg",f"{O}/target.{c}")]:
        subprocess.run([PLINK2,"--pfile",pf,"--extract",idf,"--read-freq",f"{W}/eas.pca.afreq",
                        "--score",f"{W}/eas.pca.eigenvec.allele","2","5","header-read","no-mean-imputation","variance-standardize",
                        "--score-col-nums","6","--out",out,"--threads","8","--memory","20000"],capture_output=True)
    k=pd.read_csv(f"{O}/kg.{c}.sscore",sep="\t"); d=pd.read_csv(f"{O}/target.{c}.sscore",sep="\t")
    k=k.rename(columns={"#IID":"iid"})
    chb=k[k.Population=="CHB"].PC1_AVG; chs=k[k.Population=="CHS"].PC1_AVG
    x=d.PC1_AVG.iloc[0]
    idx=(x-chs.mean())/(chb.mean()-chs.mean())
    rows.append(dict(chrom=c,n_snp=len(ids),pc1=x,index=round(idx,3),sd_from_chb=round((x-chb.mean())/chb.std(),2),sd_from_chs=round((x-chs.mean())/chs.std(),2)))
    print(rows[-1],flush=True)
    for f in [f"{O}/kg.{c}.sscore",f"{O}/target.{c}.sscore",idf]: os.remove(f)
df=pd.DataFrame(rows); df.to_csv(f"{O}/ns_index.tsv",sep="\t",index=False)
# genome-wide
chb=kg[kg.Population=="CHB"].PC1_AVG; chs=kg[kg.Population=="CHS"].PC1_AVG; x=me.PC1_AVG.iloc[0]
g=dict(index=round((x-chs.mean())/(chb.mean()-chs.mean()),3),sd_from_chb=round((x-chb.mean())/chb.std(),2),sd_from_chs=round((x-chs.mean())/chs.std(),2))
print("genome-wide:",g); pd.Series(g).to_csv(f"{O}/ns_genome.tsv",sep="\t",header=False)
print(f"per-chromosome index mean {df['index'].mean():.2f}, sd {df['index'].std():.2f}, range {df['index'].min():.2f}-{df['index'].max():.2f}")
