#!/usr/bin/env python
""f"96-context germline SNV spectrum (PASS, autosomal), pyrimidine-centred, plus the same spectrum for 1000G EAS-common vs {NAME_EN}-private SNVs."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import subprocess, collections, pysam, pandas as pd
P=str(P); fa=pysam.FastaFile(FASTA)
comp={"A":"T","C":"G","G":"C","T":"A"}
q=subprocess.run(["bcftools","query","-i",'TYPE="snp"',"-r",",".join(map(str,range(1,23))),"-f","%CHROM\t%POS\t%REF\t%ALT\t%INFO/EAS_AF\t%INFO/ALL_AF\n",f"{P}/wgs/05_clinvar/{SAMPLE}.pass.annot.vcf.gz"],capture_output=True,text=True).stdout
cnt=collections.Counter(); cnt_priv=collections.Counter(); n=0
for l in q.splitlines():
    c,p,r,a,eas,al=l.split("\t"); p=int(p); ctx=fa.fetch(c,p-2,p+1).upper()
    if len(ctx)!=3 or "N" in ctx: continue
    if r in "AG": r,a,ctx=comp[r],comp[a],"".join(comp[b] for b in reversed(ctx))
    k=f"{ctx[0]}[{r}>{a}]{ctx[2]}"; cnt[k]+=1; n+=1
    if eas=="." and al==".": cnt_priv[k]+=1
subs=["C>A","C>G","C>T","T>A","T>C","T>G"]; rows=[]
for s in subs:
    for l5 in "ACGT":
        for r3 in "ACGT":
            k=f"{l5}[{s}]{r3}"; rows.append((s,k,cnt[k],cnt_priv[k]))
df=pd.DataFrame(rows,columns=["sub","context","n","n_private"]); df["frac"]=df.n/df.n.sum(); df["frac_private"]=df.n_private/max(1,df.n_private.sum())
df.to_csv(f"{P}/wgs/17_mutspec/spectrum96.tsv",sep="\t",index=False)
print("SNVs:",n," private (absent from 1000G):",df.n_private.sum()); print(df.groupby("sub")[["n","n_private"]].sum().assign(frac=lambda d:(d.n/d.n.sum()).round(3),frac_priv=lambda d:(d.n_private/d.n_private.sum()).round(3)))
cpg=df[df.context.str.contains(r"\[C>T\]G")].n.sum(); print("CpG C>T share of all:",round(cpg/n,3))
