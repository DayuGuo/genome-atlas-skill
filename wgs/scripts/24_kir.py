#!/usr/bin/env python
"""KIR gene content from T1K plus the HLA class I ligand groups (C1/C2, Bw4/Bw6) that pair with them."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import pandas as pd
P=str(P); W=f"{P}/wgs/16_panels"
k=pd.read_csv(f"{P}/wgs/06_pgx/t1k/kir_genotype.tsv",sep="\t",header=None,
              names=["gene","n","a1","ab1","q1","a2","ab2","q2","other"])
k["present"]=(k.n>0)&(k.ab1>=5)
k["call"]=[(f"{r.a1}" + (f" / {r.a2}" if r.n>1 and r.ab2>=5 else "")) if r.present else "absent" for r in k.itertuples()]
k[["gene","present","call","ab1","ab2"]].to_csv(f"{W}/kir.tsv",sep="\t",index=False)
pres=set(k[k.present].gene)
Ahap={"KIR3DL3","KIR2DL3","KIR2DP1","KIR2DL1","KIR3DP1","KIR2DL4","KIR3DL1","KIR2DS4","KIR3DL2"}
Bonly={"KIR2DL2","KIR2DL5A","KIR2DL5B","KIR2DS1","KIR2DS2","KIR2DS3","KIR2DS5","KIR3DS1"}
hap="AA (two A haplotypes)" if not (pres & Bonly) else "Bx (at least one B haplotype)"
print("KIR genes present:", ", ".join(sorted(pres)))
print("B-haplotype-specific genes present:", ", ".join(sorted(pres & Bonly)) or "none")
print("KIR haplotype:", hap)
# HLA ligands
hla={"A":["02:07","24:02"],"B":["13:01","54:01"],"C":["03:04","01:02"]}
C1={"01","03","07","08","12","14","16"}; C2={"02","04","05","06","15","17","18"}
cg=[("C1" if a.split(":")[0] in C1 else "C2" if a.split(":")[0] in C2 else "?") for a in hla["C"]]
BW4={"13:01","27:05","44:02","44:03","51:01","52:01","53:01","57:01","58:01","38:01","37:01","47:01","49:01","57:03","59:01","63:01","77:01"}
bw=["Bw4" if a in BW4 else "Bw6" for a in hla["B"]]
print(f"HLA-C ligands: C*{hla['C'][0]} = {cg[0]}, C*{hla['C'][1]} = {cg[1]}  ->  {'C1/C1' if cg==['C1','C1'] else '/'.join(cg)}")
print(f"HLA-B epitopes: B*{hla['B'][0]} = {bw[0]}, B*{hla['B'][1]} = {bw[1]}  ->  {'/'.join(bw)}")
print("KIR2DL3 (inhibitory, C1-specific) has its ligand; KIR3DL1 has one Bw4 ligand (B*13:01); no activating KIR beyond 2DS4.")
open(f"{W}/kir_summary.txt","w").write(f"haplotype\t{hap}\npresent\t{','.join(sorted(pres))}\nC_ligands\t{'/'.join(cg)}\nB_epitopes\t{'/'.join(bw)}\n")
