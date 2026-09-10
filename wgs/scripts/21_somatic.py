#!/usr/bin/env python
"""Blood-derived somatic signals: mtDNA copy number, mosaic loss of Y (mLOY), CHIP hotspot pileups and low-VAF coding variants in CHIP genes."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pathlib
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import subprocess, io, re, pandas as pd, numpy as np
P=str(P); W=f"{P}/wgs/13_somatic"; CRAM=f"{P}/wgs/00_input/{SAMPLE}.cram"; REF=FASTA
AUTO=28.88
out=[]
# 1. mtDNA copy number
mt=4850.22; out.append(("mtDNA_copies_per_cell", round(2*mt/AUTO,1), "2 x MT depth / autosomal depth"))
# 2. mLOY: median 1-kb bin depth in male-specific single-copy Y (X-degenerate) vs autosomes; expected 0.5
def bins(reg):
    t=subprocess.run(["tabix",f"{P}/wgs/01_qc/depth.regions.bed.gz",reg],capture_output=True,text=True).stdout
    return np.array([float(l.split("\t")[3]) for l in t.splitlines()])
y=np.concatenate([bins("Y:2781480-7000000"),bins("Y:14000000-16000000"),bins("Y:20000000-22000000")]); y=y[y>2]
x=np.concatenate([bins("X:20000000-40000000"),bins("X:80000000-100000000")]); x=x[x>2]
out.append(("Y_depth_ratio", round(float(np.median(y))/AUTO,3), "median mappable Y bin / autosomal mean; 0.5 = intact Y, <0.45 suggests mosaic loss"))
out.append(("X_depth_ratio", round(float(np.median(x))/AUTO,3), "median X bin / autosomal mean; 0.5 expected in a male"))
# 3. CHIP hotspots
PANEL = pathlib.Path(__file__).resolve().parents[1] / "panel"
HOT = {r.hotspot: (str(r.chrom), int(r.pos), r.ref, r.alt) for r in pd.read_csv(PANEL / "chip_hotspots.tsv", sep="\t").itertuples()}
rows=[]
for name,(c,p,r,a) in HOT.items():
    t=subprocess.run(f"samtools mpileup -f {REF} -r {c}:{p}-{p} -q 20 -Q 20 -d 1000 {CRAM} 2>/dev/null",shell=True,capture_output=True,text=True).stdout.split("\t")
    if len(t)<5: rows.append((name,c,p,0,0,0.0)); continue
    b=t[4]; b=re.sub(r"\^.","",b).replace("$","")
    ins=len(re.findall(r"\+\d+",b)); dele=len(re.findall(r"-\d+",b))
    b2=re.sub(r"[+-](\d+)[ACGTNacgtn]+",lambda m:"",b)
    dp=int(t[3]); alt=b2.upper().count(a) if len(a)==1 else ins
    rows.append((name,c,p,dp,alt,round(alt/dp,3) if dp else 0.0))
hot=pd.DataFrame(rows,columns=["hotspot","chrom","pos","depth","alt_reads","vaf"]); hot.to_csv(f"{W}/chip_hotspots.tsv",sep="\t",index=False)
# 4. low-VAF coding variants in CHIP genes (PASS + NO_PASS), VAF 0.05-0.35, depth>=20
genes = {r.gene: (str(r.chrom), int(r.start), int(r.end)) for r in pd.read_csv(PANEL / "chip_genes.tsv", sep="\t").itertuples()}
reg=",".join(f"{c}:{s}-{e}" for c,s,e in genes.values())
t=subprocess.run(["bcftools","query","-r",reg,"-f","%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/BCSQ\t[%GT\t%DP\t%AD]\n",f"{P}/wgs/00_input/{SAMPLE}.norm.vcf.gz"],capture_output=True,text=True).stdout
csq=subprocess.run(["bcftools","query","-r",reg,"-f","%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\n",f"{P}/wgs/05_clinvar/{SAMPLE}.pass.annot.vcf.gz"],capture_output=True,text=True).stdout
cmap={tuple(l.split("\t")[:4]):l.split("\t")[4] for l in csq.splitlines()}
lv=[]
for l in t.splitlines():
    c,p,r,a,flt,_,gt,dp,ad=l.split("\t")
    try: ref,alt=map(int,ad.split(",")[:2]); dp=int(dp)
    except: continue
    if dp>=20 and 0.05<=alt/dp<=0.35:
        g=[k for k,(cc,s,e) in genes.items() if cc==c and s<=int(p)<=e][0]
        lv.append((g,c,p,r,a,flt,dp,alt,round(alt/dp,3),cmap.get((c,p,r,a),"(NO_PASS: unannotated)")[:80]))
low=pd.DataFrame(lv,columns=["gene","chrom","pos","ref","alt","filter","depth","alt_reads","vaf","consequence"]); low.to_csv(f"{W}/chip_lowvaf.tsv",sep="\t",index=False)
pd.DataFrame(out,columns=["metric","value","note"]).to_csv(f"{W}/summary.tsv",sep="\t",index=False)
pd.set_option("display.width",220); print(pd.DataFrame(out,columns=["metric","value","note"]).to_string(index=False)); print(); print(hot.to_string(index=False)); print(); print("low-VAF (5-35%) variants in CHIP genes:",len(low)); print(low[low.consequence.str.contains("missense|stop|frameshift|splice")].to_string(index=False) if len(low) else "")
