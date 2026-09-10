#!/usr/bin/env python
"""Extract a chosen set of AADR Human Origins individuals (TGENO packed format) plus {NAME_EN}'s genotypes
at the same SNPs, and write a plink1 .bed/.bim/.fam set for PCA."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import os, sys, io, subprocess, numpy as np, pandas as pd
P=str(P); A=str(pathlib.Path(AADR).parent); W=f"{P}/wgs/11_aadr"; os.makedirs(W,exist_ok=True)
PREF=f"{A}/v66.p1_HO.aadr.patch.PUB"
ind=pd.read_csv(f"{PREF}.ind",sep=r"\s+",header=None,names=["iid","sex","pop"])
snp=pd.read_csv(f"{PREF}.snp",sep=r"\s+",header=None,names=["rsid","chrom","cm","pos","a1","a2"],dtype={"chrom":str})
n_ind,n_snp=len(ind),len(snp); rlen=(n_snp+3)//4
assert 48+rlen*n_ind==os.path.getsize(f"{PREF}.geno"), "layout mismatch"

MODERN=["Han","Japanese","Korean","She","Miao","Tujia","Dai","Yi","Naxi","Lahu","Tu","Xibo","Oroqen","Hezhen","Daur","Mongola","Tibetan","Kinh_Vietnamese","Thai","Cambodian","Ami","Atayal","Zhuang","Dong","Mulam","Maonan","Gelao","Li","Qiang","Yugur","Salar","Bonan","Dongxiang","Ulchi","Nanai","Nivh","Evenk_Transbaikal","Buryat","Mongol","Uygur","Hmong","Burmese","Lao","Khmer","Malay","Igorot","CHB","CHS","CDX","Sherpa"]
ANC_KEY=["China_AmurRiverBasin","China_YellowRiver","China_Baligang","China_Henan","China_Shandong","China_Qingdao_BeiQian","China_Jinan","China_Jining","China_Shanxi","China_Shaanxi","China_InnerMongolia","China_LN_Xiaoheyan","China_Fujian","China_Taiwan","Taiwan_IA","Taiwan_EN","China_Guangxi","China_Sichuan","China_Qinghai","China_Tibet_","China_TianyuanCave","China_EBA","China_MBA","China_LBA","China_IA","China_MLBA","Mongolia_N","Mongolia_East_N","Mongolia_MLBA","Mongolia_XiongnuPeriod","Vietnam_N","Vietnam_BA","Laos_Hoabinhian","Thailand_BA","Japan_Jomon","Japan_","Korea_","Russia_DevilsCave","Russia_Boisman","Russia_AmurRiver","Russia_Shamanka","Russia_UstIda"]
anno=pd.read_csv(f"{PREF.replace('.patch.PUB','.PUB')}.anno",sep="\t",dtype=str,low_memory=False)
anno.columns=[c.strip() for c in anno.columns]
gcol=[c for c in anno.columns if c.startswith("Group ID")][0]; dcol=[c for c in anno.columns if c.startswith("Date mean in BP")][0]
locc=[c for c in anno.columns if c.startswith("Locality")][0]; idc=anno.columns[0]
anno["date"]=pd.to_numeric(anno[dcol],errors="coerce").fillna(0)
meta=anno[[idc,gcol,locc,"date"]].rename(columns={idc:"iid",gcol:"pop",locc:"loc"})
meta["iid"]=meta.iid.str.strip()
ind=ind.merge(meta.drop_duplicates("iid"),on="iid",how="left",suffixes=("","_a"))
ind["pop"]=ind.pop_a.fillna(ind["pop"])
bad=ind["pop"].str.contains("Ignore|QCremove|DontUse|_dup|Discovery|_o$|-o$|outlier|contam|_lc$|_1d|_2d",case=False,na=False)
is_mod=ind["pop"].isin(MODERN)&(ind.date<=50)
is_anc=ind["pop"].str.startswith(tuple(ANC_KEY),na=False)&(ind.date>=500)
keep=ind[(is_mod|is_anc)&~bad].copy()
keep["kind"]=np.where(keep["pop"].isin(MODERN),"modern","ancient")
# province label for Han
prov=keep["loc"].astype(str)
keep["label"]=np.where((keep["pop"]=="Han")&~prov.str.contains("modern"),"Han_"+prov.str.replace(r'[",]','',regex=True).str.strip(),keep["pop"])
print(f"keeping {len(keep)} individuals: modern {(keep.kind=='modern').sum()}, ancient {(keep.kind=='ancient').sum()}",file=sys.stderr)
print(keep[keep.kind=="modern"].label.value_counts().head(20).to_dict(),file=sys.stderr)

# --- SNP selection: autosomal, biallelic ACGT, present in {NAME_EN}'s complete set
snp["idx"]=np.arange(n_snp)
sel=snp[(snp.chrom.isin([str(i) for i in range(1,23)]))&snp.a1.isin(list("ACGT"))&snp.a2.isin(list("ACGT"))].copy()
sel[["chrom","pos"]].to_csv(f"{W}/ho_sites.tsv",sep="\t",header=False,index=False)
q=subprocess.run(["bcftools","query","-R",f"{W}/ho_sites.tsv","-f","%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n",f"{P}/wgs/02_complete/{SAMPLE}.1kg_sites.vcf.gz"],capture_output=True,text=True).stdout
dg=pd.read_csv(io.StringIO(q),sep="\t",header=None,names=["chrom","pos","ref","alt","gt"],dtype={"chrom":str}).drop_duplicates(["chrom","pos"])
sel=sel.merge(dg,on=["chrom","pos"],how="inner")
ok=((sel.a1==sel.ref)&(sel.a2==sel.alt))|((sel.a1==sel.alt)&(sel.a2==sel.ref))
sel=sel[ok].copy()
sel["dayu_a1"]=np.where(sel.a1==sel.alt,sel["gt"].str.count("1"),2-sel["gt"].str.count("1"))  # count of a1 (EIGENSTRAT convention)
print(f"SNPs usable: {len(sel)} of {n_snp}",file=sys.stderr)

# --- read packed rows for the kept individuals
G=np.zeros((len(keep),len(sel)),dtype=np.int8)
cols=sel.idx.values
with open(f"{PREF}.geno","rb") as f:
    for k,(row_i) in enumerate(keep.index.values):
        f.seek(48+rlen*row_i); buf=np.frombuffer(f.read(rlen),dtype=np.uint8)
        bits=np.unpackbits(buf).reshape(-1,2); vals=(bits[:,0]*2+bits[:,1])[:n_snp]
        v=vals[cols].astype(np.int8); v[v==3]=-1  # 3 = missing
        G[k]=v
        if k%200==0: print(f"  {k}/{len(keep)}",file=sys.stderr,flush=True)
miss=(G<0).mean(axis=1)
keep["call_rate"]=1-miss
keep2=keep[keep.call_rate>=0.3].copy(); G=G[keep.call_rate.values>=0.3]
print(f"after call-rate filter: {len(keep2)} individuals",file=sys.stderr)

# --- write plink1 bed (SNP-major): 00=hom a1a1? plink codes: 00 hom A1, 01 missing, 10 het, 11 hom A2
# We set A1 = eigen a1 (counted allele). geno value g = count of a1. plink 2-bit: g==2 -> 00, g==1 -> 10, g==0 -> 11, missing -> 01
allg=np.vstack([G,sel.dayu_a1.values.astype(np.int8)[None,:]])
labels=list(keep2.label)+[f"{NAME_EN}"]; kinds=list(keep2.kind)+["target"]; iids=list(keep2.iid)+[f"{NAME_EN}"]
code=np.select([allg==2,allg==1,allg==0],[0b00,0b10,0b11],default=0b01).astype(np.uint8)
n=len(labels); pad=(4-n%4)%4
with open(f"{W}/aadr.bed","wb") as f:
    f.write(bytes([0x6c,0x1b,0x01]))
    for j in range(allg.shape[1]):
        col=np.concatenate([code[:,j],np.zeros(pad,dtype=np.uint8)]).reshape(-1,4)
        f.write((col[:,0]|(col[:,1]<<2)|(col[:,2]<<4)|(col[:,3]<<6)).astype(np.uint8).tobytes())
pd.DataFrame({"chrom":sel.chrom,"rsid":sel.rsid,"cm":sel.cm,"pos":sel.pos,"a1":sel.a1,"a2":sel.a2}).to_csv(f"{W}/aadr.bim",sep="\t",header=False,index=False)
pd.DataFrame({"fid":labels,"iid":iids,"pat":0,"mat":0,"sex":0,"phe":-9}).to_csv(f"{W}/aadr.fam",sep=" ",header=False,index=False)
pd.DataFrame({"iid":iids,"label":labels,"kind":kinds,"date":list(keep2.date)+[0],"call_rate":list(keep2.call_rate)+[1.0]}).to_csv(f"{W}/samples.tsv",sep="\t",index=False)
print("wrote",W,"n_ind",n,"n_snp",allg.shape[1])
