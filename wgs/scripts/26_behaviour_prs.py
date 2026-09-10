#!/usr/bin/env python
"""Behaviour and psychiatric polygenic scores, split into East-Asian-derived and European-derived sets so the
transferability gap is visible. Same machinery as the disease PRS: EAS MAF>=0.01, 504 1000G East Asians as reference."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pathlib
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import gzip, io, json, os, pathlib, subprocess, urllib.request
import numpy as np, pandas as pd
ROOT = P; D=PGS; FV=ROOT/"data/ref/prs"; W=ROOT/"wgs/20_behaviour"
W.mkdir(parents=True,exist_ok=True); P=str(ROOT/"tools/plink2"); MAF=0.01

_sc = pd.read_csv(pathlib.Path(__file__).resolve().parents[1] / "panel" / "pgs_scores.tsv", sep="\t")
_sc = _sc[_sc.group == "behaviour"]
EAS_SET = {r.pgs: (r.trait_zh, r.trait_en, r.source) for r in _sc[_sc.panel == "EAS"].itertuples()}
EUR_SET = {r.pgs: (r.trait_zh, r.trait_en, r.source) for r in _sc[_sc.panel != "EAS"].itertuples()}
ALL={**EAS_SET,**EUR_SET}
def fetch(pid):
    f=D/f"{pid}_hmPOS_GRCh37.txt.gz"
    if not f.exists():
        url=f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pid}/ScoringFiles/Harmonized/{pid}_hmPOS_GRCh37.txt.gz"
        print("downloading",pid,flush=True); urllib.request.urlretrieve(url,f)
    return f
idx=pd.read_csv(FV/"kg_all.pvar",sep="\t",comment="#",header=None,names=["chrom","pos","id","ref","alt"],dtype={"chrom":str})
fr=pd.read_csv(FV/"kg_all.afreq",sep="\t",usecols=["ID","ALT_FREQS"]); fr["maf"]=np.minimum(fr.ALT_FREQS,1-fr.ALT_FREQS)
common=set(fr[fr.maf>=MAF].ID)
dayu_ids=set(l.split("\t")[2] for l in open(ROOT/f"wgs/07_prs/{SAMPLE}_pos.pvar") if not l.startswith("#"))
usable=common&dayu_ids
print(f"usable variants: {len(usable):,}")
def run(*a):
    r=subprocess.run([P,*map(str,a),"--threads","16","--memory","30000"],capture_output=True,text=True)
    if r.returncode: raise SystemExit(r.stdout[-2000:]+r.stderr[-2000:])
rows=[]
for pid,(zh,en,src) in ALL.items():
    f=fetch(pid)
    s=pd.read_csv(f,sep="\t",comment="#",dtype={"hm_chr":str},low_memory=False).dropna(subset=["hm_chr","hm_pos"])
    s["hm_pos"]=s.hm_pos.astype(int)
    if "other_allele" not in s: s["other_allele"]=np.nan
    n_total=len(s)
    mm=s.merge(idx,left_on=["hm_chr","hm_pos"],right_on=["chrom","pos"])
    ok=((mm.effect_allele==mm.ref)&((mm.other_allele==mm.alt)|mm.other_allele.isna()))|((mm.effect_allele==mm.alt)&((mm.other_allele==mm.ref)|mm.other_allele.isna()))
    mm=mm[ok].copy(); mm["sid"]=mm.chrom+":"+mm.pos.astype(str)
    keep=mm[mm.sid.isin(usable)]
    sf=W/f"{pid}.score"; keep[["sid","effect_allele","effect_weight"]].to_csv(sf,sep="\t",index=False,header=False)
    for who,pf in [("eas",FV/"kg_all"),("target",ROOT/f"wgs/07_prs/{SAMPLE}_pos")]:
        run("--pfile",pf,"--score",sf,"1","2","3","cols=+scoresums","no-mean-imputation","--out",W/f"{pid}.{who}")
    eas=pd.read_csv(W/f"{pid}.eas.sscore",sep="\t"); d=pd.read_csv(W/f"{pid}.target.sscore",sep="\t")
    x=d.SCORE1_SUM.iloc[0]; han=eas[eas.Population.isin(["CHB","CHS"])].SCORE1_SUM
    z=(x-eas.SCORE1_SUM.mean())/eas.SCORE1_SUM.std(); pct=(eas.SCORE1_SUM<x).mean()*100
    rows.append(dict(pgs=pid,zh=zh,en=en,source=src,panel="EAS" if pid in EAS_SET else "EUR",
                     variants=n_total,used=len(keep),coverage=round(100*len(keep)/n_total,1),
                     z=round(float(z),2),pct_EAS=round(float(pct),1),pct_Han=round(float((han<x).mean()*100),1)))
    print(f"{pid} {zh:<14} {rows[-1]['panel']}  used {len(keep):>8,}/{n_total:<9,} ({rows[-1]['coverage']:5.1f}%)  pct={pct:5.1f} z={z:+.2f}",flush=True)
df=pd.DataFrame(rows); df.to_csv(W/"behaviour_prs.tsv",sep="\t",index=False)
pd.set_option("display.width",240)
print("\n== 东亚人群训练的评分"); print(df[df.panel=="EAS"][["zh","en","used","coverage","pct_EAS","z"]].to_string(index=False))
print("\n== 欧洲人群训练的评分（迁移性差，仅作对照）"); print(df[df.panel=="EUR"][["zh","en","used","coverage","pct_EAS","z"]].to_string(index=False))
