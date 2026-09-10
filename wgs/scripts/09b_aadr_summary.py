#!/usr/bin/env python
""f"Nearest present-day and ancient populations to {NAME_EN} in the Human Origins PCA space."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import pandas as pd, numpy as np
W=f"{P}/wgs/11_aadr"
s=pd.read_csv(f"{W}/proj.sscore",sep="\t").rename(columns={"#FID":"label","IID":"iid"})
meta=pd.read_csv(f"{W}/samples.tsv",sep="\t")
s=s.merge(meta[["iid","kind","date","call_rate"]],on="iid",how="left")
pcs=[f"PC{i}_AVG" for i in range(1,5)]
me=s[s.label==f"{NAME_EN}"][pcs].values[0]
s["d"]=np.sqrt(((s[pcs].values-me)**2).sum(1))
mod=s[(s.kind=="modern")]; anc=s[(s.kind=="ancient")&(s.call_rate>=0.5)]
def pop_tab(df,minn=2):
    g=df.groupby("label").agg(n=("d","size"),d=("d","mean"),date=("date","mean")).query("n>=@minn").sort_values("d")
    return g
print("== nearest present-day populations (PC1-4)"); print(pop_tab(mod).head(12).round(4).to_string())
print("\n== Han by province"); print(pop_tab(mod[mod.label.str.startswith("Han")],1).round(4).to_string())
print("\n== nearest ancient groups (>=2 individuals, call rate >=50%)"); print(pop_tab(anc).head(20).round(4).to_string())
print("\n== nearest ancient individuals"); print(anc.nsmallest(12,"d")[["label","iid","date","d"]].round(3).to_string(index=False))
pop_tab(mod).to_csv(f"{W}/near_modern.tsv",sep="\t"); pop_tab(anc).to_csv(f"{W}/near_ancient.tsv",sep="\t")
s.to_csv(f"{W}/proj_annotated.tsv",sep="\t",index=False)
