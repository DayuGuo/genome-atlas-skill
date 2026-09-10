#!/usr/bin/env python3
"""Polygenic risk scores from PGS Catalog (harmonized GRCh37), scored for the sample and for all
2504 1000 Genomes samples on the same variant subset, so the sample's score can be placed as a
percentile within the East Asian reference (EAS, n=504) and overall.
Only variants present on the array AND biallelic in 1000G are used; coverage is reported.
"""
import gzip, io, pathlib, subprocess, sys, requests, numpy as np, pandas as pd
from config import *
ROOT = pathlib.Path(__file__).resolve().parents[1]; A = RESULTS/"ancestry"; W = RESULTS/"prs"; W.mkdir(exist_ok=True)
D = DATA/"pgs"; P = pathlib.Path(PLINK2)

import os
_tab = pd.read_csv(PANEL / "pgs_scores.tsv", sep="\t")
GROUP = os.environ.get("PRS_GROUP", "disease")
scores = {r.pgs_id: (r.label_zh, r.source_note) for r in _tab[_tab.group == GROUP].itertuples()}
outname = "prs_summary.tsv" if GROUP == "disease" else "prs_biomarkers.tsv"
kg = pd.read_csv(A/"kg.common.pvar", sep="\t", comment="#", header=None, usecols=[0,1,2,3,4],
                 names=["chrom","pos","id","ref","alt"], dtype={"chrom":str})
def fetch(pid):
    f = D/f"{pid}_hmPOS_GRCh37.txt.gz"
    if not f.exists():
        url = f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pid}/ScoringFiles/Harmonized/{pid}_hmPOS_GRCh37.txt.gz"
        r = requests.get(url, timeout=600); r.raise_for_status(); f.write_bytes(r.content)
    return f
rows = []
for pid, (label, note) in scores.items():
    f = fetch(pid)
    with gzip.open(f, "rt") as fh:
        hdr = [l for l in fh if l.startswith("#")]
    s = pd.read_csv(f, sep="\t", comment="#", dtype={"hm_chr":str, "chr_name":str}, low_memory=False)
    n_total = len(s)
    s = s.dropna(subset=["hm_chr","hm_pos"]); s["hm_pos"] = s.hm_pos.astype(int)
    if "other_allele" not in s: s["other_allele"] = np.nan
    m = s.merge(kg.reset_index(drop=True), left_on=["hm_chr","hm_pos"], right_on=["chrom","pos"])
    ok = ((m.effect_allele==m.ref)&((m.other_allele==m.alt)|m.other_allele.isna())) | ((m.effect_allele==m.alt)&((m.other_allele==m.ref)|m.other_allele.isna()))
    m = m[ok].drop_duplicates("id")
    sf = W/f"{pid}.score"; m[["id","effect_allele","effect_weight"]].to_csv(sf, sep="\t", index=False, header=False)
    for tag, pf in [("kg", A/"kg.common"), ("sample", A/"sample")]:
        subprocess.run([P,"--pfile",pf,"--score",sf,"1","2","3","cols=+scoresums","no-mean-imputation","--out",W/f"{pid}.{tag}","--threads", THREADS],
                       check=True, capture_output=True)
    k = pd.read_csv(W/f"{pid}.kg.sscore", sep="\t").rename(columns={"#IID":"IID"})
    d = pd.read_csv(W/f"{pid}.sample.sscore", sep="\t")
    x = d.SCORE1_SUM.iloc[0]
    eas = k[k.SuperPop==REF_SUPERPOP].SCORE1_SUM; han = k[k.Population.isin(SUBPOPS)].SCORE1_SUM
    pct_eas = (eas < x).mean()*100; pct_han = (han < x).mean()*100; pct_all = (k.SCORE1_SUM < x).mean()*100
    z_eas = (x - eas.mean())/eas.std()
    rows.append(dict(pgs=pid, trait=label, source=note, variants_in_score=n_total, variants_used=len(m),
                     coverage_pct=round(100*len(m)/n_total,1), sample_score=x, z_vs_EAS=round(z_eas,2),
                     pct_EAS=round(pct_eas,1), pct_Han=round(pct_han,1), pct_all1000G=round(pct_all,1)))
    print(f"{pid} {label:10s} used {len(m):>7d}/{n_total:<8d} ({100*len(m)/n_total:4.1f}%)  z_EAS={z_eas:+.2f}  pct_EAS={pct_eas:5.1f}  pct_Han={pct_han:5.1f}", flush=True)
    k[["IID","SuperPop","Population","SCORE1_SUM"]].to_csv(W/f"{pid}.kg.scores.tsv", sep="\t", index=False)
out = pd.DataFrame(rows); out.to_csv(RESULTS/outname, sep="\t", index=False)
