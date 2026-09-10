#!/usr/bin/env python3
"""PRS validation: for each PGS, score the 504 1000G EAS samples with (a) the FULL score matched to all
1000G biallelic variants and (b) the array-subset score already computed. Report coverage, Spearman
correlation, and the distribution of |percentile_subset - percentile_full| so the partial-score error is
measured rather than asserted. Output: audit/04_prs_validation.tsv"""
import gzip, json, os, pathlib, subprocess, numpy as np, pandas as pd
from scipy.stats import spearmanr
from config import *
ROOT = pathlib.Path(__file__).resolve().parents[1]; A = RESULTS/"ancestry"; W = RESULTS/"prs"; AU = AUDIT; FV = AU/"prs_full"; FV.mkdir(parents=True, exist_ok=True)
P = str(pathlib.Path(PLINK2)); R = REF
meta = {s["id"]: s for s in json.load(open(DATA/"pgs"/"all_scores.json"))}
# full 1000G biallelic SNP index (chr:pos -> ref, alt); built once (~80M rows -> keep only SNPs)
idx_f = FV/"kg_all_snps.parquet"
if not idx_f.exists():
    chunks = []
    for ch in pd.read_csv(R/"all_phase3.pvar", sep="\t", comment="#", header=None, usecols=[0,1,3,4], names=["chrom","pos","ref","alt"], dtype={"chrom":str}, chunksize=5_000_000):
        ch = ch[(ch.ref.str.len()==1)&(ch.alt.str.len()==1)&ch.chrom.isin([str(i) for i in range(1,23)])]
        chunks.append(ch)
    idx = pd.concat(chunks); idx = idx[~idx.duplicated(["chrom","pos"], keep=False)]  # drop multi-allelic split positions
    idx.to_parquet(idx_f)
idx = pd.read_parquet(idx_f)
if not (FV/"kg_all.pgen").exists():
    idx.assign(id=idx.chrom+":"+idx.pos.astype(str))[["chrom","pos","pos","id"]].to_csv(FV/"snp.range", sep="\t", header=False, index=False)
    subprocess.run([P,"--pfile",R/"all_phase3","--keep",A/"eas.ids","--snps-only","just-acgt","--max-alleles","2","--autosome","--extract","range",FV/"snp.range","--rm-dup","exclude-all","--set-all-var-ids","@:#","--make-pgen","--out",FV/"kg_all","--threads", THREADS],check=True,capture_output=True)
raw_pos = pd.read_parquet(DATA/"genome.parquet"); raw_pos["chrom"]=raw_pos.chrom.astype(str); raw_set = set(zip(raw_pos.chrom, raw_pos.pos))
sub_ids = set(l.strip() for l in open(A/"common.ids"))
rows = []
for pid in sorted(x.name.split(".")[0] for x in W.glob("PGS*.score") if "noAPOE" not in x.name):
    m = meta.get(pid, {}); f = DATA/"pgs"/f"{pid}_hmPOS_GRCh37.txt.gz"
    hdr = {}
    with gzip.open(f, "rt") as fh:
        for l in fh:
            if not l.startswith("#"): break
            if "=" in l: k, v = l[1:].strip().split("=", 1); hdr[k] = v
    s = pd.read_csv(f, sep="\t", comment="#", dtype={"hm_chr":str}, low_memory=False).dropna(subset=["hm_chr","hm_pos"]); s["hm_pos"] = s.hm_pos.astype(int)
    if "other_allele" not in s: s["other_allele"] = np.nan
    n_total = len(s); n_snv = int((s.effect_allele.str.len()==1).sum())
    on_array = s[[ (c,p) in raw_set for c,p in zip(s.hm_chr, s.hm_pos)]]
    pal = on_array[(on_array.effect_allele.isin(["A","T"]) & on_array.other_allele.isin(["A","T"])) | (on_array.effect_allele.isin(["C","G"]) & on_array.other_allele.isin(["C","G"]))]
    used = sum(1 for _ in open(W/f"{pid}.score"))
    # full score: match to all 1000G biallelic SNPs
    mm = s.merge(idx, left_on=["hm_chr","hm_pos"], right_on=["chrom","pos"])
    ok = ((mm.effect_allele==mm.ref)&((mm.other_allele==mm.alt)|mm.other_allele.isna()))|((mm.effect_allele==mm.alt)&((mm.other_allele==mm.ref)|mm.other_allele.isna()))
    mm = mm[ok]; mm["id"] = mm.chrom+":"+mm.pos.astype(str)
    sf = FV/f"{pid}.full.score"; mm[["id","effect_allele","effect_weight"]].to_csv(sf, sep="\t", index=False, header=False)
    if not (FV/f"{pid}.full.sscore").exists():
        subprocess.run([P,"--pfile",FV/"kg_all","--score",sf,"1","2","3","cols=+scoresums","no-mean-imputation","--out",FV/f"{pid}.full","--threads", THREADS],check=True,capture_output=True)
    full = pd.read_csv(FV/f"{pid}.full.sscore", sep="\t").rename(columns={"#IID":"IID"})
    sub = pd.read_csv(W/f"{pid}.kg.sscore", sep="\t").rename(columns={"#IID":"IID"}); sub = sub[sub.SuperPop==REF_SUPERPOP]
    j = full[["IID","SCORE1_SUM"]].merge(sub[["IID","SCORE1_SUM"]], on="IID", suffixes=("_full","_sub"))
    rho = spearmanr(j.SCORE1_SUM_full, j.SCORE1_SUM_sub).correlation
    pf = j.SCORE1_SUM_full.rank(pct=True)*100; ps = j.SCORE1_SUM_sub.rank(pct=True)*100; dev = (pf-ps).abs()
    # the sample subset percentile and its plausible full-score range from the empirical deviation
    d = pd.read_csv(W/f"{pid}.sample.sscore", sep="\t").SCORE1_SUM.iloc[0]; pct_sub = (sub.SCORE1_SUM < d).mean()*100
    # top/bottom decile retention: of samples in top 10% by full score, share still in top 20% by subset
    top = pf >= 90; ret = (ps[top] >= 80).mean()*100
    dev_anc = m.get("ancestry_distribution",{}).get("dev",{}).get("dist",{})
    rows.append(dict(pgs=pid, trait=m.get("trait_reported",""), publication=(m.get("publication",{}).get("firstauthor","")+" "+str(m.get("publication",{}).get("date_publication",""))[:4]),
                     dev_ancestry=json.dumps(dev_anc), weight_type=hdr.get("weight_type",""), build=hdr.get("HmPOS_build",""),
                     variants_total=n_total, variants_snv=n_snv, at_array_positions=len(on_array), palindromic_at_array=len(pal), used_direct=used,
                     coverage_pct=round(100*used/n_total,1), full_matched_1000G=len(mm), full_match_pct=round(100*len(mm)/n_total,1),
                     spearman_full_vs_subset=round(rho,3), pct_dev_median=round(dev.median(),1), pct_dev_p90=round(dev.quantile(.9),1), pct_dev_max=round(dev.max(),1),
                     top10_full_retained_in_top20_subset_pct=round(ret,0), sample_pct_subset=round(pct_sub,1),
                     reliability=("do_not_present_as_formal_PRS" if used/n_total < 0.5 else "caution")))
    print(f"{pid} cov {100*used/n_total:4.1f}%  full-match {100*len(mm)/n_total:4.1f}%  rho={rho:.3f}  |Δpct| median {dev.median():.1f} p90 {dev.quantile(.9):.1f}  top10→top20 {ret:.0f}%", flush=True)
    pd.DataFrame(rows).to_csv(AU/"04_prs_validation.tsv", sep="\t", index=False)
