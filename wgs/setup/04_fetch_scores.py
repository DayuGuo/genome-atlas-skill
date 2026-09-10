#!/usr/bin/env python
"""Download the PGS Catalog scoring files listed in panel/pgs_scores.tsv and match them to the PRS reference.
Each score becomes data/ref/prs/<PGS>.full.score with chrom:pos / effect allele / weight."""
import sys, pathlib, urllib.request
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
from wgsconfig import *  # noqa
import numpy as np, pandas as pd

PANEL = pathlib.Path(__file__).resolve().parents[1] / "panel" / "pgs_scores.tsv"
FV = pathlib.Path(PRS_REF).parent; FV.mkdir(parents=True, exist_ok=True)
idx = pd.read_csv(f"{PRS_REF}.pvar", sep="\t", comment="#", header=None,
                  names=["chrom", "pos", "id", "ref", "alt"], dtype={"chrom": str})
scores = pd.read_csv(PANEL, sep="\t", comment="#")
meta = {}
for r in scores.itertuples():
    pid = r.pgs
    f = PGS / f"{pid}_hmPOS_GRCh37.txt.gz"
    if not f.exists():
        url = f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pid}/ScoringFiles/Harmonized/{pid}_hmPOS_GRCh37.txt.gz"
        print("downloading", pid, flush=True); urllib.request.urlretrieve(url, f)
    s = pd.read_csv(f, sep="\t", comment="#", dtype={"hm_chr": str}, low_memory=False).dropna(subset=["hm_chr", "hm_pos"])
    s["hm_pos"] = s.hm_pos.astype(int)
    if "other_allele" not in s: s["other_allele"] = np.nan
    mm = s.merge(idx, left_on=["hm_chr", "hm_pos"], right_on=["chrom", "pos"])
    ok = ((mm.effect_allele == mm.ref) & ((mm.other_allele == mm.alt) | mm.other_allele.isna())) | \
         ((mm.effect_allele == mm.alt) & ((mm.other_allele == mm.ref) | mm.other_allele.isna()))
    mm = mm[ok].copy(); mm["sid"] = mm.chrom + ":" + mm.pos.astype(str)
    mm[["sid", "effect_allele", "effect_weight"]].to_csv(FV / f"{pid}.full.score", sep="\t", header=False, index=False)
    meta[pid] = {"id": pid, "trait_reported": r.trait_en, "trait_zh": r.trait_zh,
                 "variants_number": len(s), "panel": r.panel, "source": r.source}
    print(f"{pid} {r.trait_en[:34]:34s} {len(mm):>9,} / {len(s):<9,} matched", flush=True)
import json; json.dump(list(meta.values()), open(PGS / "all_scores.json", "w"), ensure_ascii=False)
print("wrote", PGS / "all_scores.json")
