#!/usr/bin/env python
"""PRS from WGS complete genotype set (hard calls; hom-ref filled in callable regions), same score files,
scored against the reference super-population (built by setup/03_build_prs_reference.sh) with a reference-MAF>=0.01 filter."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import json, os, pathlib, subprocess, numpy as np, pandas as pd
ROOT = P; FV = pathlib.Path(PRS_REF).parent; W = ROOT/"wgs"/"07_prs"; W.mkdir(parents=True, exist_ok=True)
PLINK = PLINK2; MAF_MIN = 0.01
def run(*a):
    r = subprocess.run([PLINK, *map(str, a), "--threads", "16", "--memory", "30000"], capture_output=True, text=True)
    if r.returncode: raise SystemExit(r.stdout[-3000:] + r.stderr[-3000:])
if not (W/f"{SAMPLE}_pos.pgen").exists():
    run("--pfile", ROOT/f"wgs/02_complete/{SAMPLE}.1kg", "--autosome", "--set-all-var-ids", "@:#", "--rm-dup", "exclude-all", "--make-pgen", "--out", W/f"{SAMPLE}_pos")
ids = set(l.split("\t")[2] for l in open(W/f"{SAMPLE}_pos.pvar") if not l.startswith("#"))
fr = pd.read_csv(FV/"kg_all.afreq", sep="\t", usecols=["ID", "ALT_FREQS"]); fr["maf"] = np.minimum(fr.ALT_FREQS, 1-fr.ALT_FREQS)
common = set(fr[fr.maf >= MAF_MIN].ID)
usable = ids & common
print(f"WGS complete-set autosomal SNVs: {len(ids):,}; EAS-common usable: {len(usable):,}")
meta = {s["id"]: s for s in json.load(open(ROOT/"data/pgs/all_scores.json"))}
rows = []
for sf in sorted(FV.glob("PGS*.full.score")):
    pid = sf.name.split(".")[0]; s = pd.read_csv(sf, sep="\t", header=None, names=["id", "ea", "w"]); n_total = meta[pid]["variants_number"]
    keep = s[s.id.isin(usable)]; keep.to_csv(W/f"{pid}.score", sep="\t", header=False, index=False)
    for who, pf in [("eas", FV/"kg_all"), ("target", W/f"{SAMPLE}_pos")]:
        run("--pfile", pf, "--score", W/f"{pid}.score", "1", "2", "3", "cols=+scoresums", "no-mean-imputation", "--out", W/f"{pid}.{who}")
    eas = pd.read_csv(W/f"{pid}.eas.sscore", sep="\t"); d = pd.read_csv(W/f"{pid}.target.sscore", sep="\t")
    x = d.SCORE1_SUM.iloc[0]; han = eas[eas.Population.isin(SUBPOPS)].SCORE1_SUM
    z = (x-eas.SCORE1_SUM.mean())/eas.SCORE1_SUM.std(); pct = (eas.SCORE1_SUM < x).mean()*100
    rows.append(dict(pgs=pid, trait=meta[pid]["trait_reported"], variants_in_score=n_total, used_wgs=len(keep), coverage_pct=round(100*len(keep)/n_total, 1),
                     weight_retained_pct=round(100*keep.w.abs().sum()/s.w.abs().sum(), 1), sample_score=x, z_vs_EAS=round(z, 2), pct_EAS=round(pct, 1), pct_sub=round((han < x).mean()*100, 1)))
    print(f"{pid} {meta[pid]['trait_reported'][:28]:28s} used {len(keep):>8,}/{n_total:<9,} ({100*len(keep)/n_total:5.1f}%) pct={pct:5.1f} z={z:+.2f}", flush=True)
df = pd.DataFrame(rows); df.to_csv(W/"prs_wgs.tsv", sep="\t", index=False)
print(df[["pgs","trait","coverage_pct","pct_EAS","z_vs_EAS"]].to_string(index=False))
