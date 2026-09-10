#!/usr/bin/env python3
"""PRS v2 on imputed dosages (Beagle 5.4, 1000G phase3 reference).

Variant selection is REFERENCE-based only (1000G EAS MAF >= 0.01), never on the sample's own DR2/GP: with a single
target sample Beagle's DR2/AF are functions of the sample's genotype (DR2 ~ 0 at homozygous-reference sites), so
filtering on them would preferentially keep alt-carrying sites and bias every score upward.
Main analysis: dosages, all reference-common variants. Sensitivity: additionally restrict to the sample max-GP >= 0.9.
Reference: the 504 1000G EAS samples scored on the identical variant set. Output: results_v2/prs_imputed.tsv"""
import json, os, pathlib, subprocess, numpy as np, pandas as pd
from config import *
ROOT = pathlib.Path(__file__).resolve().parents[1]; A = RESULTS/"ancestry"; FV = AUDIT/"prs_full"; O = RESULTS2; I = O/"imputed"; W = O/"prs"; W.mkdir(parents=True, exist_ok=True)
P = str(pathlib.Path(PLINK2)); MAF_MIN = 0.01
def run(*a):
    r = subprocess.run([P, *map(str, a), "--threads", THREADS], capture_output=True, text=True)
    if r.returncode: raise SystemExit(r.stdout[-3000:] + r.stderr[-3000:])

if not (W/"sample_imp.pgen").exists():
    parts = []
    for ch in range(1, 23):
        run("--vcf", I/f"chr{ch}.vcf.gz", "dosage=DS", "--import-max-alleles", "2", "--snps-only", "just-acgt", "--max-alleles", "2", "--set-all-var-ids", "@:#", "--rm-dup", "exclude-all", "--make-pgen", "--out", W/f"imp.chr{ch}")
        parts.append(str(W/f"imp.chr{ch}"))
    with open(W/"pmerge.list", "w") as f: f.write("\n".join(parts) + "\n")
    run("--pmerge-list", W/"pmerge.list", "--make-pgen", "--out", W/"sample_imp")
    for ch in range(1, 23):
        for ext in (".pgen", ".pvar", ".psam", ".log"): (W/f"imp.chr{ch}{ext}").unlink(missing_ok=True)
imp_ids = set(l.split("\t")[2] for l in open(W/"sample_imp.pvar") if not l.startswith("#"))
# reference-based frequency filter (EAS 504)
if not (FV/"kg_all.afreq").exists(): run("--pfile", FV/"kg_all", "--freq", "--out", FV/"kg_all")
fr = pd.read_csv(FV/"kg_all.afreq", sep="\t", usecols=["ID","ALT_FREQS"]); fr["maf"] = np.minimum(fr.ALT_FREQS, 1-fr.ALT_FREQS)
common = set(fr[fr.maf >= MAF_MIN].ID)
if not (O/"gp90.ids").exists():  # variant IDs where the sample's max posterior GP >= 0.9 (sensitivity analysis only)
    cmd = " ".join(f"(zcat {I}/chr{c}.vcf.gz | awk -F'\\t' '!/^#/ {{ n=split($10,a,\":\"); split(a[3],g,\",\"); m=g[1]; if(g[2]>m)m=g[2]; if(g[3]>m)m=g[3]; if(m>=0.9) print $1\":\"$2 }}') &" for c in range(1,23)) + " wait"
    subprocess.run(f"({cmd}) > {O}/gp90.ids", shell=True, check=True)
gp90 = set(l.strip() for l in open(O/"gp90.ids"))
typed = set(l.strip() for l in open(A/"common.ids"))
print(f"imputed biallelic SNPs: {len(imp_ids):,}; EAS-common (MAF>={MAF_MIN}): {len(imp_ids & common):,}; of which GP>=0.9: {len(imp_ids & common & gp90):,}; directly typed: {len(imp_ids & typed):,}")
json.dump({"imputed": len(imp_ids), "common": len(imp_ids & common), "common_gp90": len(imp_ids & common & gp90), "typed": len(imp_ids & typed), "maf_min": MAF_MIN}, open(O/"imputed_counts.json","w"))

meta = {s["id"]: s for s in json.load(open(DATA/"pgs"/"all_scores.json"))}
rows = []
for sf in sorted(FV.glob("PGS*.full.score")):
    pid = sf.name.split(".")[0]; s = pd.read_csv(sf, sep="\t", header=None, names=["id","ea","w"]); n_total = meta[pid]["variants_number"]
    res = {}
    for tag, keepset in [("main", imp_ids & common), ("gp90", imp_ids & common & gp90)]:
        keep = s[s.id.isin(keepset)]; keep.to_csv(W/f"{pid}.{tag}.score", sep="\t", header=False, index=False)
        for who, pf in [("eas", FV/"kg_all"), ("sample", W/"sample_imp")]:
            run("--pfile", pf, "--score", W/f"{pid}.{tag}.score", "1", "2", "3", "cols=+scoresums", "no-mean-imputation", "--out", W/f"{pid}.{tag}.{who}")
        eas = pd.read_csv(W/f"{pid}.{tag}.eas.sscore", sep="\t"); d = pd.read_csv(W/f"{pid}.{tag}.sample.sscore", sep="\t")
        x = d.SCORE1_SUM.iloc[0]; han = eas[eas.Population.isin(SUBPOPS)].SCORE1_SUM
        res[tag] = dict(n=len(keep), w=keep.w.abs().sum()/s.w.abs().sum(), typed=int(keep.id.isin(typed).sum()), x=x, z=(x-eas.SCORE1_SUM.mean())/eas.SCORE1_SUM.std(), pct=(eas.SCORE1_SUM < x).mean()*100, pct_han=(han < x).mean()*100)
    m, g = res["main"], res["gp90"]
    rows.append(dict(pgs=pid, trait=meta[pid]["trait_reported"], variants_in_score=n_total, matched_1000G=len(s), used_imputed=m["n"], coverage_of_score_pct=round(100*m["n"]/n_total,1),
                     weight_retained_pct=round(100*m["w"],1), n_typed=m["typed"], n_imputed_only=m["n"]-m["typed"], sample_score=m["x"], z_vs_EAS=round(m["z"],2), pct_EAS=round(m["pct"],1), pct_Han=round(m["pct_han"],1),
                     gp90_used=g["n"], gp90_pct_EAS=round(g["pct"],1), gp90_z=round(g["z"],2)))
    print(f"{pid} {meta[pid]['trait_reported'][:26]:26s} used {m['n']:>8,}/{n_total:<9,} ({100*m['n']/n_total:5.1f}%, w {100*m['w']:5.1f}%)  pct_EAS={m['pct']:5.1f} z={m['z']:+.2f} | GP>=0.9: pct={g['pct']:5.1f}", flush=True)
pd.DataFrame(rows).to_csv(O/"prs_imputed.tsv", sep="\t", index=False)
