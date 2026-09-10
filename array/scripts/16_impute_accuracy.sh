#!/usr/bin/env bash
# Imputation accuracy check: mask 5% of typed sites on chr20 and chr22, impute, compare to the array calls.
set -euo pipefail
source "$(dirname "$0")/env.sh"
R=$REF/imp; O=$RESULTS2/impute_check; mkdir -p $O; cd "$WORK"
for ch in 20 22; do
python3 - $ch <<'PY'
import sys, random
ch=sys.argv[1]; random.seed(20260908)
lines=[l for l in open(f"results_v2/imputed/target.chr{ch}.vcf")]
hdr=[l for l in lines if l.startswith("#")]; body=[l for l in lines if not l.startswith("#")]
mask=set(random.sample(range(len(body)), int(0.05*len(body))))
with open(f"results_v2/impute_check/target.chr{ch}.masked.vcf","w") as f:
    f.writelines(hdr); f.writelines(l for i,l in enumerate(body) if i not in mask)
with open(f"results_v2/impute_check/truth.chr{ch}.tsv","w") as f:
    for i in sorted(mask):
        t=body[i].rstrip("\n").split("\t"); f.write(f"{t[0]}\t{t[1]}\t{t[3]}\t{t[4]}\t{t[9]}\n")
print(ch, "masked", len(mask), "of", len(body))
PY
java -Xmx${BEAGLE_MEM:-32g} -jar $R/beagle.jar gt=$O/target.chr$ch.masked.vcf ref=$R/ALL.chr$ch.vcf.gz map=$R/maps/plink.chr$ch.GRCh37.map out=$O/imp.chr$ch nthreads=$THREADS gp=true > $O/chr$ch.log 2>&1
done
python3 - <<'PY'
import gzip, os, pandas as pd
W=os.environ["WORK"]; rows=[]
for ch in ["20","22"]:
    truth={ (l.split("\t")[1]): l.rstrip("\n").split("\t") for l in open(f"{W}/results_v2/impute_check/truth.chr{ch}.tsv")}
    with gzip.open(f"{W}/results_v2/impute_check/imp.chr{ch}.vcf.gz","rt") as fh:
        for l in fh:
            if l[0]=="#": continue
            t=l.rstrip("\n").split("\t")
            if t[1] not in truth: continue
            tr=truth[t[1]]
            if (t[3],t[4])!=(tr[2],tr[3]): continue
            info=dict(kv.split("=") for kv in t[7].split(";") if "=" in kv); dr2=float(info.get("DR2",0))
            fmt=t[8].split(":"); val=t[9].split(":"); gt=val[fmt.index("GT")].replace("|","/"); gp=max(float(x) for x in val[fmt.index("GP")].split(","))
            tg=sorted(tr[4].split("/")); ig=sorted(gt.split("/"))
            rows.append(dict(chrom=ch,pos=int(t[1]),dr2=dr2,gp=gp,true="/".join(tg),imp="/".join(ig),correct=tg==ig,het=tg[0]!=tg[1]))
d=pd.DataFrame(rows); d.to_csv(f"{W}/results_v2/impute_check/masked_site_results.tsv",sep="\t",index=False)
out=[]
for lo,hi in [(0,0.7),(0.7,0.9),(0.9,0.99),(0.99,1.01)]:
    s=d[(d.gp>=lo)&(d.gp<hi)]
    if len(s): out.append(dict(gp_bin=f"[{lo},{hi})",n=len(s),concordance=round(s.correct.mean(),4),het_concordance=round(s[s.het].correct.mean(),4) if s.het.any() else None))
o=pd.DataFrame(out); o.to_csv(f"{W}/results_v2/impute_check/summary.tsv",sep="\t",index=False); print(o.to_string(index=False))
print("overall n", len(d), "concordance", round(d.correct.mean(),4), "het", round(d[d.het].correct.mean(),4), "| GP>=0.99: n", (d.gp>=0.99).sum(), "concordance", round(d[d.gp>=0.99].correct.mean(),4))
print("NOTE: single-sample DR2 tracks the sample's own genotype (≈0 at hom-ref sites); never filter imputed variants on it.")
PY
