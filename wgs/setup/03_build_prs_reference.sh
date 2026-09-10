#!/usr/bin/env bash
# Build the scoring reference: the chosen super-population from 1000G, with chrom:pos variant IDs and an
# allele-frequency table. Every polygenic score is computed for the target and for these individuals on the
# identical variant set, so a percentile means "among these N people".
set -euo pipefail
source "$(dirname "$0")/../scripts/env.sh"
mkdir -p "$REF_DIR/prs"; cd "$REF_DIR/prs"
awk -v sp="$SUPERPOP" 'NR==1{print "#IID"} NR>1 && $5==sp {print $1}' "$KG_PFILE.psam" > ref.ids
"$PLINK2" --pfile "$KG_PFILE" --keep ref.ids --max-alleles 2 --snps-only just-acgt \
  --set-all-var-ids '@:#' --rm-dup exclude-all --make-pgen --out kg_all --threads "$THREADS" --memory $((MEM_GB*1000)) >/dev/null
"$PLINK2" --pfile kg_all --freq --out kg_all --threads "$THREADS" --memory $((MEM_GB*1000)) >/dev/null
# carry the population labels through so percentiles can be split by sub-population
python3 - <<PY
import pandas as pd
src = pd.read_csv("$KG_PFILE.psam", sep="\t").rename(columns={"#IID":"IID"})
dst = pd.read_csv("kg_all.psam", sep="\t").rename(columns={"#IID":"IID"})
m = dst[["IID"]].merge(src[["IID","SuperPop","Population"]], on="IID", how="left")
m.insert(1,"SEX",0); m.columns = ["#IID","SEX","SuperPop","Population"]
m.to_csv("kg_all.psam", sep="\t", index=False)
PY
echo "PRS reference: $(grep -vc '^#' kg_all.psam) individuals from $SUPERPOP"
