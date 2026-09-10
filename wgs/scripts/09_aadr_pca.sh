#!/bin/bash
# PCA built on present-day East Asians (Human Origins panel); ancient individuals and the target projected in.
set -euo pipefail
source "$(dirname "$0")/env.sh"
W=$WGS/11_aadr; cd $W
awk -F'\t' 'NR>1 && $3=="modern"{print $2"\t"$1}' samples.tsv > modern.keep
$PLINK2 --bfile aadr --keep modern.keep --geno 0.02 --maf 0.01 --indep-pairwise 200 25 0.4 --out prune --threads 16 --memory 40000 >/dev/null
$PLINK2 --bfile aadr --keep modern.keep --extract prune.prune.in --freq --pca 10 allele-wts --out pca --threads 16 --memory 40000 >/dev/null
$PLINK2 --bfile aadr --extract prune.prune.in --read-freq pca.afreq \
  --score pca.eigenvec.allele 2 6 header-read no-mean-imputation variance-standardize \
  --score-col-nums 7-16 --out proj --threads 16 --memory 40000 >/dev/null
echo "pruned SNPs: $(wc -l < prune.prune.in)"; head -2 proj.sscore | cut -f1-8
