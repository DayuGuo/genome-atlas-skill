#!/bin/bash
# 1000G PCA on LD-pruned common autosomal SNVs; project the target (complete genotype set) with --score.
set -euo pipefail
source "$(dirname "$0")/env.sh"
W=$WGS/04_ancestry; mkdir -p $W; cd $W
R=$KG_PFILE
# candidate sites: autosomal biallelic SNVs, global MAF>=0.05, present (non-missing) in the target's complete set
awk 'NR>1 && $1!="X" && $3!="." {print $3}' $WGS/02_complete/target.1kg.pvar > target.present.ids   # ~80M ids
$PLINK2 --pfile $R --autosome --snps-only just-acgt --max-alleles 2 --maf 0.05 --extract target.present.ids \
  --rm-dup exclude-all --make-pgen --out kg.common --threads 16 --memory 40000 >/dev/null
# drop strand-ambiguous not needed (same reference build, same alleles); thin then LD-prune
$PLINK2 --pfile kg.common --bp-space 2000 --indep-pairwise 200 50 0.2 --out prune --threads 16 --memory 40000 >/dev/null
$PLINK2 --pfile kg.common --extract prune.prune.in --freq --pca 10 allele-wts --out kg.pca --threads 16 --memory 40000 >/dev/null
for s in kg:kg.common target:$WGS/02_complete/target.1kg; do n=${s%%:*}; f=${s#*:}
  $PLINK2 --pfile $f --extract prune.prune.in --read-freq kg.pca.afreq \
     --score kg.pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
     --score-col-nums 6-15 --out $n.proj --threads 16 --memory 40000 >/dev/null
done
echo "pruned SNPs: $(wc -l < prune.prune.in)"; head -2 target.proj.sscore
# EAS-only PCA
awk 'NR>1 && $5=="EAS"{print $1}' $R.psam > eas.ids; sed -i '1i #IID' eas.ids
$PLINK2 --pfile kg.common --keep eas.ids --maf 0.05 --bp-space 2000 --indep-pairwise 200 50 0.2 --out prune.eas --threads 16 --memory 40000 >/dev/null
$PLINK2 --pfile kg.common --keep eas.ids --extract prune.eas.prune.in --freq --pca 10 allele-wts --out eas.pca --threads 16 --memory 40000 >/dev/null
$PLINK2 --pfile kg.common --keep eas.ids --extract prune.eas.prune.in --read-freq eas.pca.afreq --score eas.pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize --score-col-nums 6-15 --out eas.proj --threads 16 --memory 40000 >/dev/null
$PLINK2 --pfile $WGS/02_complete/target.1kg --extract prune.eas.prune.in --read-freq eas.pca.afreq --score eas.pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize --score-col-nums 6-15 --out eas.target.proj --threads 16 --memory 40000 >/dev/null
echo "EAS pruned SNPs: $(wc -l < prune.eas.prune.in)"
echo PCA_DONE
