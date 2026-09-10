#!/bin/bash
# Normalize the delivered VCF: split multiallelics, left-align, keep PASS, drop decoy contigs.
set -euo pipefail
source "$(dirname "$0")/env.sh"
cd $WGS/00_input
# 1. callable BED: depth>=8 & MQ>=20 (from mosdepth quantized), main contigs only
zcat $WGS/01_qc/depth.quantized.bed.gz | awk -F'\t' '$1!~/^GL/ && ($4=="8:20"||$4=="20:60"||$4=="60:inf")' \
  | bedtools merge -i - > callable.bed
awk '{s+=$3-$2} END{printf "callable bp: %d\n", s}' callable.bed
# 2. normalized full VCF (all filters kept, FILTER column retained) on main contigs
MAIN=$(seq 1 22 | tr '\n' ',')X,Y,MT
bcftools view -r $MAIN "$VCF_RAW" -Ou \
 | bcftools norm -f $REF -m -both -c w --threads $THREADS -Oz -o target.norm.vcf.gz 2> norm.log
tabix -f -p vcf target.norm.vcf.gz
# 3. PASS-only
bcftools view -f PASS --threads $THREADS -Oz -o target.pass.vcf.gz target.norm.vcf.gz
tabix -f -p vcf target.pass.vcf.gz
# 4. stats
bcftools stats -F $REF -s - target.norm.vcf.gz > $WGS/01_qc/stats.norm.txt
bcftools stats -F $REF -s - target.pass.vcf.gz > $WGS/01_qc/stats.pass.txt
tail -3 norm.log
echo NORMALIZE_DONE
