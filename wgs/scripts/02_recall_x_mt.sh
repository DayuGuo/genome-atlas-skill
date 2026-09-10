#!/bin/bash
# Re-call chrX with correct ploidy (PAR diploid, non-PAR haploid for a male) and pile up MT for heteroplasmy.
set -euo pipefail
source "$(dirname "$0")/env.sh"
mkdir -p $WGS/03_haplo $WGS/00_input
cd $WGS/00_input
echo -e "41240811001690\tM" > sample_sex.txt
# X: parallel over 4 chunks
for r in X:1-40000000 X:40000001-80000000 X:80000001-120000000 X:120000001-155270560; do
  ( bcftools mpileup -f $REF -r $r -q 20 -Q 20 -d 500 -a FORMAT/AD,FORMAT/DP -Ou "$CRAM" 2>/dev/null \
    | bcftools call -m -v --ploidy GRCh37 -S sample_sex.txt -Oz -o X.$(echo $r|tr ':-' '__').vcf.gz 2>/dev/null ) &
done
wait
bcftools concat -a -Oz -o X.recall.raw.vcf.gz X.X_*.vcf.gz 2>/dev/null || { for f in X.X_*.vcf.gz; do tabix -f $f; done; bcftools concat -a -Oz -o X.recall.raw.vcf.gz X.X_*.vcf.gz; }
tabix -f X.recall.raw.vcf.gz
bcftools norm -f $REF -m -both -Ou X.recall.raw.vcf.gz | bcftools filter -i 'QUAL>=30 && FORMAT/DP>=8' -Oz -o X.recall.vcf.gz
tabix -f X.recall.vcf.gz
rm -f X.X_*.vcf.gz*
echo "PAR1 het/hom:"; bcftools query -r X:60001-2699520 -f '[%GT]\n' X.recall.vcf.gz | sort | uniq -c
echo "nonPAR GTs:"; bcftools query -r X:2699521-154931043 -f '[%GT]\n' X.recall.vcf.gz | sort | uniq -c
# MT pileup: all positions, allele depths
bcftools mpileup -f $REF -r MT -q 20 -Q 20 -d 20000 -a FORMAT/AD -Ou "$CRAM" 2>/dev/null \
  | bcftools query -f '%POS\t%REF\t%ALT\t[%AD]\n' > $WGS/03_haplo/mt_pileup_ad.tsv
wc -l $WGS/03_haplo/mt_pileup_ad.tsv
echo RECALL_DONE
