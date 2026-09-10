#!/bin/bash
# Statistical phasing of the target's complete genotype set (EAS MAF>=0.01 sites) against 1000G phase3 with Beagle 5.4.
set -euo pipefail
source "$(dirname "$0")/env.sh"
W=$WGS/10_phase; mkdir -p $W; cd $W
IMP=$KG_VCF_DIR
awk 'NR>1 && $1!="X" {m=($6<0.5?$6:1-$6); if(m>=0.01) print $1"\t"$2}' $WGS/02_complete/eas_all.afreq > common.sites.tsv
wc -l common.sites.tsv
for c in $(seq 1 22); do
  awk -v c=$c '$1==c' common.sites.tsv > sites.$c.tsv
  bcftools view -T sites.$c.tsv -r $c $WGS/02_complete/target.1kg_sites.vcf.gz -Oz -o gt.$c.vcf.gz
  java -Xmx30g -jar $IMP/beagle.jar gt=gt.$c.vcf.gz ref=$IMP/ALL.chr$c.vcf.gz map=$IMP/maps/plink.chr$c.GRCh37.map impute=false nthreads=8 out=ph.$c > $WGS/logs/beagle_ph.$c.log 2>&1
  tabix -f ph.$c.vcf.gz; rm -f gt.$c.vcf.gz sites.$c.tsv
  echo "chr$c phased: $(bcftools view -H ph.$c.vcf.gz | wc -l) sites"
done
echo BEAGLE_PHASE_DONE
