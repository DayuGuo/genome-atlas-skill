#!/bin/bash
# Calibration: score 20 CHB and 20 CHS reference individuals with the same north/south panels, holding them out.
set -euo pipefail
source "$(dirname "$0")/env.sh"
JAVA=$HOME/miniforge3/envs/wgs/bin/java
W=$WGS/12_localanc; cd $W
R=$REF_DIR
awk 'NR>1 && $6=="CHB"{print $1}' $R/all_phase3.psam | head -20 > holdout.ids
awk 'NR>1 && $6=="CHS"{print $1}' $R/all_phase3.psam | head -20 >> holdout.ids
grep -vwFf holdout.ids ref.panel > ref.calib.panel
cut -f2 ref.calib.panel | sort | uniq -c
for c in 1 2 6 22; do
  bcftools view -S holdout.ids -r $c $R/imp/ALL.chr$c.vcf.gz -Oz -o /dev/null 2>/dev/null || true
  zcat $R/imp/ALL.chr$c.vcf.gz | bcftools view -S holdout.ids -Oz -o ho.$c.vcf.gz --threads 4
  tabix -f ho.$c.vcf.gz
  $JAVA -Xmx40g -jar $FLARE ref=$R/imp/ALL.chr$c.vcf.gz ref-panel=ref.calib.panel gt=ho.$c.vcf.gz \
    map=$R/imp/maps/plink.chr$c.GRCh37.map out=calib.$c gen=100 nthreads=12 > $WGS/logs/flare_calib.$c.log 2>&1
  rm -f ho.$c.vcf.gz*
  echo "chr$c calibration done"
done
echo CALIB_DONE
