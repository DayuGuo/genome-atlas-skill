#!/bin/bash
# Local ancestry along the genome with FLARE 0.6: northern vs southern East Asian panels from 1000G,
# plus European and South Asian panels as a noise floor.
set -euo pipefail
source "$(dirname "$0")/env.sh"
JAVA=$HOME/miniforge3/envs/wgs/bin/java
W=$WGS/12_localanc; mkdir -p $W; cd $W
R=$REF_DIR
# panel map: sample <tab> panel
awk 'NR>1 {p=""; if($6=="CHB"||$6=="JPT") p="NorthEA"; else if($6=="CDX"||$6=="KHV") p="SouthEA"; else if($6=="CEU"||$6=="GBR") p="European"; else if($6=="GIH"||$6=="PJL") p="SouthAsian"; if(p!="") print $1"\t"p}' $R/all_phase3.psam > ref.panel
cut -f2 ref.panel | sort | uniq -c
for c in $(seq 1 22); do
  [ -f $WGS/10_phase/ph.$c.vcf.gz ] || { echo "chr$c not phased yet, skipping"; continue; }
  $JAVA -Xmx40g -jar $FLARE ref=$R/imp/ALL.chr$c.vcf.gz ref-panel=ref.panel gt=$WGS/10_phase/ph.$c.vcf.gz \
    map=$R/imp/maps/plink.chr$c.GRCh37.map out=la.$c probs=true gen=100 nthreads=12 > $WGS/logs/flare.$c.log 2>&1
  echo "chr$c done: $(bcftools view -H la.$c.anc.vcf.gz 2>/dev/null | wc -l) markers"
done
echo FLARE_DONE
