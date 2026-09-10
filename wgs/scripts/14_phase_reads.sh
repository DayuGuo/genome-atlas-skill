#!/bin/bash
# Read-backed phasing of PASS autosomal variants with WhatsHap, 8 chromosomes at a time.
set -euo pipefail
source "$(dirname "$0")/env.sh"
WH=$HOME/miniforge3/envs/whatshap/bin/whatshap; W=$WGS/10_phase; mkdir -p $W; cd $W
run(){ c=$1; bcftools view -r $c $WGS/00_input/target.pass.vcf.gz -Oz -o in.$c.vcf.gz && tabix -f in.$c.vcf.gz && \
  $WH phase --ignore-read-groups --indels --reference $REF -o wh.$c.vcf.gz in.$c.vcf.gz $CRAM > $WGS/logs/wh.$c.log 2>&1 && tabix -f wh.$c.vcf.gz && rm -f in.$c.vcf.gz*; }
export -f run; export WGS REF CRAM WH
printf "%s\n" $(seq 1 22) | xargs -P 8 -I{} bash -c 'run {}'
bcftools concat -Oz -o target.pass.readphased.vcf.gz $(for c in $(seq 1 22); do echo wh.$c.vcf.gz; done); tabix -f target.pass.readphased.vcf.gz
$WH stats --tsv=readphased.stats.tsv target.pass.readphased.vcf.gz > /dev/null 2>&1 || true
echo WHATSHAP_DONE
