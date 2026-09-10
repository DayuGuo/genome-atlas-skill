#!/usr/bin/env bash
source "$(dirname "$0")/env.sh"; cd "$REF/imp"
for c in $(seq 1 22); do
  ( zcat ALL.chr$c.vcf.gz | awk -F'\t' '/^#/ {print; next} !seen[$2":"$4":"$5]++' | gzip -1 > ALL.chr$c.dedup.vcf.gz && mv ALL.chr$c.dedup.vcf.gz ALL.chr$c.vcf.gz && echo "chr$c deduped" ) &
done
wait; echo ALL_DEDUP_DONE
