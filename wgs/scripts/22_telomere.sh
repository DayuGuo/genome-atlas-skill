#!/bin/bash
# Telomere length from the telomeric-read fraction (TelSeq-style), using fixed-string matching for speed.
# A read counts as telomeric when it contains k or more tandem TTAGGG (or CCCTAA) hexamers.
set -euo pipefail
source "$(dirname "$0")/env.sh"
W=$WGS/14_telomere; mkdir -p $W
# TelSeq counts the number of TTAGGG (or CCCTAA) occurrences anywhere in a read, not a contiguous run.
samtools view -@ 16 -T "$REF" -F 0x100 -F 0x400 -F 0x200 "$CRAM" \
 | awk '
   {n++; s=$10; bp+=length(s); m=0
    if (index(s,"TTAGGG")) { x=s; m=gsub(/TTAGGG/,"",x) }
    if (index(s,"CCCTAA")) { y=s; q=gsub(/CCCTAA/,"",y); if (q>m) m=q }
    if (m>=7)  k7++
    if (m>=10) k10++
    if (m>=12) k12++
    if (m>=14) k14++}
   END{printf "total_reads\t%d\ntotal_bp\t%d\nmean_readlen\t%.1f\ntel_reads_k7\t%d\ntel_reads_k10\t%d\ntel_reads_k12\t%d\ntel_reads_k14\t%d\n", n, bp, bp/n, k7, k10, k12, k14}' > "$W/counts.tsv"
cat "$W/counts.tsv"
awk -v cov=28.88 -v ends=92 'BEGIN{FS="\t"} {v[$1]=$2} END{
  rl=v["mean_readlen"];
  split("tel_reads_k7 tel_reads_k10 tel_reads_k12 tel_reads_k14",ks," ");
  for (i=1;i<=4;i++) { k=ks[i]; tl=v[k]*rl/cov/ends;
    printf "%s\tmean_telomere_len_bp\t%.0f\ttelomeric_read_fraction\t%.4g\n", k, tl, v[k]/v["total_reads"] }
  }' "$W/counts.tsv" | tee "$W/estimate.tsv"
echo TELOMERE_DONE
