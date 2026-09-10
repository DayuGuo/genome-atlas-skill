#!/usr/bin/env bash
# Optional: FASTQ -> coordinate-sorted CRAM. Skip this if the provider already delivered a CRAM or BAM.
# Needs bwa-mem2 (or bwa) and samtools; set fastq1 / fastq2 in config.yaml.
set -euo pipefail
source "$(dirname "$0")/../scripts/env.sh"
[ -n "${FASTQ1:-}" ] || { echo "set fastq1/fastq2 in config.yaml first"; exit 1; }
command -v bwa-mem2 >/dev/null || "$TOOLS/env/bin/mamba" install -y -p "$TOOLS/env" -c bioconda bwa-mem2
[ -f "$REF.bwt.2bit.64" ] || bwa-mem2 index "$REF"
bwa-mem2 mem -t "$THREADS" -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:illumina" "$REF" "$FASTQ1" "$FASTQ2" \
 | samtools fixmate -@ 4 -m - - \
 | samtools sort -@ "$THREADS" -T "$RAW/sort" - \
 | samtools markdup -@ 4 - -O cram --reference "$REF" "$RAW/$SAMPLE.cram"
samtools index -@ 4 "$RAW/$SAMPLE.cram"
samtools flagstat -@ 4 "$RAW/$SAMPLE.cram" > "$RAW/$SAMPLE.flagstat.txt"
echo "wrote $RAW/$SAMPLE.cram — set reads: to this path in config.yaml"
