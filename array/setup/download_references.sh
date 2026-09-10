#!/usr/bin/env bash
# Download every reference the pipeline needs into $REF (config: ref_dir). ~45 GB. Re-runnable (skips existing files).
set -euo pipefail
source "$(dirname "$0")/../scripts/env.sh"
mkdir -p "$REF/imp" "$REF/chain" "$TOOLS"; cd "$REF"
get(){ [ -s "$2" ] && { echo "have $2"; return; }; echo "get $2"; curl -sL -o "$2" "$1"; }
# 1000 Genomes phase 3 plink2 bundle (GRCh37) — https://www.cog-genomics.org/plink/2.0/resources
get "https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1" all_phase3.pgen.zst
get "https://www.dropbox.com/s/c95n8quqwqww4s0/all_phase3_noannot.pvar.zst?dl=1" all_phase3_noannot.pvar.zst
get "https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6&dl=1" all_phase3.psam
[ -s all_phase3.pgen ] || "$PLINK2" --zst-decompress all_phase3.pgen.zst all_phase3.pgen
[ -s all_phase3.pvar ] || "$PLINK2" --zst-decompress all_phase3_noannot.pvar.zst all_phase3.pvar
# HGDP (Bergström 2020, GRCh38) plink2 bundle
get "https://www.dropbox.com/s/hppj1g1gzygcocq/hgdp_all.pgen.zst?dl=1" hgdp_all.pgen.zst
get "https://www.dropbox.com/s/mypl0wgpmhxwgjg/hgdp_all_noannot.pvar.zst?dl=1" hgdp_all_noannot.pvar.zst
get "https://www.dropbox.com/s/0zg57558fqpj3w1/hgdp.psam?dl=1" hgdp_all.psam
[ -s hgdp_all.pgen ] || "$PLINK2" --zst-decompress hgdp_all.pgen.zst hgdp_all.pgen
[ -s hgdp_all.pvar ] || "$PLINK2" --zst-decompress hgdp_all_noannot.pvar.zst hgdp_all.pvar
# ClinVar GRCh37, rCRS, UCSC chain
get "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz" clinvar_grch37.vcf.gz
get "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text" rCRS.fasta
get "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" chain/hg19ToHg38.over.chain.gz
# Beagle 5.4 + genetic maps + 1000G phase3 VCFs (imputation reference, ~15 GB)
get "https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar" imp/beagle.jar
get "https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip" imp/plink.GRCh37.map.zip
[ -d imp/maps ] || unzip -oq imp/plink.GRCh37.map.zip -d imp/maps
for c in $(seq 1 22); do get "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$c.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" "imp/ALL.chr$c.vcf.gz"; done
# haplogrep3
if [ ! -x "$TOOLS/haplogrep3" ]; then cd "$TOOLS" && curl -sL -o hg3.zip "https://github.com/genepi/haplogrep3/releases/download/v3.2.2/haplogrep3-3.2.2-linux.zip" && unzip -oq hg3.zip && rm hg3.zip; fi
echo "references ready in $REF"
