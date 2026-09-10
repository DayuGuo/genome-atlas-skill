#!/bin/bash
# Annotate PASS VCF with ClinVar (pos+ref+alt), Ensembl consequences (bcftools csq), 1000G EAS AF.
set -euo pipefail
source "$(dirname "$0")/env.sh"
A=$REF_DIR/annot; W=$WGS/05_clinvar; mkdir -p $W; cd $W
CV=$REF_DIR/clinvar_grch37.vcf.gz
[ -f $CV.tbi ] || tabix -p vcf $CV
# 1. EAS AF annotation file from plink afreq
paste <(awk 'NR>1{print $1"\t"$2"\t"$4"\t"$5"\t"$6}' $WGS/02_complete/eas_allvar.afreq) <(awk 'NR>1{print $6}' $WGS/02_complete/all_allvar.afreq) | bgzip -@4 > eas_af.tsv.gz; tabix -f -s1 -b2 -e2 eas_af.tsv.gz
printf '##INFO=<ID=EAS_AF,Number=1,Type=Float,Description="1000G phase3 EAS alt allele frequency">\n##INFO=<ID=ALL_AF,Number=1,Type=Float,Description="1000G phase3 global alt allele frequency">\n' > eas_af.hdr
# 2. ClinVar + EAS AF + csq
bcftools annotate -a $CV -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/GENEINFO,INFO/CLNSIGCONF,INFO/CLNVC,INFO/ALLELEID --threads 8 -Ou $WGS/00_input/target.pass.vcf.gz \
 | bcftools annotate -a eas_af.tsv.gz -h eas_af.hdr -c CHROM,POS,REF,ALT,INFO/EAS_AF,INFO/ALL_AF -Ou \
 | bcftools csq -f $REF -g $A/Homo_sapiens.GRCh37.87.gff3.gz -p a --ncsq 32 -l --threads 8 -Oz -o target.pass.annot.vcf.gz
tabix -f -p vcf target.pass.annot.vcf.gz
# 3. tables
bcftools query -i 'INFO/CLNSIG!=""' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/GENEINFO\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\t%INFO/CLNSIGCONF\t%INFO/EAS_AF\t%INFO/ALL_AF\t%INFO/BCSQ\t[%GT\t%DP\t%GQ\t%AD]\n' target.pass.annot.vcf.gz > clinvar_all_hits.tsv
bcftools query -i 'INFO/BCSQ~"stop_gained" || INFO/BCSQ~"frameshift" || INFO/BCSQ~"splice_acceptor" || INFO/BCSQ~"splice_donor" || INFO/BCSQ~"start_lost"' \
  -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%INFO/EAS_AF\t%INFO/ALL_AF\t%INFO/CLNSIG\t%INFO/BCSQ\t[%GT\t%DP\t%GQ\t%AD]\n' target.pass.annot.vcf.gz > lof_all.tsv
wc -l clinvar_all_hits.tsv lof_all.tsv
echo ANNOT_DONE
