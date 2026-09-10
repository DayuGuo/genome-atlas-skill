#!/bin/bash
# PharmCAT 3.4 on a GRCh37 VCF: explicitly fill hom-ref at every PharmCAT position IN GRCh37 (using the b37 base),
# then lift to GRCh38 with Picard (RECOVER_SWAPPED_REF_ALT handles sites whose reference allele changed between builds,
# e.g. CYP3A5 rs776746), then run PharmCAT with Cyrius' CYP2D6 call as an outside call.
set -euo pipefail
source "$(dirname "$0")/env.sh"
W=$WGS/06_pgx/pharmcat; mkdir -p $W; cd $W
PC=$PHARMCAT_DIR
[ -f $REF_DIR/chain/hg38ToHg19.over.chain.gz ] || curl -sSL -o $REF_DIR/chain/hg38ToHg19.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' $PC/pharmcat_positions_3.4.0.vcf.bgz > pos38.tsv
for c in $(seq 1 22) X Y; do echo -e "$c\tchr$c"; done > rename.txt
# the target's PASS calls within +-5kb of PharmCAT positions (b37), chr-renamed
python - <<'PY'
import pysam, subprocess, collections, bisect
from pyliftover import LiftOver
P="$PROJ"; lo=LiftOver(f"{P}/data/ref/chain/hg38ToHg19.over.chain.gz"); fa=pysam.FastaFile(f"{P}/data/ref/b37/human_g1k_v37.fasta")
call=collections.defaultdict(list)
for l in open(f"{P}/wgs/00_input/callable.bed"):
    c,s,e=l.split()[:3]; call[c].append((int(s),int(e)))
def callable_(c,p):
    v=call.get(c,[]); i=bisect.bisect_right([x[0] for x in v],p-1)-1
    return i>=0 and v[i][0]<=p-1<v[i][1]
rows=[]; stats=collections.Counter()
for l in open("pos38.tsv"):
    c38,p38,rid,r38,a38=l.rstrip("\n").split("\t"); r=lo.convert_coordinate(c38,int(p38)-1)
    if not r: stats["unliftable"]+=1; continue
    c37=r[0][0][3:]; p37=r[0][1]+1; b37=fa.fetch(c37,p37-1,p37).upper()
    alts=a38.split(",")
    if b37==r38: alt=alts[0]; stats["same_ref"]+=1
    elif b37 in alts: alt=r38; stats["ref_changed"]+=1
    else: stats["ref_mismatch_skip"]+=1; continue
    if not callable_(c37,p37): stats["uncallable"]+=1; continue
    rows.append((c37,p37,rid,b37,alt))
print("PharmCAT positions ->", dict(stats))
regions=open("regions37.bed","w")
with open("fill37.vcf","w") as f:
    f.write("##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    for c in [str(i) for i in range(1,23)]+["X","Y"]: f.write(f"##contig=<ID={c}>\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t41240811001690\n")
    for c,p,rid,ref,alt in rows:
        f.write(f"{c}\t{p}\t{rid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t0/0\n"); regions.write(f"{c}\t{max(0,p-5001)}\t{p+5000}\n")
regions.close()
PY
sort -k1,1V -k2,2n regions37.bed | bedtools merge -i - > regions37.merged.bed
bcftools view -R regions37.merged.bed $WGS/00_input/target.pass.vcf.gz -Oz -o target_pgx37.vcf.gz; tabix -f target_pgx37.vcf.gz
# drop fill records at positions where the target already has a call, then combine
bcftools query -f '%CHROM\t%POS\n' target_pgx37.vcf.gz > have.tsv
bcftools sort fill37.vcf -Oz -o fill37.s.vcf.gz 2>/dev/null; tabix -f fill37.s.vcf.gz
bcftools view -T ^have.tsv fill37.s.vcf.gz -Oz -o fill37.only.vcf.gz; tabix -f fill37.only.vcf.gz
bcftools concat -a target_pgx37.vcf.gz fill37.only.vcf.gz -Ou | bcftools annotate -x INFO,^FORMAT/GT -Ou | bcftools sort -Ou | bcftools annotate --rename-chrs rename.txt -Oz -o pgx37.chr.vcf.gz
tabix -f pgx37.chr.vcf.gz
echo "records to lift: $(bcftools view -H pgx37.chr.vcf.gz | wc -l) (calls $(bcftools view -H target_pgx37.vcf.gz | wc -l) + fills $(bcftools view -H fill37.only.vcf.gz | wc -l))"
java -Xmx8g -jar $PICARD LiftoverVcf -I pgx37.chr.vcf.gz -O pgx38.vcf.gz -CHAIN $CHAIN -REJECT rejected.vcf.gz -R $REF38 --RECOVER_SWAPPED_REF_ALT true --WARN_ON_MISSING_CONTIG true > liftover.log 2>&1
echo "lifted: $(bcftools view -H pgx38.vcf.gz | wc -l)  rejected: $(bcftools view -H rejected.vcf.gz | wc -l)  swapped: $(bcftools view -H -i 'INFO/SwappedAlleles=1' pgx38.vcf.gz | wc -l)"
bcftools view -H rejected.vcf.gz | cut -f1-5,7 > rejected.tsv
# outside call for CYP2D6 from Cyrius
CYP2D6=$(awk 'NR==2{print $2}' $WGS/06_pgx/cyrius/target.tsv); printf "CYP2D6\t%s\n" "$CYP2D6" > outside_calls.tsv
rm -rf prep out
python $PC/preprocessor/pharmcat_vcf_preprocessor -vcf pgx38.vcf.gz -refFna $REF38 -refVcf $PC/pharmcat_positions_3.4.0.vcf.bgz -o prep -bf target > preprocess.log 2>&1 || { tail -20 preprocess.log; exit 1; }
java -jar $PC/pharmcat-3.4.0-all.jar -vcf prep/target.preprocessed.vcf.bgz -po outside_calls.tsv -o out -reporterJson -matcherHtml > pharmcat.log 2>&1 || { tail -20 pharmcat.log; exit 1; }
echo "missing PGx positions after fill: $(grep -vc '^#' prep/target.missing_pgx_var.vcf || true)"
echo PHARMCAT_DONE
