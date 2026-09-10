#!/usr/bin/env bash
# Reference data, about 60 GB. Everything is public and re-downloadable; nothing here is personal.
set -euo pipefail
source "$(dirname "$0")/../scripts/env.sh"
mkdir -p "$REF_DIR"/{b37,hg38,chain,annot,imp,aadr,archaic,ytree,prs}
cd "$REF_DIR"

# GRCh37 primary reference (the build every step below assumes)
[ -f b37/human_g1k_v37.fasta ] || {
  curl -sSL -o b37/human_g1k_v37.fasta.gz https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
  gzip -dc b37/human_g1k_v37.fasta.gz > b37/human_g1k_v37.fasta || true   # the archive has trailing bytes; harmless
  samtools faidx b37/human_g1k_v37.fasta; }

# GRCh38, only needed because PharmCAT is GRCh38-only
B38=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
for f in Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.fai Homo_sapiens_assembly38.dict; do
  [ -f "hg38/$f" ] || curl -sSL -o "hg38/$f" "$B38/$f"; done
for c in hg19ToHg38 hg38ToHg19; do
  [ -f "chain/$c.over.chain.gz" ] || curl -sSL -o "chain/$c.over.chain.gz" "https://hgdownload.soe.ucsc.edu/goldenPath/${c%%To*}/liftOver/$c.over.chain.gz"; done

# 1000 Genomes phase 3: plink2 bundle for PCA/PRS, per-chromosome VCFs for phasing and local ancestry
for e in pgen pvar.zst psam; do [ -f "all_phase3.$e" ] || curl -sSL -o "all_phase3.$e" "https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.$e?dl=1"; done
[ -f all_phase3.pvar ] && [ ! -f all_phase3.pvar.zst ] || "$PLINK2" --zst-decompress all_phase3.pvar.zst > all_phase3.pvar
K=https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf
for c in $(seq 1 22) X; do [ -f "imp/ALL.chr$c.vcf.gz" ] || curl -sSL -o "imp/ALL.chr$c.vcf.gz" "$K/chr$c.1kg.phase3.v5a.vcf.gz"; done
[ -f imp/beagle.jar ] || curl -sSL -o imp/beagle.jar https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar
[ -d imp/maps ] || { curl -sSL -o imp/plink.GRCh37.map.zip https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip; unzip -q -o imp/plink.GRCh37.map.zip -d imp/maps; }

# annotation
[ -f clinvar_grch37.vcf.gz ] || { curl -sSL -o clinvar_grch37.vcf.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
  curl -sSL -o clinvar_grch37.vcf.gz.tbi https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi; }
[ -f annot/Homo_sapiens.GRCh37.87.gff3.gz ] || curl -sSL -o annot/Homo_sapiens.GRCh37.87.gff3.gz https://ftp.ensembl.org/pub/grch37/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz
[ -f annot/Homo_sapiens.GRCh37.87.gtf.gz ] || curl -sSL -o annot/Homo_sapiens.GRCh37.87.gtf.gz https://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
[ -f annot/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz ] || curl -sSL -o annot/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

# Y tree (YFull) and the ybrowse SNP position index
[ -f ytree/current_tree.json ] || curl -sSL -o ytree/current_tree.json https://raw.githubusercontent.com/YFullTeam/YTree/master/current_tree.json
[ -f ytree/current_version.txt ] || curl -sSL -o ytree/current_version.txt https://raw.githubusercontent.com/YFullTeam/YTree/master/current_version.txt
[ -f ytree/snps_hg19.csv ] || curl -sSL -o ytree/snps_hg19.csv http://ybrowse.org/gbrowse2/gff/snps_hg19.csv

# optional panels: ancient DNA (AADR Human Origins) and archaic introgression (Sprime, Browning 2018)
if [ "${WITH_AADR:-1}" = "1" ]; then
  DV=https://dataverse.harvard.edu/api/access/datafile
  for id_name in "13994526:v66.p1_HO.aadr.patch.PUB.ind" "13994527:v66.p1_HO.aadr.patch.PUB.snp" \
                 "13994528:v66.p1_HO.aadr.PUB.anno" "13994808:v66.p1_HO.aadr.patch.PUB.geno"; do
    id=${id_name%%:*}; n=${id_name#*:}; [ -f "aadr/$n" ] || curl -sSL -o "aadr/$n" "$DV/$id"; done
fi
if [ "${WITH_SPRIME:-1}" = "1" ]; then
  python3 - "$REF_DIR/archaic" <<'PY'
import json,sys,urllib.request,pathlib,subprocess
d=pathlib.Path(sys.argv[1]); d.mkdir(parents=True,exist_ok=True)
meta=json.load(urllib.request.urlopen("https://data.mendeley.com/public-api/datasets/y7hyt83vxr/files?folder_id=root&version=1"))
want={f"{p}_sprime_results.tar.gz" for p in ("CHB","CHS","JPT","CEU","GBR","CDX","KHV","BEB","GIH")}
for f in meta:
    if f["filename"] in want and not (d/f["filename"]).exists():
        urllib.request.urlretrieve(f["content_details"]["download_url"], d/f["filename"])
        subprocess.run(["tar","xzf",str(d/f["filename"]),"-C",str(d)],check=True)
PY
fi
echo "references ready under $REF_DIR"
