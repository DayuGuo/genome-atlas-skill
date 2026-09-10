#!/usr/bin/env bash
# Genotype imputation of the WeGene array against 1000G phase 3 (GRCh37) with Beagle 5.4.
# Input: array genotypes re-coded on 1000G REF/ALT (all autosomal SNPs, including A/T & C/G sites: the
# harmonization audit found 0/505,634 strand flips, so the array is plus-strand). Output: results_v2/imputed/chrN.vcf.gz
set -euo pipefail
source "$(dirname "$0")/env.sh"
R=$REF/imp; O=$RESULTS2/imputed; mkdir -p $O
python3 - <<'PY'
import os, pandas as pd
W=os.environ["WORK"]; S=os.environ["SAMPLE"]
kg = pd.read_csv(f"{W}/results/ancestry/kg.overlap.pvar", sep="\t", comment="#", header=None, usecols=[0,1,2,3,4], names=["chrom","pos","id","ref","alt"], dtype={"chrom":str})
df = pd.read_parquet(f"{W}/data/genome.parquet"); df["chrom"]=df.chrom.astype(str)
m = kg.merge(df[["chrom","pos","a1","a2"]], on=["chrom","pos"])
ok = [set(g) <= {r,a} for g,r,a in zip(m.a1+m.a2, m.ref, m.alt)]
m = m[ok]
print("input sites:", len(m))
for ch in [str(i) for i in range(1,23)]:
    s = m[m.chrom==ch].sort_values("pos")
    with open(f"{W}/results_v2/imputed/target.chr{ch}.vcf","w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=%s>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % (ch, S))
        for r in s.itertuples():
            gt = "/".join("0" if x==r.ref else "1" for x in (r.a1, r.a2))
            f.write(f"{ch}\t{r.pos}\t.\t{r.ref}\t{r.alt}\t.\tPASS\t.\tGT\t{gt}\n")
PY
for ch in $(seq 1 22); do
  if [ -s $O/chr$ch.vcf.gz ]; then echo "chr$ch done"; continue; fi
  java -Xmx${BEAGLE_MEM:-32g} -jar $R/beagle.jar gt=$O/target.chr$ch.vcf ref=$R/ALL.chr$ch.vcf.gz map=$R/maps/plink.chr$ch.GRCh37.map out=$O/chr$ch nthreads=$THREADS gp=true > $O/chr$ch.beagle.log 2>&1
  echo "chr$ch finished $(date)"
done
