#!/usr/bin/env bash
# Ancestry PCA: the sample + 1000 Genomes phase 3 (GRCh37), PCA on 1000G with the sample projected.
# Steps: (1) pull 1000G biallelic SNPs at the sample's autosomal array positions, (2) build a
# the sample VCF using 1000G REF/ALT (strand-flip fix, drop A/T & C/G ambiguous sites),
# (3) merge, LD-prune, PCA. One extra sample among 2504 barely perturbs the PCs.
set -euo pipefail
source "$(dirname "$0")/env.sh"
P=$PLINK2; R=$REF; W=$RESULTS/ancestry; mkdir -p $W

python3 - <<'PY'
import os, pandas as pd
W=os.environ["WORK"]; S=os.environ["SAMPLE"]
df = pd.read_parquet(f"{W}/data/genome.parquet")
d = df[(~df.is_missing) & df.chrom.isin([str(i) for i in range(1,23)]) & df.a1.isin(list("ACGT")) & df.a2.isin(list("ACGT"))]
d[["chrom","pos"]].assign(end=d.pos, id=["r%d"%i for i in range(len(d))]).to_csv(f"{W}/results/ancestry/sample.autosomal.range", sep="\t", header=False, index=False)
print("autosomal SNPs:", len(d))
PY
$P --pfile $R/all_phase3 --extract range $W/sample.autosomal.range --snps-only just-acgt --max-alleles 2 \
   --autosome --rm-dup exclude-all --set-all-var-ids '@:#' --make-pgen --out $W/kg.overlap --threads $THREADS >/dev/null

python3 - <<'PY'
import pandas as pd
COMP = str.maketrans("ACGT","TGCA")
import os
W=os.environ["WORK"]; S=os.environ["SAMPLE"]
kg = pd.read_csv(f"{W}/results/ancestry/kg.overlap.pvar", sep="\t", comment="#", header=None, usecols=[0,1,2,3,4], names=["chrom","pos","id","ref","alt"], dtype={"chrom":str})
kg = kg[kg.ref.str.len().eq(1) & kg.alt.str.len().eq(1)]
df = pd.read_parquet(f"{W}/data/genome.parquet"); df["chrom"]=df.chrom.astype(str)
m = kg.merge(df[["chrom","pos","genotype","a1","a2"]], on=["chrom","pos"])
amb = ((m.ref=="A")&(m.alt=="T"))|((m.ref=="T")&(m.alt=="A"))|((m.ref=="C")&(m.alt=="G"))|((m.ref=="G")&(m.alt=="C"))
m = m[~amb].copy()

def code(r):
    al = {r.ref, r.alt}
    g = [r.a1, r.a2]
    if set(g) <= al: pass
    elif set(x.translate(COMP) for x in g) <= al: g = [x.translate(COMP) for x in g]
    else: return None
    return "/".join("0" if x==r.ref else "1" for x in g)
m["gt_"] = [code(r) for r in m.itertuples()]
n_bad = m.gt_.isna().sum(); m = m[m.gt_.notna()]
with open(f"{W}/results/ancestry/sample.vcf","w") as f:
    f.write("##fileformat=VCFv4.2\n" + "".join(f"##contig=<ID={c}>\n" for c in [str(i) for i in range(1,23)]))
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+S+'\n')
    for r in m.itertuples(): f.write(f"{r.chrom}\t{r.pos}\t{r.id}\t{r.ref}\t{r.alt}\t.\t.\t.\tGT\t{r.gt_}\n")
m.id.to_csv(f"{W}/results/ancestry/common.ids", index=False, header=False)
print(f"1000G overlap: {len(kg)}  ambiguous dropped: {int(amb.sum())}  allele-mismatch dropped: {n_bad}  used: {len(m)}")
PY
$P --vcf $W/sample.vcf --make-pgen --out $W/sample --threads $THREADS >/dev/null
printf "#IID\tSEX\tSuperPop\tPopulation\n%s\t%s\tME\t%s\n" "$SAMPLE" "$([ "${SEX:-male}" = male ] && echo 1 || echo 2)" "$SAMPLE" > $W/sample.psam
$P --pfile $W/kg.overlap --extract $W/common.ids --make-pgen --out $W/kg.common --threads $THREADS >/dev/null
$P --pfile $W/kg.common --maf 0.01 --indep-pairwise 200 50 0.2 --out $W/prune --threads $THREADS >/dev/null
$P --pfile $W/kg.common --extract $W/prune.prune.in --freq --pca 10 allele-wts --out $W/kg.pca --threads $THREADS >/dev/null
for s in kg.common sample; do
  $P --pfile $W/$s --extract $W/prune.prune.in --read-freq $W/kg.pca.afreq \
     --score $W/kg.pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
     --score-col-nums 6-15 --out $W/$s.proj --threads $THREADS >/dev/null
done
echo "pruned SNPs: $(wc -l < $W/prune.prune.in)"; head -2 $W/sample.proj.sscore
