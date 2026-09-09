#!/usr/bin/env python3
"""Round 3 'fun' analyses on the existing 1000G-aligned data.
 A) per-chromosome north/south index on the CHB-CHS axis (East Asian PCA, PC1 per chromosome)
 B) X-chromosome ancestry (a male's X is entirely maternal): PCA of 1000G X, project the sample
 C) rarest alleles the sample carries (1000G EAS frequency < 0.5%), annotated with Ensembl VEP
 D) reference-genome similarity: share of sites homozygous for the GRCh37 reference allele,
    compared with every 1000G individual on the same sites
Outputs: results/fun_*.tsv / .txt
"""
import json, pathlib, subprocess, requests, numpy as np, pandas as pd
from config import *
ROOT = pathlib.Path(__file__).resolve().parents[1]; A = RESULTS/"ancestry"; W = RESULTS/"fun"; W.mkdir(exist_ok=True)
P = str(pathlib.Path(PLINK2)); R = REF
def run(*a):
    r = subprocess.run([P, *map(str, a), "--threads", THREADS], capture_output=True, text=True)
    if r.returncode: raise SystemExit(r.stdout[-2000:] + r.stderr[-2000:])

# ---------- A) per-chromosome N/S index ----------
rows = []
import os
for ch in ([] if os.environ.get("SKIP_A") else range(1, 23)):
    for s in ["kg.common", "sample"]:
        keep = ["--keep", A/"eas.ids"] if s == "kg.common" else []
        run("--pfile", A/s, *keep, "--chr", str(ch), "--extract", A/"eas.prune.prune.in", "--read-freq", A/"eas.pca.afreq",
            "--score", A/"eas.pca.eigenvec.allele", "2", "5", "header-read", "no-mean-imputation", "variance-standardize",
            "--score-col-nums", "6", "--out", W/f"chr{ch}.{s}")
    k = pd.read_csv(W/f"chr{ch}.kg.common.sscore", sep="\t"); d = pd.read_csv(W/f"chr{ch}.sample.sscore", sep="\t")
    chb, chs = k[k.Population==AXIS[1]].PC1_AVG, k[k.Population==AXIS[0]].PC1_AVG
    x = d.PC1_AVG.iloc[0]
    idx = (x - chs.mean()) / (chb.mean() - chs.mean())          # 1 = CHB centroid, 0 = CHS centroid
    pooled_sd = np.sqrt((chb.var() + chs.var()) / 2)
    rows.append(dict(chrom=ch, n_snps=int(d.ALLELE_CT.iloc[0]//2), sample_pc1=x, chb_mean=chb.mean(), chs_mean=chs.mean(),
                     ns_index=round(idx, 3), sep_sd=round((chb.mean()-chs.mean())/pooled_sd, 2)))
ns = pd.DataFrame(rows) if rows else pd.read_csv(RESULTS/"fun_chrom_north_south.tsv", sep="\t")
if rows: ns.to_csv(RESULTS/"fun_chrom_north_south.tsv", sep="\t", index=False)
print("A) per-chromosome north/south index (1 = CHB, 0 = CHS):"); print(ns[["chrom","n_snps","ns_index","sep_sd"]].to_string(index=False))
print("genome-wide mean index:", round(ns.ns_index.mean(), 2))

# ---------- B) X-chromosome ancestry ----------
df = pd.read_parquet(DATA/"genome.parquet"); df["chrom"] = df.chrom.astype(str)
x = df[(df.chrom=="X") & (~df.is_missing) & df.a1.isin(list("ACGT"))]
x[["chrom","pos"]].assign(e=x.pos, i=[f"x{i}" for i in range(len(x))]).to_csv(W/"x.range", sep="\t", header=False, index=False)
run("--pfile", R/"all_phase3", "--chr", "X", "--extract", "range", W/"x.range", "--snps-only", "just-acgt", "--max-alleles", "2",
    "--rm-dup", "exclude-all", "--set-all-var-ids", "@:#", "--make-pgen", "--out", W/"kgX")
kx = pd.read_csv(W/"kgX.pvar", sep="\t", comment="#", header=None, usecols=[0,1,2,3,4], names=["chrom","pos","id","ref","alt"], dtype={"chrom":str})
m = kx.merge(x[["chrom","pos","a1"]], on=["chrom","pos"])
amb = ((m.ref=="A")&(m.alt=="T"))|((m.ref=="T")&(m.alt=="A"))|((m.ref=="C")&(m.alt=="G"))|((m.ref=="G")&(m.alt=="C")); m = m[~amb]
COMP = str.maketrans("ACGT","TGCA")
def code(r):
    a = r.a1 if r.a1 in (r.ref, r.alt) else r.a1.translate(COMP)
    return None if a not in (r.ref, r.alt) else ("0/0" if a == r.ref else "1/1")
m["gt_"] = [code(r) for r in m.itertuples()]; m = m[m.gt_.notna()]
with open(W/"sampleX.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n##contig=<ID=X>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDayu\n")
    for r in m.itertuples(): f.write(f"X\t{r.pos}\t{r.id}\t{r.ref}\t{r.alt}\t.\t.\t.\tGT\t{r.gt_}\n")
m.id.to_csv(W/"x.ids", index=False, header=False)
open(W/"sampleX.sex.txt","w").write("#IID\tSEX\nDayu\t1\n")
run("--vcf", W/"sampleX.vcf", "--psam", W/"sampleX.sex.txt", "--make-pgen", "--out", W/"sampleX")
# plink2 PCA is autosome-only: relabel X as a pseudo-autosome "23"->"1" via VCF round-trip (haploid male calls become homozygous)
def relabel(src, dst):
    run("--pfile", W/src, "--export", "vcf", "--out", W/f"{src}.tmp")
    subprocess.run(f"sed -e 's/^X\t/1\t/' -e 's/ID=X/ID=1/' {W}/{src}.tmp.vcf | sed -E 's#\t([01])(\t|$)#\t\\1/\\1\\2#g' > {W}/{dst}.vcf", shell=True, check=True)
    run("--vcf", W/f"{dst}.vcf", "--make-pgen", "--out", W/dst)
run("--pfile", W/"kgX", "--extract", W/"x.ids", "--make-pgen", "--out", W/"kgX.raw")
relabel("kgX.raw", "kgX.common"); relabel("sampleX", "sampleXa")
import shutil; shutil.copy(W/"kgX.raw.psam", W/"kgX.common.psam")
run("--pfile", W/"kgX.common", "--maf", "0.01", "--indep-pairwise", "200", "50", "0.2", "--out", W/"xprune")
run("--pfile", W/"kgX.common", "--extract", W/"xprune.prune.in", "--freq", "--pca", "6", "allele-wts", "--out", W/"xpca")
for s in ["kgX.common", "sampleXa"]:
    run("--pfile", W/s, "--extract", W/"xprune.prune.in", "--read-freq", W/"xpca.afreq", "--score", W/"xpca.eigenvec.allele", "2", "5",
        "header-read", "no-mean-imputation", "variance-standardize", "--score-col-nums", "6-11", "--out", W/f"{s}.proj")
kxp = pd.read_csv(W/"kgX.common.proj.sscore", sep="\t"); dxp = pd.read_csv(W/"sampleXa.proj.sscore", sep="\t")
pcs = ["PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG"]; v = dxp[pcs].iloc[0]
cent = kxp.groupby("Population")[pcs].mean(); dist = np.sqrt(((cent - v)**2).sum(axis=1)).sort_values()
sp = kxp.groupby("SuperPop")[pcs].mean(); dsp = np.sqrt(((sp - v)**2).sum(axis=1)).sort_values()
xs = f"B) X chromosome ({sum(1 for _ in open(W/'xprune.prune.in'))} LD-pruned SNPs, maternal lineage only)\n  nearest super-populations: " + ", ".join(f"{k} {x:.3f}" for k, x in dsp.items()) + "\n  nearest populations: " + ", ".join(f"{k} {x:.3f}" for k, x in dist.head(6).items())
print(xs); open(W/"x_summary.txt","w").write(xs)
kxp[["#IID","SuperPop","Population"]+pcs].to_csv(W/"x_pcs.tsv", sep="\t", index=False); dxp[pcs].to_csv(W/"x_sample.tsv", sep="\t", index=False)

# ---------- C) rarest alleles ----------
run("--pfile", A/"kg.common", "--keep", A/"eas.ids", "--freq", "--out", W/"eas")
run("--pfile", A/"kg.common", "--freq", "--out", W/"all")
fe = pd.read_csv(W/"eas.afreq", sep="\t").rename(columns={"#CHROM":"chrom","ALT_FREQS":"f_eas"})[["ID","REF","ALT","f_eas"]]
fa = pd.read_csv(W/"all.afreq", sep="\t").rename(columns={"ALT_FREQS":"f_all"})[["ID","f_all"]]
gt = pd.read_csv(A/"sample.vcf", sep="\t", comment="#", header=None, usecols=[0,1,2,3,4,9], names=["chrom","pos","ID","ref","alt","geno"], dtype={"chrom":str})
g = gt.merge(fe, on="ID").merge(fa, on="ID")
g["n_alt"] = g.geno.str.count("1")
carried_rare = g[(g.n_alt > 0) & (g.f_eas < 0.005)].copy(); carried_rare["allele"] = carried_rare.alt; carried_rare["freq_eas"] = carried_rare.f_eas; carried_rare["freq_all"] = carried_rare.f_all
ref_rare = g[(g.n_alt < 2) & (g.f_eas > 0.995)].copy(); ref_rare["allele"] = ref_rare.ref; ref_rare["freq_eas"] = 1-ref_rare.f_eas; ref_rare["freq_all"] = 1-ref_rare.f_all
rare = pd.concat([carried_rare, ref_rare]).sort_values("freq_eas")
rs = df[df.chrom.isin([str(i) for i in range(1,23)])][["chrom","pos","rsid"]]
rare = rare.merge(rs, on=["chrom","pos"], how="left")
print(f"C) alleles carried with EAS freq < 0.5%: {len(rare)}  (of which absent in all 504 EAS: {(rare.freq_eas==0).sum()})")
# annotate top 40 with Ensembl VEP (GRCh37)
top = rare.head(40).copy()
regions = [f"{r.chrom} {r.pos} . {r.ref} {r.alt}" if r.allele == r.alt else f"{r.chrom} {r.pos} . {r.alt} {r.ref}" for r in top.itertuples()]
ann = {}
try:
    rr = requests.post("https://grch37.rest.ensembl.org/vep/homo_sapiens/region", headers={"Content-Type":"application/json","Accept":"application/json"},
                       json={"variants": regions}, timeout=180); rr.raise_for_status()
    for v in rr.json():
        key = v["input"]; tc = v.get("transcript_consequences") or []
        genes = sorted({t.get("gene_symbol","") for t in tc if t.get("gene_symbol")})
        ann[key] = (";".join(genes[:3]), v.get("most_severe_consequence",""), ";".join(sorted({c.get("clin_sig_allele","") for t in tc for c in []})))
except Exception as e:
    print("VEP unavailable:", e)
top["gene"] = [ann.get(k, ("","",""))[0] for k in regions]; top["consequence"] = [ann.get(k, ("","",""))[1] for k in regions]
out = top[["rsid","chrom","pos","allele","geno","freq_eas","freq_all","gene","consequence"]]
out.to_csv(RESULTS/"fun_rarest_variants.tsv", sep="\t", index=False)
rare[["rsid","chrom","pos","allele","geno","freq_eas","freq_all"]].to_csv(W/"rare_all.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); print(out.head(25).to_string(index=False))

# ---------- D) reference similarity ----------
run("--pfile", A/"kg.common", "--sample-counts", "cols=homref,homalt,het", "--out", W/"kgcounts")
kc = pd.read_csv(W/"kgcounts.scount", sep="\t")
n_sites = len(g); sample_homref = int((g.n_alt == 0).sum()); sample_het = int((g.n_alt == 1).sum())
kc["homref_frac"] = kc.HOM_REF_CT / (kc.HOM_REF_CT + kc.HOM_ALT_CT + kc.HET_CT)
kc = kc.rename(columns={"#IID":"IID"}).merge(pd.read_csv(R/"all_phase3.psam", sep="\t").rename(columns={"#IID":"IID"})[["IID","SuperPop","Population"]], on="IID")
byp = kc.groupby("SuperPop").homref_frac.agg(["mean","min","max"]).round(4)
ds = f"D) reference similarity on {n_sites:,} shared autosomal SNPs\n  {SAMPLE} homozygous-reference share: {sample_homref/n_sites:.4f}  (het {sample_het/n_sites:.4f})\n" + byp.to_string()
print(ds)
kc[["IID","SuperPop","Population","homref_frac"]].to_csv(W/"kg_homref.tsv", sep="\t", index=False)
with open(RESULTS/"fun_summary.txt", "w") as f:
    f.write("A) per-chromosome north/south index (1 = CHB centroid, 0 = CHS centroid)\n" + ns[["chrom","n_snps","ns_index","sep_sd"]].to_string(index=False) + f"\n  genome-wide mean: {ns.ns_index.mean():.2f}\n\n" + xs + "\n\n" + f"C) rare alleles (EAS freq < 0.5%): {len(rare)}; absent in EAS: {(rare.freq_eas==0).sum()}\n\n" + ds + "\n")
json.dump({"sample_homref": sample_homref/n_sites, "sample_het": sample_het/n_sites, "n_sites": n_sites}, open(W/"refsim.json","w"))
