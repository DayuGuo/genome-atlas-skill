#!/usr/bin/env python3
"""Fine-scale East Asian ancestry with HGDP (Bergström 2020 WGS callset, GRCh38).
1) liftover the sample's array positions hg19 -> hg38 (pyliftover, UCSC chain)
2) pull HGDP biallelic SNPs at those positions, build the sample VCF with HGDP REF/ALT
3) PCA on HGDP East Asian samples (+Uygur excluded), project the sample
4) nearest populations + NNLS mixture on population centroids (rough "calculator" style)
"""
import pathlib, subprocess, numpy as np, pandas as pd
from pyliftover import LiftOver
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from config import *
ROOT = pathlib.Path(__file__).resolve().parents[1]; R = REF; W = RESULTS/"hgdp"; W.mkdir(exist_ok=True)
P = str(pathlib.Path(PLINK2))
import os
def run(*a):
    out = str(a[list(a).index("--out")+1])
    if os.environ.get("SKIP_DONE") and any(pathlib.Path(out+ext).exists() for ext in (".pgen",".prune.in",".eigenvec",".sscore")): return
    subprocess.run([P, *map(str, a), "--threads", THREADS], check=True, capture_output=True)

# 1) liftover
df = pd.read_parquet(DATA/"genome.parquet"); df["chrom"] = df.chrom.astype(str)
d = df[(~df.is_missing) & df.chrom.isin([str(i) for i in range(1,23)]) & df.a1.isin(list("ACGT")) & df.a2.isin(list("ACGT"))].copy()
lo = LiftOver(str(R/"chain"/"hg19ToHg38.over.chain.gz"))
new = []
for c, p in zip(d.chrom, d.pos):
    r = lo.convert_coordinate("chr"+c, p-1)
    new.append((r[0][0][3:], r[0][1]+1) if r and r[0][0] == "chr"+c else (None, None))
d["chrom38"], d["pos38"] = zip(*new); d = d.dropna(subset=["pos38"]); d["pos38"] = d.pos38.astype(int)
print("lifted to hg38:", len(d))
d[["chrom38","pos38"]].assign(e=d.pos38, i=[f"r{i}" for i in range(len(d))]).to_csv(W/"sample.hg38.range", sep="\t", header=False, index=False)

# 2) HGDP extraction (East Asia region only) and the sample VCF
psam = pd.read_csv(R/"hgdp_all.psam", sep="\t"); eas = psam[psam.region==HGDP_REGION]
eas[["#IID"]].to_csv(W/"eas.ids", sep="\t", index=False, header=False)
run("--pfile", R/"hgdp_all", "--keep", W/"eas.ids", "--extract", "range", W/"sample.hg38.range", "--snps-only", "just-acgt",
    "--max-alleles", "2", "--autosome", "--rm-dup", "exclude-all", "--set-all-var-ids", "@:#", "--make-pgen", "--out", W/"hgdp.eas")
kg = pd.read_csv(W/"hgdp.eas.pvar", sep="\t", comment="#", header=None, usecols=[0,1,2,3,4], names=["chrom38","pos38","id","ref","alt"], dtype={"chrom38":str})
m = kg.merge(d[["chrom38","pos38","a1","a2"]], on=["chrom38","pos38"])
amb = ((m.ref=="A")&(m.alt=="T"))|((m.ref=="T")&(m.alt=="A"))|((m.ref=="C")&(m.alt=="G"))|((m.ref=="G")&(m.alt=="C"))
m = m[~amb]
COMP = str.maketrans("ACGT","TGCA")
def code(r):
    g = [r.a1, r.a2]; al = {r.ref, r.alt}
    if not set(g) <= al:
        g = [x.translate(COMP) for x in g]
        if not set(g) <= al: return None
    return "/".join("0" if x==r.ref else "1" for x in g)
m["gt_"] = [code(r) for r in m.itertuples()]; bad = m.gt_.isna().sum(); m = m[m.gt_.notna()]
print(f"HGDP overlap {len(kg)}, ambiguous dropped {int(amb.sum())}, mismatch dropped {bad}, used {len(m)}")
with open(W/"sample.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n" + "".join(f"##contig=<ID={c}>\n" for c in map(str, range(1,23))))
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDayu\n')
    for r in m.itertuples(): f.write(f"{r.chrom38}\t{r.pos38}\t{r.id}\t{r.ref}\t{r.alt}\t.\t.\t.\tGT\t{r.gt_}\n")
m.id.to_csv(W/"common.ids", index=False, header=False)
run("--vcf", W/"sample.vcf", "--make-pgen", "--out", W/"sample")
run("--pfile", W/"hgdp.eas", "--extract", W/"common.ids", "--make-pgen", "--out", W/"hgdp.common")

# 3) PCA + projection
run("--pfile", W/"hgdp.common", "--maf", "0.02", "--indep-pairwise", "200", "50", "0.2", "--out", W/"prune")
run("--pfile", W/"hgdp.common", "--extract", W/"prune.prune.in", "--freq", "--pca", "10", "allele-wts", "--out", W/"pca")
for s in ["hgdp.common", "sample"]:
    run("--pfile", W/s, "--extract", W/"prune.prune.in", "--read-freq", W/"pca.afreq", "--score", W/"pca.eigenvec.allele", "2", "5",
        "header-read", "no-mean-imputation", "variance-standardize", "--score-col-nums", "6-15", "--out", W/f"{s}.proj")
ref = pd.read_csv(W/"hgdp.common.proj.sscore", sep="\t").rename(columns={"#IID":"IID"})
me = pd.read_csv(W/"sample.proj.sscore", sep="\t")
pcs = [c for c in ref.columns if c.startswith("PC")]; v = me[pcs].iloc[0].to_numpy()
ev = np.loadtxt(W/"pca.eigenval"); pve = ev/ev.sum()*100

# 4) nearest populations: Euclidean distance to centroids on PC1-6 (weighted by sqrt eigenvalue),
#    and population make-up of the 15 nearest individuals (kNN)
k = 6; wts = np.sqrt(ev[:k])
X = ref[pcs[:k]].to_numpy()*wts; b = v[:k]*wts
cent = pd.DataFrame(X, index=ref.population).groupby(level=0).mean()
dist = pd.Series(np.linalg.norm(cent.to_numpy()-b, axis=1), index=cent.index).sort_values()
ind = np.linalg.norm(X-b, axis=1); order = np.argsort(ind)[:15]
knn = ref.population.iloc[order].value_counts()
with open(RESULTS/"hgdp_ancestry_summary.txt", "w") as f:
    f.write(f"HGDP East Asia PCA ({len(ref)} samples, {sum(1 for _ in open(W/'prune.prune.in'))} LD-pruned SNPs); PC variance %: " + ", ".join(f"{x:.1f}" for x in pve[:6]) + "\n")
    f.write("{SAMPLE} PC1-4: " + ", ".join(f"{x:.4f}" for x in v[:4]) + "\n\nNearest population centroids (eigenvalue-weighted Euclidean distance, PC1-6):\n")
    for p_, x in dist.head(8).items(): f.write(f"  {p_:12s} {x:7.4f}\n")
    f.write("\n15 nearest individuals by population:\n")
    for p_, n in knn.items(): f.write(f"  {p_:12s} {n}\n")
print(open(RESULTS/"hgdp_ancestry_summary.txt").read())
ref.to_csv(W/"hgdp.eas.pcs.tsv", sep="\t", index=False)

# plot PC1/2 and PC3/4 with labels at centroids
fig, axes = plt.subplots(1, 2, figsize=(14, 6.5), facecolor="#fcfcfb")
pops = sorted(ref.population.unique())
markers = ["o","s","^","D","v","P","X","<",">","*","h","p","d","8","H","o","s","^"]
for ax, (x, y) in zip(axes, [("PC1_AVG","PC2_AVG"),("PC3_AVG","PC4_AVG")]):
    for i, p in enumerate(pops):
        g = ref[ref.population==p]; col = "#2a78d6" if p in ("Han","NorthernHan") else ("#eb6834" if p in ("Japanese",) else "#8a8984")
        ax.scatter(g[x], g[y], s=26, c=col, marker=markers[i%len(markers)], alpha=0.7, linewidths=0)
        ax.annotate(p, (g[x].mean(), g[y].mean()), fontsize=8, color="#52514e", ha="center", va="bottom", xytext=(0,4), textcoords="offset points")
    ax.scatter(me[x], me[y], s=280, c="#0b0b0b", marker="*", edgecolors="#fcfcfb", linewidths=1.5, zorder=5)
    ax.annotate(SAMPLE, (me[x].iloc[0], me[y].iloc[0]), fontsize=10, fontweight="bold", xytext=(6,6), textcoords="offset points")
    i1, i2 = int(x[2])-1, int(y[2])-1
    ax.set_xlabel(f"{x[:3]} ({pve[i1]:.1f}%)"); ax.set_ylabel(f"{y[:3]} ({pve[i2]:.1f}%)")
    ax.set_facecolor("#fcfcfb"); ax.grid(True, color="#e6e5e0", linewidth=0.6)
    for s_ in ["top","right"]: ax.spines[s_].set_visible(False)
axes[0].set_title("HGDP East Asia: PC1 vs PC2 (blue = Han / NorthernHan, orange = Japanese)", fontsize=11, loc="left")
axes[1].set_title("PC3 vs PC4", fontsize=11, loc="left")
plt.tight_layout(); plt.savefig(RESULTS/"hgdp_pca.png", dpi=150); print("saved results/hgdp_pca.png")
