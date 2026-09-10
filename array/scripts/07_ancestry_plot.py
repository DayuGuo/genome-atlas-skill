#!/usr/bin/env python3
"""Summarise the PCA projection: global placement and East-Asian-only PCA, with plots."""
import pathlib, subprocess, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from config import *
ROOT = pathlib.Path(__file__).resolve().parents[1]; W = RESULTS/"ancestry"; P = pathlib.Path(PLINK2)
PAL = {"EAS":"#2a78d6","EUR":"#eb6834","AFR":"#1baf7a","SAS":"#eda100","AMR":"#e87ba4"}
MRK = {"EAS":"o","EUR":"s","AFR":"^","SAS":"D","AMR":"v"}

def load(prefix):
    kg = pd.read_csv(W/f"{prefix}kg.common.proj.sscore" if prefix=="" else W/f"{prefix}.kg.proj.sscore", sep="\t")
    me = pd.read_csv(W/f"{prefix}sample.proj.sscore" if prefix=="" else W/f"{prefix}.sample.proj.sscore", sep="\t")
    return kg, me

def report(kg, me, label, fh):
    pcs = [c for c in kg.columns if c.startswith("PC")]
    cent = kg.groupby("Population")[pcs].mean()
    sd = kg.groupby("Population")[pcs].std()
    v = me[pcs].iloc[0]
    # Mahalanobis-like (diagonal) distance to each population using first 4 PCs
    use = pcs[:4]
    d = (((cent[use] - v[use]) / sd[use])**2).sum(axis=1).pow(0.5).sort_values()
    fh.write(f"\n== {label} ==\nDayu PCs: " + ", ".join(f"{c}={v[c]:.4f}" for c in use) + "\n")
    fh.write("Nearest populations (standardised distance on PC1-4, lower = closer):\n")
    for pop, val in d.head(8).items():
        fh.write(f"  {pop:4s} ({kg[kg.Population==pop].SuperPop.iloc[0]})  {val:6.2f}\n")
    return d

out = open(RESULTS/"ancestry_summary.txt", "w")
out.write(f"Ancestry PCA: {SAMPLE} projected onto 1000 Genomes phase 3 PCs (LD-pruned autosomal SNPs)\n")
kg, me = load("")
d_glob = report(kg, me, "Global (26 populations)", out)

# ---- East Asian only PCA
eas = kg[kg.SuperPop==REF_SUPERPOP][["#IID"]]; eas.to_csv(W/"eas.ids", sep="\t", index=False, header=False)
pfx = "eas"
subprocess.run([P,"--pfile",W/"kg.common","--keep",W/"eas.ids","--maf","0.01","--indep-pairwise","200","50","0.2","--out",W/f"{pfx}.prune","--threads", THREADS],check=True,capture_output=True)
subprocess.run([P,"--pfile",W/"kg.common","--keep",W/"eas.ids","--extract",W/f"{pfx}.prune.prune.in","--freq","--pca","10","allele-wts","--out",W/f"{pfx}.pca","--threads", THREADS],check=True,capture_output=True)
for s,o in [("kg.common",f"{pfx}.kg"),("sample",f"{pfx}.sample")]:
    keep = ["--keep", W/"eas.ids"] if s=="kg.common" else []
    subprocess.run([P,"--pfile",W/s,*keep,"--extract",W/f"{pfx}.prune.prune.in","--read-freq",W/f"{pfx}.pca.afreq",
                    "--score",W/f"{pfx}.pca.eigenvec.allele","2","5","header-read","no-mean-imputation","variance-standardize",
                    "--score-col-nums","6-15","--out",W/f"{o}.proj","--threads", THREADS],check=True,capture_output=True)
kge, mee = load(pfx)
d_eas = report(kge, mee, "East Asia only (CHB=北京汉族, CHS=南方汉族, CDX=西双版纳傣族, JPT=日本, KHV=越南京族)", out)
out.close(); print(open(RESULTS/"ancestry_summary.txt").read())

# ---- plots
fig, axes = plt.subplots(1, 2, figsize=(13, 6), facecolor="#fcfcfb")
ax = axes[0]
for sp, g in kg.groupby("SuperPop"):
    ax.scatter(g.PC1_AVG, g.PC2_AVG, s=14, c=PAL[sp], marker=MRK[sp], alpha=0.55, linewidths=0, label=sp)
ax.scatter(me.PC1_AVG, me.PC2_AVG, s=260, c="#0b0b0b", marker="*", edgecolors="#fcfcfb", linewidths=1.5, zorder=5, label=SAMPLE)
ax.set_title("Global: 1000 Genomes super-populations", fontsize=12, loc="left"); ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
ax.legend(frameon=False, fontsize=9)
ax = axes[1]
EPAL = {"CHB":"#2a78d6","CHS":"#eb6834","CDX":"#1baf7a","JPT":"#eda100","KHV":"#e87ba4"}
EMRK = {"CHB":"o","CHS":"s","CDX":"^","JPT":"D","KHV":"v"}
for pop, g in kge.groupby("Population"):
    ax.scatter(g.PC1_AVG, g.PC2_AVG, s=22, c=EPAL[pop], marker=EMRK[pop], alpha=0.65, linewidths=0, label=pop)
ax.scatter(mee.PC1_AVG, mee.PC2_AVG, s=260, c="#0b0b0b", marker="*", edgecolors="#fcfcfb", linewidths=1.5, zorder=5, label=SAMPLE)
ax.set_title("East Asia only: CHB / CHS / CDX / JPT / KHV", fontsize=12, loc="left"); ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
ax.legend(frameon=False, fontsize=9)
for a in axes:
    a.set_facecolor("#fcfcfb"); a.grid(True, color="#e6e5e0", linewidth=0.6); a.tick_params(colors="#52514e")
    for s in ["top","right"]: a.spines[s].set_visible(False)
    for s in ["left","bottom"]: a.spines[s].set_color("#c3c2b7")
plt.tight_layout(); plt.savefig(RESULTS/"ancestry_pca.png", dpi=150); print("saved results/ancestry_pca.png")
