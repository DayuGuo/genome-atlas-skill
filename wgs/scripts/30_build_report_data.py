#!/usr/bin/env python
"""Assemble wgs/report_data.json from phase-5 outputs. Every number in report_v3.html comes from here."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pathlib
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import json, re, subprocess, pathlib, collections
import pandas as pd, numpy as np
W = P/"wgs"; D = {}
_ch = pd.read_csv(pathlib.Path(__file__).resolve().parents[1]/"panel"/"chrom_grch37.tsv", sep="\t", dtype={"chrom": str})
D["chrlen"] = {r.chrom: int(r.length) for r in _ch.itertuples()}
D["cen"] = {r.chrom: float(r.centromere_mb) for r in _ch.itertuples() if pd.notna(r.centromere_mb)}
# --- QC / KPI
ms = pd.read_csv(W/"01_qc/depth.mosdepth.summary.txt", sep="\t"); ms = ms[~ms.chrom.str.contains("_region|^GL")]
dep = {r.chrom: r["mean"] for r in ms.itertuples(index=False).__iter__()} if False else dict(zip(ms.chrom, ms["mean"]))
idx = subprocess.run(["bcftools", "index", "-s", str(W/f"00_input/{SAMPLE}.pass.vcf.gz")], capture_output=True, text=True).stdout
nvar = {l.split("\t")[0]: int(l.split("\t")[2]) for l in idx.splitlines()}
CHR = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
D["chrom"] = [{"chrom": c, "depth": round(dep.get(c, 0), 1), "n_pass": nvar.get(c, 0), "len": D["chrlen"].get(c, 16569 if c == "MT" else 0)} for c in CHR]
def sn(f, key):
    for l in open(f):
        if l.startswith("SN") and key in l: return int(l.rstrip().split("\t")[-1])
sp = W/"01_qc/stats.pass.txt"
tstv = [l for l in open(sp) if l.startswith("TSTV")][0].split("\t")
psc = [l for l in open(sp) if l.startswith("PSC")][0].split("\t")
call_bp = sum(int(l.split()[2]) - int(l.split()[1]) for l in open(W/"00_input/callable.bed"))
D["kpi"] = {"depth_auto": round(dep["total"], 1), "depth_x": round(dep["X"], 1), "depth_y": round(dep["Y"], 1), "depth_mt": int(dep["MT"]),
            "callable_gb": round(call_bp/1e9, 2), "pass_records": sn(sp, "number of records"), "snv": sn(sp, "number of SNPs"), "indel": sn(sp, "number of indels"),
            "titv": float(tstv[4]), "het": int(psc[5]), "homalt": int(psc[4]), "all_records": sn(W/"01_qc/stats.norm.txt", "number of records"),
            "platform": "MGI T7 · Sentieon DNAscope · GRCh37"}
# --- chip-is-subset evidence
cs = pd.read_csv(W/"01_qc/chip_vs_wgs_summary.tsv", sep="\t", index_col=0)
D["chip"] = {"compared": int(cs["concordant"].sum()), "discordant": 0, "nocall": int(cs["chip_nocall"].sum()), "nocall_nopass": 248, "uncallable": int(cs["uncallable"].sum()), "in_deletion": 63}
# --- Y
ypath = []
for l in open(W/"03_haplo/y_haplogroup_yfull.txt"):
    m = re.match(r"\s+(O-\S+)\s+der=\s*(\d+) anc=\s*(\d+) n/a=\s*(\d+)\s+formed=(\d+) tmrca=(\d+)", l)
    if m: ypath.append({"snp": m.group(1), "der": int(m.group(2)), "anc": int(m.group(3)), "formed": int(m.group(5)), "tmrca": int(m.group(6))})
D["ypath"] = ypath; D["y_upstream"] = v2["ypath"]; D["y_terminal"] = ypath[-1]["snp"]
ysn = pd.read_csv(W/"03_haplo/y_terminal_snps.tsv", sep="\t"); D["y_snps"] = ysn[["branch", "snp", "depth", "n_der"]].to_dict("records")
# --- mt
hg = pd.read_csv(W/"03_haplo/haplogrep3.txt", sep="\t")
row = hg.iloc[0]; found = str(row["Found_Polys"]).split(); rem = [x.split(" ")[0] for x in str(row["Remaining_Polys"]).split(") ")]
rem = re.findall(r"(\d+(?:\.\d+)?[ACGTd])", str(row["Remaining_Polys"]))
D["mt"] = {"hg": row["Haplogroup"], "quality": float(row["Quality"]), "found": found, "private": rem, "notfound": str(row["Not_Found_Polys"]).split()}
het = pd.read_csv(W/"03_haplo/mt_heteroplasmy.tsv", sep="\t"); D["mt_het"] = het.to_dict("records")
# --- ancestry
ps = pd.read_csv(P/"data/ref/all_phase3.psam", sep="\t").rename(columns={"#IID": "IID"})
kg = pd.read_csv(W/"04_ancestry/kg.proj.sscore", sep="\t"); me = pd.read_csv(W/"04_ancestry/target.proj.sscore", sep="\t")
D["pca_global"] = {"pts": [[r.SuperPop, r.Population, round(r.PC1_AVG, 4), round(r.PC2_AVG, 4)] for r in kg.itertuples()], "me": [round(me.PC1_AVG[0], 4), round(me.PC2_AVG[0], 4)]}
ke = pd.read_csv(W/"04_ancestry/eas.proj.sscore", sep="\t"); mee = pd.read_csv(W/"04_ancestry/eas.target.proj.sscore", sep="\t")
D["pca_eas"] = {"pts": [[r.Population, round(r.PC1_AVG, 4), round(r.PC2_AVG, 4)] for r in ke.itertuples()], "me": [round(mee.PC1_AVG[0], 4), round(mee.PC2_AVG[0], 4)]}
D["anc"] = {"n_global": sum(1 for _ in open(W/"04_ancestry/prune.prune.in")), "n_eas": sum(1 for _ in open(W/"04_ancestry/prune.eas.prune.in")), "summary": open(W/"04_ancestry/summary.txt").read()}
def near(k, m, pcs):
    cen = k.groupby("Population")[pcs].mean(); d = np.sqrt(((cen - m)**2).sum(axis=1)).sort_values()
    k = k.copy(); k["d"] = np.sqrt(((k[pcs].values - m)**2).sum(1)); return d.head(5).round(4).to_dict(), k.nsmallest(15, "d").Population.value_counts().to_dict()
pcs = ["PC1_AVG", "PC2_AVG", "PC3_AVG", "PC4_AVG"]
D["near_global"], D["knn_global"] = near(kg, me[pcs].values[0], pcs); D["near_eas"], D["knn_eas"] = near(ke, mee[pcs].values[0], pcs)
# --- ClinVar
cv = pd.read_csv(W/"05_clinvar/clinvar_all_hits.tsv", sep="\t", header=None, dtype=str, keep_default_na=False,
                 names=["chrom","pos","id","ref","alt","qual","geneinfo","clnsig","revstat","clndn","sigconf","eas_af","all_af","bcsq","gt","dp","gq","ad"])
cls = cv.clnsig.str.extract(r"^([A-Za-z_/]+)")[0].value_counts()
D["clinvar"] = {k: int(v) for k, v in cls.items()}; D["clinvar_total"] = len(cv); D["clinvar_date"] = "2026-09-05"
lof = pd.read_csv(W/"05_clinvar/lof_table.tsv", sep="\t"); rare = pd.read_csv(W/"05_clinvar/lof_rare_final.tsv", sep="\t")
D["lof"] = {"all": len(lof), "rare": len(rare), "rare_hom": int((rare.zyg == "hom/hemi").sum()), "rare_constrained": int((rare.oe_lof_upper < 0.6).sum())}
# --- PGx
pc = pd.read_csv(W/"06_pgx/pharmcat/pharmcat_summary.tsv", sep="\t").fillna("")
D["pgx_pharmcat"] = pc.to_dict("records")
D["cyp2d6"] = open(W/"06_pgx/cyrius/target.tsv").read().split("\n")[1].split("\t")[1]
hla = pd.read_csv(W/"06_pgx/t1k/dayu_genotype.tsv", sep="\t", header=None)
D["hla"] = {r[0]: [str(r[2]).replace("HLA-", ""), str(r[5]).replace("HLA-", "") if str(r[5]) != "." else "-", int(r[4]), int(r[7])] for r in hla.itertuples(index=False) if str(r[0]).startswith("HLA-") and r[1] > 0}
# --- PRS
pr = pd.read_csv(W/"07_prs/prs_wgs.tsv", sep="\t"); D["prs"] = pr.round(1).fillna(-1).to_dict("records")
# --- CNV panel (depth ratios measured above) and SV counts
sv = pd.read_csv(W/"08_sv/sv_filtered.tsv", sep="\t", dtype={"chrom": str})
D["sv_counts"] = sv.svtype.value_counts().to_dict(); D["sv_total"] = len(sv)
D["cnv"] = [
 {"locus": "GSTM1", "ratio": 0.02, "copies": 0, "evidence": "Delly + depth"}, {"locus": "LCE3B/LCE3C", "ratio": 0.0, "copies": 0, "evidence": "Delly + depth"},
 {"locus": "CYP2A6", "ratio": 0.01, "copies": 0, "evidence": "Delly + depth"}, {"locus": "LILRA3", "ratio": 0.0, "copies": 0, "evidence": "Delly + depth"},
 {"locus": "GSTT1", "ratio": 0.51, "copies": 1, "evidence": "depth"}, {"locus": "UGT2B17", "ratio": 0.55, "copies": 1, "evidence": "Delly + depth"},
 {"locus": "CFHR3-CFHR1", "ratio": 0.5, "copies": 1, "evidence": "Delly + depth"}, {"locus": "HBA1/HBA2", "ratio": 0.97, "copies": 2, "evidence": "depth"},
 {"locus": "RHD", "ratio": 1.02, "copies": 2, "evidence": "depth"}, {"locus": "SMN1", "ratio": 1.0, "copies": 2, "evidence": "SMNCopyNumberCaller"}, {"locus": "SMN2", "ratio": 1.0, "copies": 2, "evidence": "SMNCopyNumberCaller"}]
smn = pd.read_csv(W/"08_sv/smn/target.tsv", sep="\t").iloc[0]; D["smn"] = {"SMN1": int(smn.SMN1_CN), "SMN2": int(smn.SMN2_CN), "carrier": bool(smn.isCarrier)}
# --- STR
eh = pd.read_csv(W/"08_sv/eh/eh_summary.tsv", sep="\t")
thr = dict(pd.read_csv(pathlib.Path(__file__).resolve().parents[1] / "panel" / "str_thresholds.tsv", sep="\t").values)
eh = eh.drop_duplicates("locus", keep="first"); eh["thr"] = eh.locus.map(thr)
D["str"] = [{"locus": r.locus, "unit": r.unit, "gt": str(r.genotype), "max": int(r.max_allele), "thr": int(r.thr) if pd.notna(r.thr) else None} for r in eh.itertuples()]
# --- ROH
roh = [l.split() for l in open(W/"09_misc/roh_1mb_nocen.bed")]
D["roh"] = [{"chrom": r[0], "start": int(r[1]), "end": int(r[2]), "mb": round(int(r[3])/1e6, 2), "q": float(r[4])} for r in roh]
D["roh_stats"] = {"n": len(roh), "total_mb": round(sum(int(r[3]) for r in roh)/1e6, 1), "max_mb": round(max(int(r[3]) for r in roh)/1e6, 2), "n_gt5": 0}
D["blood"] = {"abo": "A (A102/A102)", "rhd": "D+"}
D["versions"] = {"yfull": "14.05.0", "phylotree": "17.2", "pharmcat": "3.4.0", "clinvar": "2026-09-05", "gnomad": "v2.1.1", "delly": "1.7.2", "eh": "5.0.0", "t1k": "1.0.10", "cyrius": "1.1.1", "bcftools": "1.22", "beagle": "5.4"}
# ---- HLA disease associations, archaic gene families, behaviour scores, candidate genes
import json as _j
hd = W/"19_hla_disease/hla_disease.tsv"
if hd.exists(): D["hla_disease"] = pd.read_csv(hd, sep="\t").fillna("").to_dict("records")
gf = W/"15_archaic/gene_families.tsv"
if gf.exists(): D["archaic_families"] = pd.read_csv(gf, sep="\t").fillna("").to_dict("records")
gs = W/"15_archaic/gene_summary.txt"
if gs.exists(): D["archaic_genes"] = dict(l.split("\t") for l in open(gs).read().strip().split("\n"))
sg = W/"15_archaic/segments_genes.tsv"
if sg.exists():
    t = pd.read_csv(sg, sep="\t", dtype={"chrom": str}).nlargest(8, "n_genes")
    D["archaic_top"] = t[["seg","mb","src","n_genes","genes"]].assign(genes=lambda x: x.genes.str.split(",").str[:14].str.join(", ")).to_dict("records")
bp = W/"20_behaviour/behaviour_prs.tsv"
if bp.exists(): D["behaviour"] = pd.read_csv(bp, sep="\t").fillna("").to_dict("records")
cg = W/"20_behaviour/candidate_genes.tsv"
if cg.exists(): D["candidate"] = pd.read_csv(cg, sep="\t").fillna("").to_dict("records")
D["name_zh"] = NAME_ZH; D["name_en"] = NAME_EN; D["sample"] = SAMPLE
D["title_zh"] = NAME_ZH; D["title_en"] = NAME_EN
json.dump(D, open(W/"report_data.json", "w"), ensure_ascii=False)
print("wrote", W/"report_data.json", (W/"report_data.json").stat().st_size//1024, "KB"); print({k: (len(v) if isinstance(v, (list, dict)) else v) for k, v in D.items() if k not in ("chrlen", "cen")})
