#!/usr/bin/env python3
"""Assemble results_v2/report_data_v2.json: every number in report v2 comes from result files, nothing hand-typed."""
import json, pathlib, numpy as np, pandas as pd
from config import *
ROOT = pathlib.Path(__file__).resolve().parents[1]; R = RESULTS; V = RESULTS2; AU = AUDIT
# ---- base bundle from result files (chromosomes, PCA coordinates, Y chain, mt sites, ROH, ClinVar, exploratory)
S = {}
summ = pd.read_csv(R/"summary_basic.tsv", sep="\t", dtype={"chrom":str}); S["chrom"] = summ.to_dict("records")
kg = pd.read_csv(R/"ancestry"/"kg.common.proj.sscore", sep="\t"); me = pd.read_csv(R/"ancestry"/"sample.proj.sscore", sep="\t")
S["pca_global"] = {"pts": [[r.SuperPop, r.Population, round(r.PC1_AVG,4), round(r.PC2_AVG,4)] for r in kg.itertuples()], "me": [round(me.PC1_AVG[0],4), round(me.PC2_AVG[0],4)]}
kge = pd.read_csv(R/"ancestry"/"eas.kg.proj.sscore", sep="\t"); mee = pd.read_csv(R/"ancestry"/"eas.sample.proj.sscore", sep="\t")
S["pca_eas"] = {"pts": [[r.Population, round(r.PC1_AVG,4), round(r.PC2_AVG,4)] for r in kge.itertuples()], "me": [round(mee.PC1_AVG[0],4), round(mee.PC2_AVG[0],4)]}
h = pd.read_csv(R/"hgdp"/"hgdp.common.proj.sscore", sep="\t"); hm = pd.read_csv(R/"hgdp"/"sample.proj.sscore", sep="\t")
S["pca_hgdp"] = {"pts": [[r.population, round(r.PC1_AVG,4), round(r.PC2_AVG,4), round(r.PC3_AVG,4), round(r.PC4_AVG,4)] for r in h.itertuples()], "me": [round(hm.PC1_AVG[0],4), round(hm.PC2_AVG[0],4), round(hm.PC3_AVG[0],4), round(hm.PC4_AVG[0],4)]}
_ev = np.loadtxt(R/"hgdp"/"pca.eigenval"); S["pve_hgdp"] = [round(float(x),1) for x in (_ev/_ev.sum()*100)[:4]]
S["prs"] = pd.read_csv(R/"prs_summary.tsv", sep="\t").to_dict("records"); S["prs_bio"] = pd.read_csv(R/"prs_biomarkers.tsv", sep="\t").to_dict("records")
_tp = pd.read_csv(R/"trait_panel_report.tsv", sep="\t").fillna(""); S["traits"] = _tp[["rsid","gene","category","trait","genotype","n_effect","interpretation"]].to_dict("records")
cv = pd.read_csv(R/"clinvar_screen.tsv", sep="\t"); S["clinvar"] = cv.significance.fillna("unclassified").astype(str).str.split(r"[|/]").str[0].replace("", "unclassified").value_counts().to_dict(); S["clinvar_total"] = int(len(cv))
S["roh"] = pd.read_csv(R/"roh.tsv", sep="\t", dtype={"chrom":str}).to_dict("records")
S["chrlen"] = {"1":249250621,"2":243199373,"3":198022430,"4":191154276,"5":180915260,"6":171115067,"7":159138663,"8":146364022,"9":141213431,"10":135534747,"11":135006516,"12":133851895,"13":115169878,"14":107349540,"15":102531392,"16":90354753,"17":81195210,"18":78077248,"19":59128983,"20":63025520,"21":48129895,"22":51304566}
S["cen"] = {"1":125.0,"2":93.3,"3":91.0,"4":50.4,"5":48.4,"6":61.0,"7":59.9,"8":45.6,"9":49.0,"10":40.2,"11":53.7,"12":35.8,"13":17.9,"14":17.6,"15":19.0,"16":36.6,"17":24.0,"18":17.2,"19":26.5,"20":27.5,"21":13.2,"22":14.7}
# Y descent chain from the called path: label each node with its first derived defining SNP, count derived sites
_ym = pd.read_csv(R/"y_isogg_marker_states.tsv", sep="\t"); _ytxt = open(R/"y_haplogroup.txt").read()
_nodes = []
for l in _ytxt.split("Path (node")[1].split("\n")[1:]:
    if l.startswith("Children"): break
    if l.strip(): _nodes.append(l.split()[0])
_nodes = [n for n in _nodes if n != "Root"]
_start = max(0, len(_nodes) - 9)  # last 9 steps of the path for the chain figure
S["ypath"] = []
for hg in _nodes[_start:]:
    d = _ym[_ym.hg.isin(hg.split("/")) & (_ym.state == "der")]
    _rank = lambda x: (0 if x.startswith("M") and x[1:2].isdigit() else 1 if x.startswith("P") and x[1:2].isdigit() else 2 if x.startswith("Page") else 3 if x.startswith("F") and x[1:2].isdigit() else 4, len(x), x)
    pref = sorted(d.name.tolist(), key=_rank)
    S["ypath"].append({"snp": (hg.split("/")[0][0] + "-" + pref[0]) if len(pref) else hg, "hg": hg, "n": int(len(d))})
_mt = pd.read_csv(R/"mt_sites_vs_rCRS.tsv", sep="\t"); _mt = _mt[_mt["diff"] & _mt.a1.isin(list("ACGT"))]
S["mt_variants"] = [f"{p}{a}" for p, a in zip(_mt.pos, _mt.a1)]
_hg = pd.read_csv(R/"mt_haplogrep3"/"haplogroups.txt", sep="\t"); S["mt_private"] = [x.split()[0] for x in str(_hg.Remaining_Polys.iloc[0]).split(")") if x.strip() and "(" in x]
S["ns"] = pd.read_csv(R/"fun_chrom_north_south.tsv", sep="\t")[["chrom","n_snps","ns_index"]].to_dict("records")
_xp = pd.read_csv(R/"fun"/"x_pcs.tsv", sep="\t"); _xd = pd.read_csv(R/"fun"/"x_sample.tsv", sep="\t")
S["xpca"] = {"pts": [[r.SuperPop, r.Population, round(r.PC1_AVG,4), round(r.PC2_AVG,4)] for r in _xp.itertuples()], "me": [round(_xd.PC1_AVG[0],4), round(_xd.PC2_AVG[0],4)]}
_rare = pd.read_csv(R/"fun_rarest_variants.tsv", sep="\t").fillna(""); S["rare"] = _rare.head(15).to_dict("records")
_ra = pd.read_csv(R/"fun"/"rare_all.tsv", sep="\t"); S["rare_n"] = int(len(_ra)); S["rare_absent"] = int((_ra.freq_eas==0).sum())
S.pop("refsim", None); S.pop("refsim_pts", None)     # module removed by audit
qc = json.load(open(AU/"qc_and_harmonization.json")); S["qc"] = qc
summ = pd.read_csv(R/"summary_basic.tsv", sep="\t", dtype={"chrom":str})
S["kpi"] = {"sites": int(qc["n_rows"]), "call_rate": qc["call_rate"], "het_auto": qc["het_autosomal"], "x_called": qc["X_called"], "x_het": qc["X_het"],
            "x_par_sites": qc["X_PAR_sites"], "x_par_het": qc["X_PAR_het"], "y_called": qc["Y_called"], "y_share": qc["Y_called_share"],
            "n_rsid": qc["n_rsid"], "n_wid": qc["n_wegene_id"], "indels": qc["indel_sites"]}
roh = pd.read_csv(R/"roh.tsv", sep="\t"); S["roh_stats"] = {"n": len(roh), "total_mb": round(roh.length_mb.sum(),1), "max_mb": float(roh.length_mb.max()), "n_gt5": int((roh.length_mb>5).sum())}
# 1000G nearest populations, eigenvalue-weighted Euclidean on PC1-6 + kNN (replaces standardised distance)
def nearest(sscore, dayu, eigenval, popcol, k=6, knn=15):
    ref = pd.read_csv(sscore, sep="\t"); me = pd.read_csv(dayu, sep="\t"); ev = np.loadtxt(eigenval)
    pcs = [c for c in ref.columns if c.startswith("PC")][:k]; w = np.sqrt(ev[:k])
    X = ref[pcs].to_numpy()*w; b = me[pcs].iloc[0].to_numpy()*w
    cent = pd.DataFrame(X, index=ref[popcol]).groupby(level=0).mean()
    d = pd.Series(np.linalg.norm(cent.to_numpy()-b, axis=1), index=cent.index).sort_values()
    ind = np.linalg.norm(X-b, axis=1); order = np.argsort(ind)[:knn]; kn = ref[popcol].iloc[order].value_counts()
    return [[p, round(float(x),4)] for p,x in d.items()], {p:int(n) for p,n in kn.items()}
S["near_global"], S["knn_global"] = nearest(R/"ancestry"/"kg.common.proj.sscore", R/"ancestry"/"sample.proj.sscore", R/"ancestry"/"kg.pca.eigenval", "Population")
S["near_eas"], S["knn_eas"] = nearest(R/"ancestry"/"eas.kg.proj.sscore", R/"ancestry"/"eas.sample.proj.sscore", R/"ancestry"/"eas.pca.eigenval", "Population")
S["near_hgdp"], S["knn_hgdp"] = nearest(R/"hgdp"/"hgdp.common.proj.sscore", R/"hgdp"/"sample.proj.sscore", R/"hgdp"/"pca.eigenval", "population")
S["n_pruned_global"] = sum(1 for _ in open(R/"ancestry"/"prune.prune.in")); S["n_pruned_eas"] = sum(1 for _ in open(R/"ancestry"/"eas.prune.prune.in")); S["n_pruned_hgdp"] = sum(1 for _ in open(R/"hgdp"/"prune.prune.in"))
S["hgdp_n"] = {p: int(n) for p, n in pd.read_csv(R/"hgdp"/"hgdp.common.proj.sscore", sep="\t").population.value_counts().items()}
# N/S index: the sample position within CHB distribution in SD units (audit: 1.15 is a coordinate, not a fraction)
ns = pd.read_csv(R/"fun_chrom_north_south.tsv", sep="\t"); S["ns_mean"] = round(float(ns.ns_index.mean()),2)
eas = pd.read_csv(R/"ancestry"/"eas.kg.proj.sscore", sep="\t"); me = pd.read_csv(R/"ancestry"/"eas.sample.proj.sscore", sep="\t")
chb = eas[eas.Population==AXIS[1]].PC1_AVG; chs = eas[eas.Population==AXIS[0]].PC1_AVG; x = me.PC1_AVG.iloc[0]
S["ns_genome"] = {"index": round(float((x-chs.mean())/(chb.mean()-chs.mean())),2), "sd_from_chb": round(float((x-chb.mean())/chb.std()),2), "sd_from_chs": round(float((x-chs.mean())/chs.std()),2)}
# Y path support counts + conflicts, from the path actually called
ym = pd.read_csv(R/"y_isogg_marker_states.tsv", sep="\t")
ytxt = open(R/"y_haplogroup.txt").read()
ypath_nodes = []
for l in ytxt.split("Path (node")[1].split("\n")[1:]:
    if l.startswith("Children"): break
    if l.strip(): ypath_nodes.append(l.split()[0])
ypath_nodes = [n for n in ypath_nodes if n != "Root"]
S["y_path_nodes"] = ypath_nodes
S["y_conflicts"] = {hg: {"der": int((ym.hg.isin(hg.split("/")) & (ym.state=="der")).sum()), "anc": int((ym.hg.isin(hg.split("/")) & (ym.state=="anc")).sum())} for hg in ypath_nodes}
S["y_terminal"] = open(R/"y_haplogroup.txt").readline().split(":")[-1].strip()
# mt
hg = pd.read_csv(R/"mt_haplogrep3"/"haplogroups.txt", sep="\t"); S["mt"] = {"hg": hg.Haplogroup.iloc[0], "quality": float(hg.Quality.iloc[0]), "sites": int(qc["MT_called"]), "coverage_pct": round(100*qc["MT_called"]/16569,1), "n_found": len(str(hg.Found_Polys.iloc[0]).split()), "n_private": len(S["mt_private"])}
# traits joined with audit verdicts
tv = pd.read_csv(AU/"05_trait_validation.tsv", sep="\t").fillna("")
S["traits_v2"] = tv[tv.verdict!="REMOVE"][["rsid","gene","category","genotype","direct_or_proxy","evidence_level","phenotype","verdict","corrected_reading_or_reason","source"]].to_dict("records")
S["traits_removed"] = tv[tv.verdict=="REMOVE"].rsid.tolist()  # verdict column may be filled by the analyst
S["pgx_v2"] = pd.read_csv(AU/"03_pgx_validation.tsv", sep="\t").fillna("").to_dict("records")
tp = pd.read_csv(R/"trait_panel_report.tsv", sep="\t").fillna(""); S["traits_zh_orig"] = dict(zip(tp.rsid, tp.interpretation))
S["pop"] = {"superpop": REF_SUPERPOP, "subpops": SUBPOPS, "axis": AXIS, "hgdp_focus": HGDP_FOCUS, "hgdp_second": HGDP_SECOND}
# PRS: validation + imputed
pv = pd.read_csv(AU/"04_prs_validation.tsv", sep="\t"); S["prs_validation"] = pv.to_dict("records")
S["prs_validation_summary"] = {"n": len(pv), "rho_min": float(pv.spearman_full_vs_subset.min()), "rho_max": float(pv.spearman_full_vs_subset.max()), "dev_median_median": float(pv.pct_dev_median.median()), "dev_p90_median": float(pv.pct_dev_p90.median()), "dev_p90_max": float(pv.pct_dev_p90.max()), "cov_min": float(pv.coverage_pct.min()), "cov_max": float(pv.coverage_pct.max())}
pi = pd.read_csv(V/"prs_imputed.tsv", sep="\t"); S["prs_imputed"] = pi.to_dict("records")
S["prs_imputed_summary"] = {"cov_min": float(pi.coverage_of_score_pct.min()), "cov_max": float(pi.coverage_of_score_pct.max()), "w_min": float(pi.weight_retained_pct.min()), "w_max": float(pi.weight_retained_pct.max())}
ic = pd.read_csv(V/"impute_check"/"summary.tsv", sep="\t"); S["impute_check"] = ic.to_dict("records")
msr = pd.read_csv(V/"impute_check"/"masked_site_results.tsv", sep="\t")
S["impute_check_overall"] = {"n_masked": len(msr), "conc_all": round(float(msr.correct.mean()),4), "het_conc_all": round(float(msr[msr.het].correct.mean()),4), "n_gp99": int((msr.gp>=0.99).sum()), "conc_gp99": round(float(msr[msr.gp>=0.99].correct.mean()),4)}
S["imputed_counts"] = json.load(open(V/"imputed_counts.json")); S["imputed_variants_hq"] = S["imputed_counts"]["common"]
S["clinvar_date"] = qc["clinvar_fileDate"]
S["lifted_sites"] = sum(1 for _ in open(R/"hgdp"/"sample.hg38.range"))
S["impute_input_sites"] = sum(1 for ch in range(1,23) for l in open(V/"imputed"/f"target.chr{ch}.vcf") if not l.startswith("#"))
json.dump(S, open(V/"report_data_v2.json","w"), ensure_ascii=False)
print("wrote report_data_v2.json", (V/"report_data_v2.json").stat().st_size//1024, "KB")
print("near_global", S["near_global"][:4], S["knn_global"]); print("near_eas", S["near_eas"][:3], S["knn_eas"]); print("ns", S["ns_genome"]); print("y_conflicts", S["y_conflicts"])
