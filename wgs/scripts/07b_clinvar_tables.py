#!/usr/bin/env python
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import pandas as pd, numpy as np, re, gzip
W = f"{P}/wgs/05_clinvar"
cols = ["chrom","pos","id","ref","alt","qual","geneinfo","clnsig","revstat","clndn","sigconf","eas_af","all_af","bcsq","gt","dp","gq","ad"]
cv = pd.read_csv(f"{W}/clinvar_all_hits.tsv", sep="\t", header=None, names=cols, dtype=str, na_values=["."], keep_default_na=False)
cv["gene"] = cv.geneinfo.str.split(":").str[0]
stars = {"practice_guideline":4, "reviewed_by_expert_panel":3, "criteria_provided,_multiple_submitters,_no_conflicts":2,
         "criteria_provided,_single_submitter":1, "criteria_provided,_conflicting_classifications":1, "criteria_provided,_conflicting_interpretations":1}
cv["stars"] = cv.revstat.map(lambda s: stars.get(s, 0))
def zyg(g):
    g = g.replace("|","/")
    if g in ("1","1/1"): return "hom/hemi"
    if g in ("0/1","1/0"): return "het"
    return g
cv["zyg"] = cv["gt"].map(zyg)
print("ClinVar-annotated PASS variants:", len(cv)); print(cv.clnsig.str.extract(r"^([A-Za-z_]+)")[0].value_counts().head(12).to_string())
plp = cv[cv.clnsig.str.contains(r"^(Pathogenic|Likely_pathogenic|Pathogenic/Likely_pathogenic)", regex=True) & ~cv.clnsig.str.contains("Conflicting")]
plp = plp.sort_values(["stars","gene"], ascending=[False,True])
ACMG = set("""ACTA2 ACTC1 ACVRL1 APC APOB ATP7B BAG3 BMPR1A BRCA1 BRCA2 BTD CACNA1S CALM1 CALM2 CALM3 CASQ2 COL3A1 DES DSC2 DSG2 DSP ENG FBN1 FLNC GAA GLA HFE HNF1A KCNH2 KCNQ1 LDLR LMNA MAX MEN1 MLH1 MSH2 MSH6 MUTYH MYBPC3 MYH11 MYH7 MYL2 MYL3 NF2 OTC PALB2 PCSK9 PKP2 PMS2 PRKAG2 PTEN RB1 RBM20 RET RPE65 RYR1 RYR2 SCN5A SDHAF2 SDHB SDHC SDHD SMAD3 SMAD4 STK11 TGFBR1 TGFBR2 TMEM127 TMEM43 TNNC1 TNNI3 TNNT2 TP53 TPM1 TRDN TSC1 TSC2 TTN TTR VHL WT1""".split())
plp["acmg_sf"] = plp.gene.isin(ACMG)
out = plp[["chrom","pos","id","ref","alt","gene","clnsig","stars","revstat","clndn","eas_af","zyg","dp","gq","ad","acmg_sf","bcsq"]]
out.to_csv(f"{W}/clinvar_PLP.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(f"\nP/LP (non-conflicting): {len(out)}  | >=1 star: {(out.stars>=1).sum()}  | >=2 stars: {(out.stars>=2).sum()} | hom/hemi: {(out.zyg=='hom/hemi').sum()} | ACMG SF genes: {out.acmg_sf.sum()}")
print(out[out.stars>=1][["chrom","pos","id","gene","clnsig","stars","clndn","eas_af","zyg","dp","acmg_sf"]].to_string(index=False))
# conflicting with P entries in genes of interest
conf = cv[cv.clnsig.str.contains("Conflicting") & cv.sigconf.str.contains("athogenic", na=False)].copy()
def counts(sc):
    d = dict(re.findall(r"([A-Za-z_]+)_\((\d+)\)", sc or "")); g = lambda k: int(d.get(k, 0))
    return g("Pathogenic")+g("Likely_pathogenic"), g("Benign")+g("Likely_benign"), g("Uncertain_significance")
conf[["n_P","n_B","n_VUS"]] = pd.DataFrame([counts(x) for x in conf.sigconf], index=conf.index)
conf["carrier_candidate"] = (conf.n_P >= 3) & (conf.n_P >= 3*conf.n_B)
conf[["chrom","pos","id","ref","alt","gene","sigconf","n_P","n_B","n_VUS","carrier_candidate","clndn","eas_af","all_af","zyg","dp","ad"]].to_csv(f"{W}/clinvar_conflicting_with_P.tsv", sep="\t", index=False)
print("\nConflicting with a Pathogenic submission:", len(conf), " carrier candidates:", conf.carrier_candidate.sum())
print(conf[conf.carrier_candidate][["chrom","pos","id","gene","sigconf","eas_af","zyg"]].to_string(index=False))
# risk factor / drug response / protective
other = cv[cv.clnsig.str.contains("risk_factor|drug_response|protective|association", regex=True)]
other[["chrom","pos","id","gene","clnsig","clndn","eas_af","zyg","stars"]].to_csv(f"{W}/clinvar_other.tsv", sep="\t", index=False)
print("risk_factor/drug_response/protective/association entries:", len(other))
# --- LoF
lc = ["chrom","pos","id","ref","alt","qual","eas_af","all_af","clnsig","bcsq","gt","dp","gq","ad"]
lof = pd.read_csv(f"{W}/lof_all.tsv", sep="\t", header=None, names=lc, dtype=str, na_values=["."], keep_default_na=False)
def parse(b):
    for c in b.split(","):
        f = c.split("|")
        if len(f) >= 4 and re.search("stop_gained|frameshift|splice_acceptor|splice_donor|start_lost", f[0]) and f[3] == "protein_coding": return f[0], f[1], f[2]
    return None, None, None
lof[["csq","gene","tx"]] = pd.DataFrame([parse(b) for b in lof.bcsq], index=lof.index)
lof = lof[lof.gene.notna()]
lof["eas_af"] = pd.to_numeric(lof.eas_af, errors="coerce"); lof["all_af"] = pd.to_numeric(lof.all_af, errors="coerce"); lof["zyg"] = lof["gt"].map(zyg)
con = pd.read_csv(GNOMAD_CONSTRAINT, sep="\t", compression="gzip", usecols=["gene","transcript","pLI","oe_lof_upper","oe_lof"])
con = con.sort_values("oe_lof_upper").drop_duplicates("gene")
lof = lof.merge(con, on="gene", how="left")
lof["rare"] = (lof.eas_af.fillna(0) < 0.01) & (lof.all_af.fillna(0) < 0.01)
lof["in_1000G"] = lof.eas_af.notna()
lof = lof.sort_values(["oe_lof_upper"], na_position="last")
lof[["chrom","pos","id","ref","alt","gene","csq","zyg","eas_af","all_af","in_1000G","pLI","oe_lof_upper","dp","gq","ad","clnsig"]].to_csv(f"{W}/lof_table.tsv", sep="\t", index=False)
print(f"\nLoF (protein_coding): {len(lof)}; rare (EAS AF<1% or absent): {lof.rare.sum()}; hom rare: {((lof.zyg=='hom/hemi')&lof.rare).sum()}; rare in constrained genes (LOEUF<0.35): {(lof.rare&(lof.oe_lof_upper<0.35)).sum()}")
print(lof[lof.rare & (lof.oe_lof_upper < 0.6)][["chrom","pos","id","ref","alt","gene","csq","zyg","eas_af","all_af","pLI","oe_lof_upper","dp","ad"]].to_string(index=False))
print("\nhomozygous/hemizygous LoF (any freq), genes:", ", ".join(sorted(set(lof[lof.zyg=="hom/hemi"].gene))))
