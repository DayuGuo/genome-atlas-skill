#!/usr/bin/env python3
"""Build the evidence-graded trait table and the PGx genotype table (genotype-agnostic panels + this sample's calls).
Outputs: WORK/audit/05_trait_validation.tsv, WORK/audit/03_pgx_validation.tsv.
Interpretation sentences are NOT generated here: the analyst/agent writes them into report_text.yaml after
reading these tables (see docs/AUDIT_CHECKLIST.md)."""
import pathlib, pandas as pd, duckdb
from config import *
tp = pd.read_csv(RESULTS / "trait_panel_report.tsv", sep="\t").fillna("")
ev = pd.read_csv(PANEL / "trait_evidence.tsv", sep="\t").fillna("")
t = tp.merge(ev, on="rsid", how="left").fillna({"direct_or_proxy": "DIRECT_GENOTYPE", "evidence_level": "D", "ancestry_note": "", "source": ""})
t["direct_or_proxy"] = t.apply(lambda r: "NOT_COVERED" if r.genotype == "" else r.direct_or_proxy, axis=1)
t["evidence_level"] = t.apply(lambda r: "—" if r.genotype == "" else r.evidence_level, axis=1)
t["verdict"] = ""; t["corrected_reading_or_reason"] = ""
t = t.rename(columns={"trait": "phenotype", "interpretation": "current_interpretation"})
t[["rsid","gene","category","genotype","strand_check","ensembl_alleles","direct_or_proxy","evidence_level","ancestry_note","phenotype","current_interpretation","verdict","corrected_reading_or_reason","source"]].to_csv(AUDIT / "05_trait_validation.tsv", sep="\t", index=False)
print("trait rows:", len(t), " levels:", t.evidence_level.value_counts().to_dict())

con = duckdb.connect(str(DATA / "genome.duckdb"), read_only=True)
pg = pd.read_csv(PANEL / "pgx_panel.tsv", sep="\t").fillna("")
rows = []
for r in pg.itertuples():
    rs = r.rsids.split(";"); maps = r.allele_map.split(";")
    g = con.execute("select rsid, genotype from snp where rsid in (select unnest(?))", [rs]).df().set_index("rsid").genotype.to_dict()
    genos = [g.get(x, "—") for x in rs]
    carried = []
    for x, gt, mp in zip(rs, genos, maps + [""] * (len(rs) - len(maps))):
        eff = mp.split("=")[0] if "=" in mp else ""
        if gt not in ("—", "--") and eff and eff in gt: carried.append(f"{x}:{mp}×{gt.count(eff)}")
    rows.append(dict(gene=r.gene, variants=" / ".join(rs), genotype=" / ".join(genos), direct_or_tag=("TAG_SNP" if "tag" in r.allele_map.lower() else "DIRECT"),
                     star_allele_coverage=r.star_allele_coverage, effect_alleles_carried="; ".join(carried) or "none of the tested effect alleles",
                     not_on_array="; ".join(x for x, gt in zip(rs, genos) if gt == "—"), phenotype_assignable="", corrected_interpretation="", source=r.source))
pd.DataFrame(rows).to_csv(AUDIT / "03_pgx_validation.tsv", sep="\t", index=False)
print("pgx genes:", len(rows)); print(pd.DataFrame(rows)[["gene","genotype","effect_alleles_carried"]].to_string(index=False))
