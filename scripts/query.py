#!/usr/bin/env python3
"""Interactive lookup tool for the genome database.

  python3 scripts/query.py rs671 rs1229984          # by rsid (+ Ensembl GRCh37 annotation)
  python3 scripts/query.py region 12:112200000-112300000
  python3 scripts/query.py gene ALDH2                # all array sites inside a gene (GRCh37 coords via Ensembl)
  python3 scripts/query.py sql "select chrom,count(*) from snp where is_het group by chrom"
Add --raw to skip the Ensembl call.
"""
import sys, json, pathlib, requests, duckdb, pandas as pd
from config import *
ROOT = pathlib.Path(__file__).resolve().parents[1]
con = duckdb.connect(str(DATA/"genome.duckdb"), read_only=True)
pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 80); pd.set_option("display.max_rows", 500)
ENS = "https://grch37.rest.ensembl.org"
H = {"Content-Type": "application/json", "Accept": "application/json"}

def annotate(rsids):
    try:
        r = requests.post(f"{ENS}/variation/homo_sapiens", headers=H, json={"ids": rsids, "phenotypes": 1}, timeout=60)
        r.raise_for_status(); return r.json()
    except Exception as e:
        print(f"[ensembl unavailable: {e}]"); return {}

def show_rsids(rsids, raw=False):
    g = con.execute("select rsid, chrom, pos, genotype from snp where rsid in (select unnest(?))", [rsids]).df()
    missing = set(rsids) - set(g.rsid)
    if not raw and len(g):
        ann = annotate(g.rsid.tolist())
        rows = []
        for r in g.itertuples():
            e = ann.get(r.rsid, {})
            mp = [m for m in e.get("mappings", []) if m.get("assembly_name") == "GRCh37"]
            phen = sorted({p.get("trait", "") for p in e.get("phenotypes", []) if p.get("trait")})
            rows.append(dict(alleles=mp[0]["allele_string"] if mp else "", ancestral=e.get("ancestral_allele"),
                             MAF=e.get("MAF"), clinical=";".join(e.get("clinical_significance", [])),
                             consequence=e.get("most_severe_consequence"), phenotypes="; ".join(phen)[:120]))
        g = pd.concat([g, pd.DataFrame(rows)], axis=1)
    print(g.to_string(index=False))
    if missing: print("not on array:", ", ".join(sorted(missing)))

def show_region(chrom, start, end):
    g = con.execute("select rsid, chrom, pos, genotype, is_het from snp where chrom=? and pos between ? and ? order by pos",
                    [chrom, start, end]).df()
    print(g.to_string(index=False)); print(f"{len(g)} sites, {int(g.is_het.sum())} heterozygous")

def show_gene(name):
    r = requests.get(f"{ENS}/lookup/symbol/homo_sapiens/{name}?expand=0", headers=H, timeout=60)
    if r.status_code != 200: print("gene not found in Ensembl GRCh37:", name); return
    j = r.json(); print(f"{name}: chr{j['seq_region_name']}:{j['start']}-{j['end']} ({j.get('biotype')}, {j.get('description','')})")
    show_region(j["seq_region_name"], j["start"], j["end"])

if __name__ == "__main__":
    a = sys.argv[1:]; raw = "--raw" in a; a = [x for x in a if x != "--raw"]
    if not a: print(__doc__); sys.exit()
    if a[0] == "region":
        c, rng = a[1].split(":"); s, e = map(int, rng.replace(",", "").split("-")); show_region(c, s, e)
    elif a[0] == "gene": show_gene(a[1])
    elif a[0] == "sql": print(con.execute(a[1]).df().to_string(index=False))
    else: show_rsids(a, raw)
