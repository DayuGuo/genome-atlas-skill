#!/usr/bin/env python
"""Screen the mitochondrial genome for the disease variants in panel/mt_disease.tsv, directly from the pileup.
Depth here is normally several thousand, so a clean zero is a strong negative; anything above about 1% of reads
is worth a second look, below that is sequencing noise."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403

import pandas as pd

W = f"{P}/wgs/03_haplo"
PANEL = pathlib.Path(__file__).resolve().parents[1] / "panel" / "mt_disease.tsv"
counts = {}
for l in open(f"{W}/mt_pileup_ad.tsv"):
    pos, ref, alt, ad = l.rstrip("\n").split("\t")
    alts = alt.split(","); ads = list(map(int, ad.split(",")))
    d = {ref: ads[0]}
    for a, n in zip(alts, ads[1:]):
        if a != "<*>": d[a] = n
    counts[int(pos)] = d
rows = []
for r in pd.read_csv(PANEL, sep="\t").itertuples():
    c = counts.get(int(r.pos), {}); dp = sum(c.values()); n = c.get(r.alt, 0)
    rows.append(dict(pos=int(r.pos), change=f"{r.ref}>{r.alt}", condition=r.condition,
                     depth=dp, alt_reads=n, fraction=round(n/dp, 5) if dp else None,
                     call="present" if dp and n/dp >= 0.01 else ("absent" if dp else "no coverage")))
df = pd.DataFrame(rows); df.to_csv(f"{W}/mt_disease.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(df.to_string(index=False))
hit = df[df.call == "present"]
print(f"\n{len(hit)} of {len(df)} screened positions carry the variant at >=1% of reads")
if len(hit): print(hit.to_string(index=False))
