#!/usr/bin/env python3
"""Y-chromosome haplogroup calling by greedy descent of the ISOGG tree.

Uses the ISOGG SNP table and YCC tree written by yhaplo (results/yhaplo/), and the
plus-strand Y genotypes from data/genome.parquet. At every node it scores each child
by (#derived, #ancestral) and descends into the child with the strongest derived
support, stopping when no child has derived-majority evidence.
"""
import re, pathlib, pandas as pd

from config import *
ROOT = pathlib.Path(__file__).resolve().parents[1]
YH = RESULTS / "yhaplo"
# Step 0: run yhaplo once to obtain the ISOGG SNP table and YCC tree files (its own haplogroup call is not used)
if not (YH / "isogg.snps.cleaned.2016.01.04.txt").exists():
    import subprocess
    _df = pd.read_parquet(DATA / "genome.parquet"); _y = _df[(_df.chrom == "Y") & (~_df.is_missing)]
    with open(DATA / "sample.Y.genos.txt", "w") as f:
        f.write("ID\t" + "\t".join(map(str, _y.pos)) + "\n"); f.write(f"{SAMPLE}\t" + "\t".join(_y.a1) + "\n")
    YH.mkdir(parents=True, exist_ok=True)
    subprocess.run(["yhaplo", "-i", str(DATA / "sample.Y.genos.txt"), "-o", str(YH), "--all_aux_output"], capture_output=True)
df = pd.read_parquet(DATA / "genome.parquet")
y = df[(df.chrom == "Y") & (~df.is_missing)].set_index("pos")

iso = pd.read_csv(YH / "isogg.snps.cleaned.2016.01.04.txt", sep=r"\s+", header=None,
                  names=["name", "hg", "pos", "mut"])
iso["anc"], iso["der"] = iso.mut.str[0], iso.mut.str[-1]
iso = iso[iso.pos.isin(y.index)].copy()
iso["geno"] = iso.pos.map(y.a1)
iso["state"] = ["der" if g == d else ("anc" if g == a else "other")
                for g, a, d in zip(iso.geno, iso.anc, iso.der)]

# --- parse newick tree of haplogroup names -> children dict
nwk = (YH / "y.tree.ycc.2016.01.04.nwk").read_text().strip().rstrip(";")
children, stack, name, i = {}, [], "", 0
tokens = re.findall(r"\(|\)|,|[^(),;]+", nwk)
# recursive descent
def parse(tokens):
    pos = 0
    def node():
        nonlocal pos
        kids = []
        if tokens[pos] == "(":
            pos += 1
            while True:
                kids.append(node())
                if tokens[pos] == ",": pos += 1; continue
                if tokens[pos] == ")": pos += 1; break
        label = tokens[pos].strip() if pos < len(tokens) and tokens[pos] not in "()," else ""
        if label: pos += 1
        children[label] = [k for k in kids]
        return label
    return node()
root = parse(tokens)

counts = iso.groupby(["hg", "state"]).size().unstack(fill_value=0)
for c in ["der", "anc", "other"]:
    if c not in counts: counts[c] = 0

def score(hg):
    d = a = 0
    for alias in hg.split("/"):
        if alias in counts.index:
            d += int(counts.loc[alias, "der"]); a += int(counts.loc[alias, "anc"])
    return d, a

path, node = [], root
while True:
    d, a = score(node)
    path.append((node, d, a))
    best, best_key = None, None
    for ch in children.get(node, []):
        d, a = score(ch)
        if d == 0: continue
        if d > a:  # derived majority
            key = (d - a, d)
            if best is None or key > best_key: best, best_key = ch, key
    if best is None: break
    node = best

terminal = node
# defining SNPs of terminal node observed derived
term_snps = iso[(iso.hg.isin(terminal.split("/"))) & (iso.state == "der")].name.tolist()
out = RESULTS / "y_haplogroup.txt"
with open(out, "w") as f:
    f.write(f"Terminal haplogroup (ISOGG 2016.01.04): {terminal}\n")
    f.write(f"Defining SNPs observed derived: {', '.join(term_snps[:20])}\n\n")
    f.write("Path (node, #derived, #ancestral):\n")
    for n, d, a in path: f.write(f"  {n:12s} der={d:4d} anc={a:4d}\n")
    f.write("\nChildren of terminal node and their evidence (unresolved below here):\n")
    for ch in children.get(terminal, []):
        d, a = score(ch); f.write(f"  {ch:12s} der={d:4d} anc={a:4d}\n")
print(out.read_text())
iso.sort_values(["hg","pos"]).to_csv(RESULTS / "y_isogg_marker_states.tsv", sep="\t", index=False)
