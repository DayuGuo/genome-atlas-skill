#!/usr/bin/env python
"""Y haplogroup by greedy descent of the YFull tree (v14) using pileup allele counts from the CRAM at every branch-defining
SNP with a known hg19 position (ybrowse). A branch is 'derived' if the majority of its typed SNPs (depth>=2) show the
derived allele. Reports the path, per-branch support, and the terminal branch; also lists private (novel) Y variants."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from wgsconfig import *  # noqa: F401,F403 -- P, W, REF, TOOLS, SAMPLE, THREADS ...

import json, csv, subprocess, collections, os, sys
import pandas as pd
P = str(P); W = f"{P}/wgs/03_haplo"; CRAM = f"{P}/wgs/00_input/{SAMPLE}.cram"; REF = FASTA
# 1. SNP index: name -> (pos, anc, der)
idx = {}
with open(f"{YTREE}/snps_hg19.csv", newline="") as f:
    for r in csv.DictReader(f):
        try: pos = int(r["start"])
        except: continue
        a, d = r["allele_anc"].upper(), r["allele_der"].upper()
        if len(a) == 1 and len(d) == 1 and a in "ACGT" and d in "ACGT":
            idx.setdefault(r["Name"], (pos, a, d))
print("ybrowse SNPs with hg19 pos:", len(idx), file=sys.stderr)
tree = json.load(open(f"{YTREE}/current_tree.json"))
def snps_of(node):
    out = []
    for tok in (node.get("snps") or "").split(","):
        for name in tok.strip().split("/"):
            name = name.strip()
            if name in idx: out.append((name,) + idx[name])
    return out
# 2. pileup at all tree SNP positions (one pass): collect positions
positions = set()
def collect(n):
    for s in snps_of(n): positions.add(s[1])
    for c in n.get("children", []): collect(c)
collect(tree)
bed = f"{W}/yfull_positions.bed"
with open(bed, "w") as f:
    for p in sorted(positions): f.write(f"Y\t{p-1}\t{p}\n")
print("tree SNP positions:", len(positions), file=sys.stderr)
PILE = f"{W}/y_pileup.tsv"
if not os.path.exists(PILE):
    # samtools mpileup over a sorted BED is a single sequential pass (much faster than bcftools -R random access)
    subprocess.run(f"samtools mpileup -f {REF} -l {bed} -q 20 -Q 20 -d 500 {CRAM} 2>/dev/null > {PILE}", shell=True, check=True)
q = open(PILE).read()
import re
counts = {}
for l in q.splitlines():
    f = l.split("\t")
    if len(f) < 5: continue
    pos, ref, bases = int(f[1]), f[2].upper(), f[4]
    bases = re.sub(r"\^.", "", bases); bases = bases.replace("$", "")
    bases = re.sub(r"[+-](\d+)(?=[ACGTNacgtn*#]+)", lambda m: "\x00" * 0, bases)  # strip indel length markers
    # remove indel sequences: after +N or -N, N bases follow
    out = []; i = 0
    b = f[4]
    while i < len(b):
        ch = b[i]
        if ch == "^": i += 2; continue
        if ch == "$": i += 1; continue
        if ch in "+-":
            j = i + 1; num = ""
            while j < len(b) and b[j].isdigit(): num += b[j]; j += 1
            i = j + int(num or 0); continue
        out.append(ch); i += 1
    c = collections.Counter()
    for ch in out:
        if ch in ".,": c[ref] += 1
        elif ch.upper() in "ACGT": c[ch.upper()] += 1
    counts[pos] = dict(c)
refbase = {}
for l in q.splitlines():
    f = l.split("\t")
    if len(f) >= 3: refbase[int(f[1])] = f[2].upper()
def state(pos, anc, der):
    """Within haplogroup O the hg19 reference (R1b) carries the ancestral allele at every branch SNP, so
    derived := non-reference allele; ybrowse anc/der labels are only used when they agree with this."""
    c = counts.get(pos, {}); dp = sum(c.values())
    if dp < 2: return "nocov", dp, c
    rb = refbase.get(pos)
    if rb and rb in (anc, der):
        anc, der = rb, (der if der != rb else anc)
    nd, na = c.get(der, 0), c.get(anc, 0)
    if nd >= max(2, 0.8 * dp): return "der", dp, c
    if na >= max(2, 0.8 * dp): return "anc", dp, c
    return "mixed", dp, c
def score(node):
    st = collections.Counter(state(p, a, d)[0] for _, p, a, d in snps_of(node))
    return st["der"], st["anc"], st["nocov"] + st["mixed"]
# 3. greedy descent
START = "O-F438"   # terminal branch from the ISOGG-2016 descent (rounds 1-4); orientation is reliable inside O
def find_node(n, i):
    if n["id"] == i: return n
    for c in n.get("children", []):
        r = find_node(c, i)
        if r: return r
node = find_node(tree, START); d0, a0, o0 = score(node)
path = [(node["id"], d0, a0, o0, node.get("formed"), node.get("tmrca"))]
# sanity: also score a few O trunk nodes
for chk in ["O", "O-M175", "O2", "O-M122", "O-M117", "O-F8"]:
    nn = find_node(tree, chk)
    if nn: print(chk, score(nn), file=sys.stderr)
while True:
    kids = [c for c in node.get("children", []) if not c["id"].endswith("*")]
    best = None
    for c in kids:
        d, a, o = score(c)
        if d > a and d >= 1 and (best is None or (d - a) > (best[1] - best[2])): best = (c, d, a, o)
    if best is None: break
    node = best[0]; path.append((node["id"], best[1], best[2], best[3], node.get("formed"), node.get("tmrca")))
lines = [f"YFull tree {open(f'{P}/data/ref/ytree/current_version.txt').read().strip()} (descent started at {START}); terminal branch: {path[-1][0]}",
         f"  formed ~{path[-1][4]} ybp, TMRCA ~{path[-1][5]} ybp", "", "Path (branch, #derived, #ancestral, #untyped/mixed, formed, tmrca):"]
for p in path: lines.append("  %-28s der=%3d anc=%3d n/a=%3d  formed=%s tmrca=%s" % p)
# children of terminal: show why we stopped
lines.append(""); lines.append("Children of terminal branch (not supported):")
for c in node.get("children", []):
    d, a, o = score(c); lines.append(f"  {c['id']:28s} der={d} anc={a} n/a={o}")
# 4. detail table for last 6 branches
rows = []
for bid, *_ in path[-6:]:
    def find(n):
        if n["id"] == bid: return n
        for c in n.get("children", []):
            r = find(c)
            if r: return r
    for name, pos, a, d in snps_of(find(tree)):
        s, dp, c = state(pos, a, d); rows.append((bid, name, pos, a, d, s, dp, c.get(a, 0), c.get(d, 0)))
pd.DataFrame(rows, columns=["branch", "snp", "pos_hg19", "anc", "der", "state", "depth", "n_anc", "n_der"]).to_csv(f"{W}/y_terminal_snps.tsv", sep="\t", index=False)
open(f"{W}/y_haplogroup_yfull.txt", "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
