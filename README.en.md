# Genome Atlas

Turn one person's genetic data into an **audited, reproducible**, bilingual (中文/English) single-file HTML
report. Two pipelines live here:

| | Input | Directory | Example report |
|---|---|---|---|
| **Array** | WeGene / 23andMe-style export, ~600k sites, 15 MB | [`array/`](array) | [open](https://dayuguo.github.io/genome-atlas-skill/array/example/report.html) |
| **Whole genome** | FASTQ / BAM / CRAM at 30× | [`wgs/`](wgs) | [open](https://dayuguo.github.io/genome-atlas-skill/wgs/example/report.html) |

Landing page: <https://dayuguo.github.io/genome-atlas-skill/>

## How the two relate

Same design, different input, **conventions do not carry over**. The array rule that "not found ≠ absent" and
its subset-score error correction are wrong on whole-genome data; the whole-genome callable mask and its
fill-before-lift-over step mean nothing on an array. Each `SKILL.md` states what its pipeline is and is not for.

## What they do

**Both**: quality control · Y and mtDNA haplogroups · 1000 Genomes ancestry PCA and projection ·
ClinVar overlap · runs of homozygosity · polygenic scores (the reference population is scored on the identical
variant set, so a percentile means something) · a bilingual single-file report

**Whole genome only**: callable mask · local ancestry with held-out calibration · ancient-DNA projection ·
star-allele pharmacogenomics including copy number · HLA typing with disease and drug associations ·
structural variants and copy number · repeat expansions · archaic introgressed segments · read-backed and
statistical phasing · fourteen blood-group systems · KIR · somatic signals in blood · 96-context mutation spectrum

## What they do not do

No diagnosis, no dosing, no inference about personality, intelligence or ethnic origin. Every reading in the
report carries an evidence grade (A/B/C/D) and its limits; anything touching medication carries
"confirm with clinical-grade testing before acting".

## Quick start

```bash
git clone https://github.com/DayuGuo/genome-atlas-skill.git
cd genome-atlas-skill/array     # or cd genome-atlas-skill/wgs
cp config.example.yaml config.yaml
# install tools and download references per that directory's README and SKILL.md, then
bash run_all.sh
```

`config.yaml` and `work/` are git-ignored. **Personal data lives only under work/ and never enters the repository.**

## For AI agents

The `SKILL.md` in each directory is the full contract: what to run, in what order, how to write the
conclusions, what must never be written, and every trap already hit. Read it before touching anything.
The writing rules are not suggestions; they came out of an audit.

## Example reports

Both examples come from one East Asian male's real data, published with his consent. **They are there for
format and tone only** — every conclusion in them belongs to that sample. The repository contains no raw data.

## Licence

Code is MIT. Both `templates/` directories follow the report layout of
[Lieflat Charts](https://github.com/larashero3-dotcom/lieflat-charts) and carry PolyForm Noncommercial.
Reference data each carry their own citation requirements: see [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).
