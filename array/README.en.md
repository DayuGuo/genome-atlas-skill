# Genome Atlas · an agent skill for consumer genotyping-array data

[中文](README.md)

Turns a consumer genotyping export (WeGene / 23andMe format, GRCh37) into an **audited, reproducible, bilingual
(中文/English)** single-file HTML report. Built for AI coding tools (Claude Code, Codex): `SKILL.md` tells the
agent how to run the pipeline, how to write conclusions, and what it must not say.

Example report: **https://dayuguo.github.io/genome-atlas-skill/** (live preview; source [`example/report.html`](example/report.html), an East Asian sample, published with consent, name and sex removed).

## What it does

| Module | Method | Output |
|---|---|---|
| QC & harmonization | call rate, heterozygosity, sex consistency, site-by-site strand/allele check against 1000G | proceed only with zero flips and zero conflicts |
| Paternal / maternal lines | ISOGG-tree descent cross-checked with YFull; haplogrep3 with coverage ranges | Y / mtDNA haplogroups with supporting-site counts |
| Ancestry | PCA space from 1000G, sample projected; HGDP regional PCA (automatic liftover) | nearest reference groups, kNN |
| Traits / pharmacogenomics | 76-site panel, Ensembl strand check, evidence grades A–D, CPIC conventions | tables + plain-language readings |
| ClinVar | position + allele matching | overlaps within array coverage |
| ROH | homozygosity scan | relatedness hint |
| Polygenic scores | Beagle imputation → reference-MAF filter → scored with 504 reference people on the same variant set; measured subset error | percentiles with error |
| Report | Lieflat-style template, every number from result files, prose written by the analyst under the rules | `work/report/report.html` |

## Quick start

```bash
git clone <this repo> && cd genome-atlas-skill
cp config.example.yaml config.yaml   # sample id, display names, sex, input path
pip install -r requirements.txt
bash setup/download_references.sh    # ≈45 GB; plink2: download an x86_64 binary, or build on ARM64 with setup/build_plink2_arm64.sh
bash run_all.sh                      # stops at step 19 until you write work/report_text.yaml
python3 scripts/18_html_report_v2.py
```

With an AI tool: load this repository as a skill (Claude Code: `.claude/skills/genome-atlas/`; Codex: its skills
directory), then say "analyse my xxx.txt with genome-atlas".

## Rules that matter

- Personal data lives only in `work/` (git-ignored). Never commit `config.yaml`, `work/`, or any genotype file.
- "Not detected" ≠ "absent": the array covers ~600k common sites.
- Pharmacogenomics reports alleles detected/not detected, never a phenotype the array cannot support; any dosing decision needs clinical-grade testing.
- PRS: imputed results only, with measured error; percentiles are relative to the 1000G reference group.
- Behaviour, cognition, athletics, longevity: "reported in some studies, no visible effect on an individual".

See `SKILL.md`, `docs/EVIDENCE_LEVELS.md`, `docs/AUDIT_CHECKLIST.md`, `docs/METHODS.md`.

## License
Code: MIT. `templates/` is adapted from [Lieflat Charts](https://github.com/larashero3-dotcom/lieflat-charts) and
keeps its PolyForm Noncommercial 1.0.0 license (see `THIRD_PARTY_NOTICES.md`). Reference data and tools have their
own licenses and are not redistributed here.

For research, education and personal curiosity only; not a medical diagnosis or medication advice.
