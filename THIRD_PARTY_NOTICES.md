# Third-party notices

Pipeline code in `array/` and `wgs/` (`scripts/`, `setup/`, `panel/`, `docs/`) is MIT — see `LICENSE`.

Both skills' `templates/` follow the report layout and chart grammar of **Lieflat Charts**
(`larashero3-dotcom/lieflat-charts`), which is PolyForm Noncommercial. Those template files therefore carry
that licence, not MIT: see `THIRD_PARTY_LIEFLAT_LICENSE.txt`. Commercial use of the report layout needs
permission from its author.

## Tools the pipelines call (installed by each skill's `setup/`, not vendored here)

| Tool | Licence | Used by | For |
|---|---|---|---|
| plink2 | GPLv3 | both | PCA, scoring, file conversion |
| Beagle 5.4 | GPLv3 | both | imputation (array) / statistical phasing (wgs) |
| haplogrep3 | MIT | both | mtDNA haplogroups |
| yhaplo | Apache-2.0 | array | ISOGG SNP index |
| samtools, bcftools, htslib | MIT/Expat | wgs | reads, variants, pileups |
| mosdepth | MIT | wgs | depth and the callable mask |
| bedtools | MIT | wgs | interval arithmetic |
| FLARE | GPLv3 | wgs | local ancestry |
| Delly | BSD-3 | wgs | structural variants |
| ExpansionHunter | Apache-2.0 / PolyForm | wgs | repeat expansions |
| WhatsHap | MIT | wgs | read-backed phasing |
| PharmCAT | MPL-2.0 | wgs | star alleles |
| Picard | MIT | wgs | lift-over between builds |
| Cyrius, SMNCopyNumberCaller | GPLv3 | wgs | CYP2D6 and SMN copy number |
| T1K | MIT | wgs | HLA and KIR typing |

## Reference and annotation data (downloaded by each skill's `setup/`)

1000 Genomes phase 3 · HGDP 2019 · GRCh37 and GRCh38 assemblies · UCSC lift-over chains ·
ClinVar (NCBI, public domain) · Ensembl release 87 (EMBL-EBI, Apache-2.0) ·
gnomAD v2.1.1 constraint (Broad, CC0) · PGS Catalog scoring files (EBI, terms per score) ·
AADR v66.p1 Human Origins (Harvard Dataverse, CC0 — **cite the AADR release and the original studies**) ·
Sprime introgression calls from Browning et al. 2018 (Mendeley Data) ·
YFull YTree and the ybrowse SNP index · IPD-IMGT/HLA and IPD-KIR · PhyloTree 17.2.

Each carries its own citation requirement. If you publish anything derived from these pipelines, cite the
data sources, not this repository alone.
