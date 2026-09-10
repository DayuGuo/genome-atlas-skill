#!/usr/bin/env bash
# Run the whole pipeline in order. Any step that fails stops the run: do not work around a failure,
# fix it. Step 30 stops and waits for you to write work/report_text.yaml.
set -euo pipefail
cd "$(dirname "$0")"
source scripts/env.sh
S=scripts
step(){ echo -e "\n=== $* ==="; }

step 01 normalise the VCF and build the callable mask;         bash  $S/01_normalize.sh
step 02 re-call X with the right ploidy and pile up MT;        bash  $S/02_recall_x_mt.sh
step 03 complete genotype set at the reference panel sites;    python3 $S/03_complete_set.py
step 04 ancestry PCA and projection;                           bash  $S/04_ancestry_pca.sh
step 04b nearest reference populations;                        python3 $S/04b_ancestry_summary.py
step 04c per-chromosome axis index;                            python3 $S/04c_axis_index.py
step 05 Y haplogroup on the YFull tree;                        python3 $S/05_y_haplogroup.py
step 06 mtDNA haplogroup and heteroplasmy;                     python3 $S/06_mtdna.py
step 06b mtDNA disease screen;                                 python3 $S/06b_mt_disease.py
step 07 annotate with ClinVar, consequences, frequencies;      bash  $S/07_annotate.sh
step 07b pathogenic and loss-of-function tables;               python3 $S/07b_clinvar_tables.py
step 07c gnomAD frequencies for the candidates;                python3 $S/07c_gnomad_lookup.py
step 08 extract the ancient-DNA panel;                         python3 $S/08_aadr_extract.py
step 09 ancient-DNA PCA and projection;                        bash  $S/09_aadr_pca.sh
step 09b nearest present-day and ancient groups;               python3 $S/09b_aadr_summary.py
step 10 PharmCAT star alleles, CYP2D6, HLA;                    bash  $S/10_pharmcat.sh
step 11 extended pharmacogenomic markers;                      python3 $S/11_pgx_extra.py
step 12 polygenic scores;                                      python3 $S/12_prs.py
step 13 structural variants and copy number;                   python3 $S/13_sv_filter.py
step 14 read-backed phasing;                                   bash  $S/14_phase_reads.sh
step 14b phase summary and cis/trans questions;                python3 $S/14b_phase_summary.py
step 15 statistical phasing;                                   bash  $S/15_phase_statistical.sh
step 16 local ancestry;                                        bash  $S/16_local_ancestry.sh
step 16b local-ancestry calibration;                           bash  $S/16b_local_ancestry_calibration.sh
step 17 local-ancestry summary;                                python3 $S/17_local_ancestry_summary.py
step 17b calibrated against held-out references;               python3 $S/17b_local_ancestry_calibrated.py
step 18 archaic introgressed segments;                         python3 $S/18_archaic.py
step 19 genes covered by those segments;                       python3 $S/19_archaic_genes.py
step 20 mutation spectrum;                                     python3 $S/20_mutation_spectrum.py
step 21 somatic signals in blood;                              python3 $S/21_somatic.py
step 22 telomere read fraction;                                bash  $S/22_telomere.sh
step 23 blood groups;                                          python3 $S/23_bloodgroups.py
step 24 KIR and HLA ligands;                                   python3 $S/24_kir.py
step 25 HLA disease and drug associations;                     python3 $S/25_hla_disease.py
step 26 behavioural polygenic scores;                          python3 $S/26_behaviour_prs.py
step 27 candidate behaviour genes;                             python3 $S/27_candidate_genes.py
step 30 assemble the report data;                              python3 $S/30_build_report_data.py

cat <<'MSG'

=== now write work/report_text.yaml ===
Copy example/report_text.yaml and rewrite every line for this sample. The scripts produce numbers;
the conclusions are yours to write, under the rules in SKILL.md. Then:

    python3 scripts/31_html_report.py

MSG
