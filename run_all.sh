#!/usr/bin/env bash
# Run the whole pipeline. Stops before rendering if work/report_text.yaml is missing (the agent must write it).
set -euo pipefail
cd "$(dirname "$0")"; S=scripts
python3 $S/01_convert.py && python3 $S/02_y_haplogroup.py && python3 $S/03_mt_haplogroup.py && python3 $S/04_trait_panel.py
bash $S/05_ancestry_pca.sh && python3 $S/06_clinvar_screen.py && python3 $S/07_ancestry_plot.py && python3 $S/08_roh.py
python3 $S/09_prs.py && PRS_GROUP=biomarker python3 $S/09_prs.py
python3 $S/10_hgdp_pca.py && python3 $S/13_fun.py
python3 $S/audit_harmonization.py && python3 $S/audit_traits_pgx.py && python3 $S/audit_prs_validation.py
bash $S/dedup_ref.sh && bash $S/14_impute.sh && bash $S/16_impute_accuracy.sh && python3 $S/15_prs_imputed.py
python3 $S/17_build_v2_data.py
source $S/env.sh
if [ ! -s "$WORK/report_text.yaml" ]; then echo "Now write $WORK/report_text.yaml (start from example/report_text.yaml, rewrite every sentence from the result files), then run: python3 $S/18_html_report_v2.py"; exit 0; fi
python3 $S/18_html_report_v2.py
