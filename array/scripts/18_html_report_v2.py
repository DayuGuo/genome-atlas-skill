#!/usr/bin/env python3
"""Render the bilingual single-file HTML report.
Inputs: WORK/results_v2/report_data_v2.json (numbers, from 17_build_v2_data.py) and WORK/report_text.yaml
(narrative sentences + per-site plain readings, written by the analyst/agent; see example/report_text.yaml).
Generic labels live in templates/ui_generic.json; every number is injected from the data file."""
import json, pathlib, sys, yaml
from config import *
sys.path.insert(0, str(PANEL))
from report_i18n import TRAIT_EN, CAT_EN, PRS_EN
D = json.load(open(RESULTS2 / "report_data_v2.json"))
txt_path = pathlib.Path(CFG.get("report_text", WORK / "report_text.yaml"))
if not txt_path.exists(): raise SystemExit(f"missing {txt_path}: copy example/report_text.yaml, then rewrite every sentence from the result files")
T = yaml.safe_load(open(txt_path))
UI = json.load(open(TEMPLATES / "ui_generic.json")); UI.update(T["text"])
LV = json.load(open(TEMPLATES / "levels.json"))
name_zh, name_en = T.get("name_zh", NAME_ZH), T.get("name_en", NAME_EN)
UI = {k: [v[0].replace("{name}", name_zh), v[1].replace("{name}", name_en)] for k, v in UI.items()}
D["name"] = [name_zh, name_en]; D["y_meta"] = T["y"]; D["mt_lines"] = [T["mt"]["lines_zh"], T["mt"]["lines_en"]]
HTML = open(TEMPLATES / "report_template.html", encoding="utf-8").read()
out = (HTML.replace("__DATA__", json.dumps(D, ensure_ascii=False, separators=(",", ":")))
           .replace("__UI__", json.dumps(UI, ensure_ascii=False)).replace("__TRAIT_EN__", json.dumps(TRAIT_EN, ensure_ascii=False))
           .replace("__CAT_EN__", json.dumps(CAT_EN, ensure_ascii=False)).replace("__PRS_EN__", json.dumps(PRS_EN, ensure_ascii=False))
           .replace("__PLAIN__", json.dumps(T.get("trait_readings", {}), ensure_ascii=False)).replace("__PGX_PLAIN__", json.dumps(T.get("pgx_readings", {}), ensure_ascii=False))
           .replace("__PGX_COV__", json.dumps(T.get("pgx_coverage", {}), ensure_ascii=False))
           .replace("__LEVEL_ZH__", json.dumps(LV["LEVEL_ZH"], ensure_ascii=False)).replace("__LEVEL_EN__", json.dumps(LV["LEVEL_EN"], ensure_ascii=False))
           .replace("__TITLE__", f"{name_zh}基因组图鉴"))
OUT = REPORT / "report.html"; OUT.write_text(out, encoding="utf-8"); print("wrote", OUT, OUT.stat().st_size // 1024, "KB")
