---
name: genome-atlas
description: Audited, reproducible analysis of a consumer genotyping-array export (WeGene / 23andMe-style, GRCh37) into a bilingual (中文/English) single-file HTML report. Covers QC, Y and mtDNA haplogroups, 1000G/HGDP ancestry PCA, evidence-graded trait and pharmacogenomic sites, ClinVar overlap, ROH, Beagle imputation and polygenic scores with measured error. Use when a user hands you a raw SNP-array genotype file and asks to analyse, explore or report on their own genome. Never for diagnosis or dosing.
---

# Genome Atlas — 个人基因组芯片数据分析 skill

你（AI agent）的任务是**运行一条已经审计过的流水线，然后诚实地写结论**。脚本负责数字，你负责把数字变成分层可信度的、大众能读懂的句子。不要跳过审计规则去"让报告更有趣"。

## 适用与不适用

- 适用：WeGene、23andMe 等消费级芯片导出的文本文件（rsid / 染色体 / 位置 / 基因型，GRCh37）。约 50–100 万位点。
- 不适用：FASTQ / BAM / WGS VCF（用专门流程）；GRCh38 坐标的导出（先 liftover 到 GRCh37 或改参考）；医学诊断、用药剂量决定。
- 输出：`work/report/report.html`（单文件、中英切换）+ `work/results*/`、`work/audit/` 的全部中间表。

## 第一步：环境与配置

```bash
cp config.example.yaml config.yaml          # 填 sample_id / name_zh / name_en / sex / input 路径
pip install -r requirements.txt
bash setup/download_references.sh           # 1000G phase3 (plink2 bundle + VCF), HGDP, ClinVar, rCRS, chain, Beagle, haplogrep3 ≈ 45 GB
# plink2: x86_64 直接下载；ARM64 见 setup/build_plink2_arm64.sh
```

`config.yaml` 与 `work/` 已在 `.gitignore`。**个人数据只能出现在 work/ 里。**

## 第二步：按顺序运行（每步失败就停，不要绕过）

| # | 命令 | 做什么 | 关键输出 |
|---|---|---|---|
| 1 | `python3 scripts/01_convert.py` | 转 parquet/duckdb，基本 QC | `work/results/summary_basic.tsv` |
| 2 | `python3 scripts/02_y_haplogroup.py` | Y 单倍群（ISOGG 树贪心下行；首次自动跑 yhaplo 取位点表） | `work/results/y_haplogroup.txt` |
| 3 | `python3 scripts/03_mt_haplogroup.py` | mtDNA 单倍群（haplogrep3，带覆盖范围） | `work/results/mt_haplogrep3/` |
| 4 | `python3 scripts/04_trait_panel.py` | 性状/PGx 位点 + Ensembl 链方向校验 | `work/results/trait_panel_report.tsv` |
| 5 | `bash scripts/05_ancestry_pca.sh` | 1000G PCA，样本投影 | `work/results/ancestry/` |
| 6 | `python3 scripts/06_clinvar_screen.py` | ClinVar 坐标+等位基因匹配 | `work/results/clinvar_screen.tsv` |
| 7 | `python3 scripts/07_ancestry_plot.py` | 参考超级人群内部 PCA + 图 | `work/results/ancestry_summary.txt` |
| 8 | `python3 scripts/08_roh.py` | 纯合片段 | `work/results/roh.tsv` |
| 9 | `python3 scripts/09_prs.py` 与 `PRS_GROUP=biomarker python3 scripts/09_prs.py` | 子集 PRS（**仅作误差对照，不是结果**） | `work/results/prs_*.tsv` |
| 10 | `python3 scripts/10_hgdp_pca.py` | hg19→hg38 liftover，HGDP 区域 PCA | `work/results/hgdp_ancestry_summary.txt` |
| 11 | `python3 scripts/13_fun.py` | 染色体级轴指数、X 祖源、低频等位基因 | `work/results/fun_*.tsv` |
| 12 | `python3 scripts/audit_harmonization.py` | 对齐统计、PAR 检查、ClinVar 一致性 | `work/audit/qc_and_harmonization.json` |
| 13 | `python3 scripts/audit_traits_pgx.py` | 证据等级表、PGx 基因型表 | `work/audit/05_*.tsv`, `03_*.tsv` |
| 14 | `python3 scripts/audit_prs_validation.py` | 子集 vs 完整评分误差实测（~1 h） | `work/audit/04_prs_validation.tsv` |
| 15 | `bash scripts/dedup_ref.sh && bash scripts/14_impute.sh` | Beagle 填补（~30 min，20 线程） | `work/results_v2/imputed/` |
| 16 | `bash scripts/16_impute_accuracy.sh` | 遮蔽自检 | `work/results_v2/impute_check/summary.tsv` |
| 17 | `python3 scripts/15_prs_imputed.py` | 填补后 PRS（参考 MAF 过滤） | `work/results_v2/prs_imputed.tsv` |
| 18 | `python3 scripts/17_build_v2_data.py` | 汇总为报告数据 | `work/results_v2/report_data_v2.json` |
| 19 | **写 `work/report_text.yaml`**（见第三步） | 结论文案 | — |
| 20 | `python3 scripts/18_html_report_v2.py` | 渲染报告 | `work/report/report.html` |

一键：`bash run_all.sh`（第 19 步会停下来等你写文案）。

## 第三步：写结论（这是你的工作，脚本不替你写）

复制 `example/report_text.yaml` 到 `work/report_text.yaml`，**逐句重写**。示例是一位东亚样本的真实文案，只用来看格式和语气，不能沿用其结论。每个值是 `[中文, English]`。

必须遵守的写法（来自审计，违反即返工）：

1. **分三层写摘要**：很可靠（直接分型、意义明确）/ 比较可靠（多位点推断、单倍群、ABO 预测）/ 仅供参考（PRS、行为/运动/长寿位点、低频等位基因、染色体级指数）。
2. **"没测到" ≠ "没有"**。凡写"未检出"，必须带"在测过的位点里"；HLA、CYP2D6、G6PD、罕见变异、CNV 一律注明未覆盖。
3. **PGx 只写等位基因，不写表型**，除非该基因所有主要等位基因都在芯片上（NUDT15、VKORC1、ALDH2 可以）。CYP2C19 缺 *17 → 不能写"代谢正常"。DPYD rs1801159 是 *5 正常功能，不是可操作变异。HLA 标签 SNP 不是 HLA 分型。
4. **单 SNP 不推单倍型**：TAS2R38 要三位点；ABO 只能写"常见等位基因模型预测 X 型，以血清学为准"；ACE 优先用 rs1799752 插缺本身。
5. **祖源只写参考面板内的相对位置**："最近的参考人群中心是 X"，不写籍贯、民族、"证明是 X 人"。轴指数是坐标，不是百分比。
6. **PRS 只用填补后结果**（`prs_imputed.tsv`），子集结果作附录并画实测误差。百分位写明"在 1000G 的 N 个参考人中"，不写"在所有中国人中"。90 以上或 10 以下才值得提；填补版、子集版、GP 敏感性三者方向不一致的评分不进正文。
7. **行为、认知、运动、长寿、睡眠位点**统一为"部分研究提过，效应小，个人看不出"；不写"你更擅长""你的性格是"。删除智力、人格、政治倾向、性取向、犯罪相关位点。
8. **APOE**：写基因型与"人群水平风险因素"，不写致病、不写固定倍数。**MTHFR**：只写基因型与酶活性，不推同型半胱氨酸、不建议补叶酸。**ROH**：写"未见提示近亲婚配的长片段"，不写"父母无亲缘"。
9. **每条读数标证据等级** A/B/C/D（`docs/EVIDENCE_LEVELS.md`）和来源 DIRECT / TAG / INDEL / HAPLOTYPE / IMPUTED / PRS。
10. **通俗**：目标读者是没学过遗传学的人。先说这意味着什么，再说可信到什么程度，一到两句。

`trait_readings` 按 rsid 写观察到的基因型的解读；`pgx_readings` 按基因写；`y` 块填 YFull 交叉核对的路径（`docs/AUDIT_CHECKLIST.md` 有做法）；`mt.lines_*` 填定义位点与私有变异。

## 第四步：自检后再交付

- `node --check` 或用 jsdom 渲染 `report.html` 两种语言无报错。
- 页面里每个数字都能在 `work/results*/` 找到来源；不要在 yaml 里手打数字，用 `{占位符}`（见 `templates/ui_generic.json` 的键）。
- 对照 `docs/AUDIT_CHECKLIST.md` 的 26 个问题逐条回答一遍，写进 `work/audit/00_executive_summary.md`。
- 报告末尾保留免责声明；涉及用药的每一句都带"用药前经临床级检测确认"。

## 已知陷阱（都踩过）

- 单样本 Beagle 输出的 DR2 / AF 是**目标样本基因型的函数**（纯合参考位点 DR2≈0）。**绝不能**按 DR2 过滤填补变异，否则 PRS 系统性虚高；只按参考人群 MAF 过滤。
- 1000G phase 3 VCF 含重复标记，Beagle 会中止；先 `dedup_ref.sh`。
- plink2 PCA 只接受常染色体（X 需改名）；`--pmerge` 不支持不同样本集合并（用投影）；带剂量的多等位 VCF 需 `--import-max-alleles 2`。
- 男性 X 在 WeGene 导出中被强制纯合（PAR 也是），X 零杂合不是独立的性别证据；看 Y 检出率。
- yhaplo 自身的搜索可能停在树根附近；只用它生成位点表，单倍群由 `02_y_haplogroup.py` 判定，并用 YFull 当前树核对命名。
- "最近群体"用特征值加权欧氏距离，不用按群体标准差标准化的距离（后者偏向离散度大的群体）。

## 视觉规范

报告模板遵循 Lieflat Charts（`larashero3-dotcom/lieflat-charts`）的报告版式与图库语法：结论式标题 + 副标题图例 + 图 + 大写来源行；可数量用 tick rows；青瓷蓝单色相；不要改成其他图表库的默认样式。改文案改 yaml，改数字改脚本，不要手改 HTML。
