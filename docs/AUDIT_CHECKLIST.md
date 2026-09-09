# 审计清单 · 交付前逐条回答

把答案写进 `work/audit/00_executive_summary.md`。不能回答的标 `REQUIRES_EXTERNAL_VERIFICATION`，不要猜。

## A. 数据与对齐
1. 输入是芯片导出还是测序？位点数、rsID 数、内部 ID 数、indel 数？
2. `audit_harmonization.py`：exact match / 需互补 / 歧义 / 冲突 各多少？冲突为 0 才能继续。
3. 性别一致性：Y 检出率多少？X 是否被强制纯合（PAR 区是否也零杂合）？
4. rsID 与 1000G 同坐标一致率（应 >99%）。按坐标匹配，不按 rsID。
5. 所有参考数据的版本、build、下载日期是否记录？

## B. 单倍群
6. Y 终端由几个位点支持？路径上有多少祖先态冲突？用 YFull 当前树（GitHub `YFullTeam/YTree`）按精确 SNP 名核对路径；报告以 SNP 命名为主，ISOGG 长名为辅并注明版本。
7. mtDNA 覆盖了 rCRS 的百分之几？haplogrep 质量分？私有变异几个？

## C. 祖源
8. PCA 是否参考建空间、样本投影（样本不参与特征向量估计）？
9. 最近群体用的是特征值加权欧氏距离（不是按群体标准差标准化）？参考子群样本量多大？
10. 结论是否限定为"参考面板内的相对位置"？

## D. 单位点
11. 每个效应等位基因是否通过 Ensembl GRCh37 allele_string 校验？A/T、C/G 位点依赖整体正链结论。
12. TAS2R38 是否用三位点？ABO 是否写"预测"？ACE 是否用 rs1799752？
13. 行为/认知/运动/长寿/睡眠位点是否全部降为 D 级并移入"仅供参考"？

## E. 药物基因组（按 CPIC 口径）
14. 每个基因写的是"检出/未检出的等位基因"还是表型？表型只在主要等位基因全覆盖时才写。
15. DPYD rs1801159 是否被误当可操作变异（它是 *5 正常功能）？
16. HLA 是否标为 TAG 并注明"不是 HLA 分型"？
17. NUDT15 / TPMT 是否一并解读，并加"临床级检测确认、不自行调药"？

## F. ClinVar 与 ROH
18. 匹配是否要求 ALT 出现在基因型中，且基因型 ⊆ {REF, ALT} 100%？
19. 是否写成"覆盖范围内未发现带评审标准的 P/LP；未覆盖不可排除"？
20. ROH 是否写"未见提示近亲婚配的长片段"（不写"父母无亲缘"）？

## G. PRS
21. 子集 vs 完整评分的 Spearman、百分位偏差中位数与第 90 分位是多少？（`04_prs_validation.tsv`）
22. 填补后变异筛选是否只用参考 MAF（绝不用样本自身 DR2/GP）？
23. 填补自检一致率？填补后覆盖率与权重保留？
24. 正文只保留填补版、子集版、GP 敏感性三者方向一致的评分？
25. 百分位是否写明"相对于 N 个参考人"？

## H. 报告
26. 首页是否分三层？每个数字是否来自结果文件？免责声明是否保留？删除了哪些位点、为什么？

## 分项评分（0–10）
Data QC · Allele harmonization · Population genetics · Y · mtDNA · Pharmacogenomics · ClinVar · PRS methodology · Trait interpretation · Statistical calibration · Reproducibility · Scientific wording · Privacy。
总体状态：ROBUST / MOSTLY ROBUST WITH CORRECTIONS / EXPLORATORY / MAJOR REANALYSIS REQUIRED。
