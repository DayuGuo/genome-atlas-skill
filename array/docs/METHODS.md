# Genome Atlas · 方法学与流程（以示例样本的数字为例）

版本：2026-09-08（含第 1–3 轮分析与第 4 轮独立审计、填补重算）
示例数字来自 `example/`（一位东亚样本的 WeGene 芯片）；你的运行会生成自己的数字。
本文档只记录**怎么做的**；结果见 `example/report.html`。

---

## 0. 一页总览

```
work/data/input.txt (WeGene 芯片导出, 596,744 位点, GRCh37)
   │
   ├─ 01 转换 ──────────► parquet / duckdb / 23andMe 文本 / 基本 QC
   ├─ 02 Y 单倍群 ──────► ISOGG 树贪心下行 (+ YFull 9.05 交叉核对)
   ├─ 03 mtDNA ─────────► HSD → haplogrep3 (PhyloTree 17.2)
   ├─ 04 性状/PGx 面板 ─► Ensembl GRCh37 校验等位基因方向
   ├─ 05 祖源 PCA ──────► 1000G 建空间, 样本投影 (全球 + 东亚)
   ├─ 06 ClinVar ───────► 坐标 + ALT 等位基因匹配 (SNV)
   ├─ 08 ROH ───────────► 自写扫描, 着丝粒断开
   ├─ 09 PRS (子集) ────► PGS Catalog, 仅直接分型位点  ← 审计判定不可作为正式结果
   ├─ 10 HGDP ──────────► hg19→hg38 liftover, 东亚 223 人 PCA 投影
   ├─ 13 探索 ──────────► 染色体南北指数 / X 祖源 / 低频等位基因
   │
   ├─ 审计 ─────────────► 对齐统计, ClinVar 一致性, 逐条结论/性状/PGx 判定, PRS 子集误差实测
   ├─ 14 填补 ──────────► Beagle 5.4 + 1000G phase3 (hg19)
   ├─ 15 PRS (填补后) ──► 参考 MAF 过滤, 剂量打分, 同一位点集给 504 名东亚参考打分
   ├─ 16 填补自检 ──────► 遮蔽 5% 已分型位点回填
   └─ 17/18 报告 v2 ────► 全部数字来自 work/results_v2/report_data_v2.json
```

---

## 1. 输入数据

| 项 | 值 |
|---|---|
| 文件 | `work/data/input.txt`（配置项 input），WeGene 导出，文件头时间戳 2026-09-08 |
| 格式 | 制表符分隔：rsid、染色体、位置、基因型（两个字母，含 I/D 表示插入/缺失，`--` 为未检出） |
| 位点数 | 596,744（578,154 个 rsID + 18,590 个 WeGene 内部 `w` 编号；末尾 3 行空记录已过滤） |
| 坐标系 | GRCh37 / hg19，正链（见 §9 审计：505,634 个可判定位点零链翻转） |
| 染色体 | 1–22、X、Y、MT；MT 坐标与 rCRS 一致 |
| 缺失 | 444 个（0.07%），chr6 的 HLA 区最多（185） |
| 性质 | 基因分型芯片数据，**不是**测序数据；不含探针强度，无法做 CNV |

---

## 2. 计算环境

| 项 | 值 |
|---|---|
| 机器 | Linux 6.17，ARM64（aarch64），20 核，119 GB 内存 |
| Python | 3.12；pandas 2.3.3、pyarrow、duckdb、scipy、matplotlib、requests、pyliftover |
| plink2 | v2.0.0-b.1（2026-09-07）；无 ARM64 预编译版，从 `plink-ng` 源码编译，链接本地 conda 前缀的 OpenBLAS/LAPACK（`tools/blas`）。二进制 `tools/plink2` |
| Beagle | 5.4（22Jul22.46e），Java 8+ |
| haplogrep3 | 3.2.2 |
| yhaplo | 23andMe GitHub 版（仅用于生成 ISOGG 位点表与树文件） |
| Ensembl REST | `grch37.rest.ensembl.org`（variation、VEP region），查询日 2026-09-08 |

参考数据（`work/data/ref`，约 30 GB，均可重新下载）：

| 资源 | 版本 | Build | 用途 |
|---|---|---|---|
| 1000 Genomes phase 3，plink2 bundle（cog-genomics） | 20130502 callset | GRCh37 | PCA、等位基因频率、PRS 参考 |
| 1000 Genomes phase 3 VCF（EBI FTP，chr1–22） | 20130502 v5b | GRCh37 | Beagle 填补参考（chr8–22 去重，见 §10.1） |
| HGDP（Bergström 2020），plink2 bundle | 2019-05-16 | **GRCh38** | 东亚精细 PCA |
| UCSC hg19ToHg38.over.chain | — | — | liftover |
| ClinVar `clinvar.vcf.gz` | fileDate 2026-09-05 | GRCh37 | 变异筛查 |
| PGS Catalog harmonized `*_hmPOS_GRCh37.txt.gz` + 全目录元数据 | 下载 2026-09-08 | GRCh37 | PRS |
| ISOGG SNP index（yhaplo 内置） | 2016.01.04 | GRCh37 | Y 树 |
| YFull YTree | 9.05.0 | — | Y 命名交叉核对 |
| PhyloTree（haplogrep3 内置） | 17.2 | rCRS | mtDNA |
| rCRS NC_012920.1 | NCBI | — | mtDNA 参考序列 |
| Beagle 遗传图谱 `plink.GRCh37.map` | — | GRCh37 | 填补 |

---

## 3. 数据转换与 QC（`scripts/01_convert.py`，`scripts/audit_harmonization.py`）

1. 读入原始文件，拆分基因型为两个等位基因，标记缺失与杂合。
2. 输出 `data/genome.parquet`、`data/genome.duckdb`（表 `snp`，rsid 索引）、`data/genome.23andme.txt`。
3. QC 指标：总位点、rsID/内部 ID 数、重复 rsID、重复坐标、字母表、检出率、按染色体的缺失数、常染色体杂合率（28.9%）、X/Y/MT 杂合数。
4. 性别一致性：Y 染色体 18,943/18,956 位点有效检出（99.9%）；X 零杂合。**注意**：X 的 PAR1/PAR2 区 307 个位点也全部纯合，男性在 PAR 应为二倍体，说明 WeGene 对男性 X 强制输出纯合；因此性别推断主要依据 Y 检出，X 零杂合只是一致而非独立证据。
5. 交互查询：`scripts/query.py`（按 rsid / 区域 / 基因 / SQL，自动附 Ensembl 注释）。

---

## 4. 父系与母系单倍群

### 4.1 Y 染色体（`scripts/02_y_haplogroup.py`）
- 输入：Y 染色体 18,943 个有效基因型（取第一个等位基因作为半合子）。
- 位点表与树：yhaplo 生成的 ISOGG 2016.01.04 清洗版位点表与 YCC 树（`work/results/yhaplo/`）。yhaplo 自身的搜索在本数据上停在 BT，故自写调用器。
- 算法：从根开始贪心下行；每个节点对其子节点统计衍生态/祖先态位点数，选择衍生态占多数且 `der − anc` 最大的子节点；无子节点满足则停止。节点别名（如 `NO/K2a`）按 `/` 拆分合并计数。
- 输出：`work/results/y_haplogroup.txt`（路径、各级支持数、终端子节点证据）、`work/results/y_isogg_marker_states.tsv`。
- 交叉核对：下载 YFull YTree 9.05.0 JSON，按精确 SNP 名定位各级节点，确认路径 O-M122 > M324 > P201 > P164 > M134 > PAGE23(=M117) > M1706(F5/F8) > A9459 > F438。
- 报告口径：以 SNP 命名 O-F438 为主；披露终端仅 1 个位点支持、路径上零星祖先态冲突数。

### 4.2 线粒体（`scripts/03_mt_haplogroup.py`）
- 从 NCBI 取 rCRS；对 4,438 个芯片 MT 位点（排除 I/D）与 rCRS 比对，得到 36 个差异。
- 生成 HSD 文件：`Range` 字段列出**实际覆盖的位置区间**，使 haplogrep3 只在覆盖范围内评估缺失的定义位点。
- `haplogrep3 classify --tree phylotree-rcrs@17.2 --extend-report`。
- 输出：`work/results/mt_haplogrep3/haplogroups.txt`、`work/results/mt_sites_vs_rCRS.tsv`。质量分 0.945 是在芯片覆盖范围（约 27% 的 rCRS）内的评分。

---

## 5. 祖源分析

### 5.1 1000G 全球与东亚 PCA（`scripts/05_ancestry_pca.sh`，`07_ancestry_plot.py`）
1. 取样本常染色体 ACGT 位点（548,438），用 `--extract range` 从 1000G bundle 抽取同坐标的双等位 SNP（522,418）。
2. 在 Python 中按 1000G 的 REF/ALT 重新编码 样本基因型为 VCF：基因型 ⊆ {REF,ALT} 直接编码；否则尝试互补；A/T、C/G 歧义位点丢弃（16,784）；冲突位点丢弃（实际为 0）。得到 505,634 个位点。
3. plink2 不支持不同样本集合并（`--pmerge` 非拼接模式未实现），因此采用**投影**：
   - 参考：`--maf 0.01 --indep-pairwise 200 50 0.2` → 136,231 个位点；`--freq --pca 10 allele-wts` 只在 2,504 名参考样本上估计特征向量。
   - 投影：对参考与样本分别 `--score <eigenvec.allele> 2 5 header-read no-mean-imputation variance-standardize --read-freq <afreq>`。样本不参与特征向量估计。
4. 东亚内部：仅 EAS 504 人重新修剪（87,019 位点）、重新 PCA 并投影。
5. 最近群体度量：**特征值加权欧氏距离**（PC1–6 各乘 √特征值）到各群体质心，另报告 15 个最近个体的群体构成。第一轮使用的按群体标准差标准化的距离对内部离散度大的群体有偏，审计后统一替换。
6. 未做：长程 LD 区域（MHC、倒位区）剔除；预期不改变结论。

### 5.2 HGDP 东亚精细 PCA（`scripts/10_hgdp_pca.py`）
1. pyliftover + UCSC chain 将样本常染色体位点 hg19 → hg38（548,139/548,438 成功；要求映射后染色体不变，多重映射取第一个）。
2. 抽取 HGDP `region == EAST_ASIA` 的 223 人（18 个群体）在这些坐标的双等位 SNP（513,417），同 §5.1 方法编码样本（497,363 位点，0 冲突）。
3. LD 修剪（MAF ≥ 0.02）→ 83,679 位点；PCA + 投影同上。
4. 输出最近质心与 kNN；图 `work/results/hgdp_pca.png`。
5. 口径：结论限定为“在当前参考面板中位于东亚/汉族连续体的偏北侧”；NorthernHan 仅 10 人。

### 5.3 探索性祖源（`scripts/13_fun.py`）
- **染色体级南北指数**：对每条常染色体单独用东亚 PCA 的 PC1 权重打分（`--chr N`），样本与参考同法；指数 = (x − mean_CHS)/(mean_CHB − mean_CHS)。它是坐标，不是祖源比例；报告同时给出相对 CHB 分布的 SD 位置。
- **X 染色体**：样本半合子写成纯合 VCF；1000G X 位点抽取后，因 plink2 PCA 只接受常染色体，用 VCF 往返把染色体名改为 `1` 做 PCA（4,044 个 LD 独立位点）。PAR 未剔除；只支持大洲级判断。

---

## 6. 单位点性状与药物基因组（`scripts/04_trait_panel.py`，`scripts/trait_panel.tsv`）

1. 面板为手工维护的 TSV：rsid、基因、类别、性状、效应等位基因（正链）、0/1/2 拷贝的解读。
2. 从 DuckDB 取基因型；批量调用 Ensembl GRCh37 `variation` 接口取 allele_string、MAF、ancestral、clinical_significance，缓存到 `work/results/ensembl_grch37_variation_cache.json`。
3. 链方向校验：效应等位基因 ∈ allele_string 记 `ok`；其互补 ∈ allele_string 记 `FLIP?`；否则 `mismatch`。对 A/T、C/G 位点该检验无判别力，依赖 §9 的整体正链结论。
4. X 连锁位点按半合子解读。
5. 审计补充（`scripts/audit_traits_pgx.py` + Ensembl VEP region 查询）：
   - TAS2R38 用三个定义位点（rs713598/rs1726866/rs10246939）核实氨基酸（Ala49/Val262/Ile296）后判定单倍型；
   - ABO 用五个位点（rs8176719 缺失、rs8176746/747、rs7853989 B 标记、rs41302905 O2、rs8176750 A2 标签）按常见等位基因模型预测；
   - ACE 直接使用芯片上的 I/D 插缺位点 rs1799752；
   - PGx 按 CPIC 口径只陈述“检测到/未检测到的等位基因”，不指定芯片无法支持的代谢表型；DPYD rs1801159 归为 *5 正常功能；UGT1A1 *28 用标签 rs887829；HLA 用标签位点并明示非分型。
6. 每条结果标注证据等级（A 直接分型且解读明确 / B 多位点或强标签推断 / C 人群统计 / D 探索性）与来源类型（DIRECT / TAG / INDEL / HAPLOTYPE / IMPUTED / PRS）。

---

## 7. ClinVar 筛查（`scripts/06_clinvar_screen.py`）

- 逐行读取 ClinVar GRCh37 VCF；仅 SNV（REF、ALT 单碱基）；按 (染色体, 位置) 匹配芯片有效基因型；要求 ALT 出现在基因型中。
- 输出 `work/results/clinvar_screen.tsv`：分类、评审状态、基因、疾病、HGVS；按 P/LP → 风险因素/药物反应 → 关联 → 冲突 → VUS → 良性排序。
- 审计补验：7,494 行基因型全部 ⊆ {REF, ALT}。
- 口径：“芯片直接覆盖且通过 QC 的 SNV 中未发现带评审标准的 P/LP”；未覆盖的罕见变异、indel、CNV、SV 不能排除。

---

## 8. 纯合片段 ROH（`scripts/08_roh.py`）

- 常染色体 ACGT 位点按坐标排序；连续纯合 run，允许 1 个杂合 call；相邻位点间隔 > 500 kb 即断开（着丝粒与组装缺口）；保留 ≥ 100 个 SNP 且 ≥ 1.5 Mb 的片段。
- 芯片密度约 190 SNP/Mb，参数与之匹配；不做 LD 修剪。
- 口径：“未观察到提示近期近亲婚配的长 ROH”。

---

## 9. 审计：等位基因对齐与结论复核（`work/audit/`）

| 步骤 | 脚本 | 输出 |
|---|---|---|
| 原始 QC、与 1000G 对齐分类（exact / 需互补 / 歧义 / 冲突 / 无记录）、PAR 检查、ClinVar 等位基因一致性、数据库日期 | `scripts/audit_harmonization.py` | `work/audit/qc_and_harmonization.json`，`work/audit/allele_harmonization.tsv` |
| rsID 与 1000G 同坐标一致性（99.76%） | 内联脚本 | `work/results_v2/rsid_mismatch_vs_1000G.tsv` |
| 性状逐条判定、PGx 逐条判定 | `scripts/audit_traits_pgx.py` | `work/audit/05_trait_validation.tsv`，`work/audit/03_pgx_validation.tsv` |
| 逐条结论审计（KEEP / KEEP_WITH_CAVEAT / REWRITE / RECALCULATE / REMOVE） | 内联脚本 | `work/audit/02_claim_audit.tsv` |
| PRS 子集误差实测 | `scripts/audit_prs_validation.py` | `work/audit/04_prs_validation.tsv` |
| 文档 | — | `work/audit/00_executive_summary.md`、`01_pipeline_audit.md`、`06_reanalysis_plan.md`、`07_corrected_report.md`、`08_change_log.md` |

核心数字：505,634 个可判定位点零链翻转、零等位基因冲突；16,784 个歧义位点在子集分析中保守丢弃，在填补输入中纳入（因整体已证明正链）。

---

## 10. 基因型填补与 PRS 重算

### 10.1 填补（`scripts/14_impute.sh`）
- 输入：样本基因型按 1000G REF/ALT 重编码的 VCF，22 条常染色体，522,418 位点（含歧义位点）。
- 参考：1000G phase 3 hg19 VCF；chr8–22 的原始文件存在重复标记（同坐标同 REF/ALT），Beagle 拒绝读取，已用 `awk '!seen[pos:ref:alt]++'` 去重（`scripts/dedup_ref.sh`）。
- 命令：`java -Xmx48g -jar beagle.jar gt=target.chrN.vcf ref=ALL.chrN.vcf.gz map=plink.chrN.GRCh37.map gp=true nthreads=20`。
- 输出：`work/results_v2/imputed/chrN.vcf.gz`（含 DS 剂量与 GP 后验）。

### 10.2 填补自检（`scripts/16_impute_accuracy.sh`）
- chr20、chr22 各随机遮蔽 5% 已分型位点（固定种子），重新填补，与芯片基因型比较。
- 1,079 个位点总体一致率 94.4%，杂合位点 91.9%；GP ≥ 0.99 的 872 个位点 98.6%。这是常见位点上的乐观估计。

### 10.3 关键发现：单样本 Beagle 的 DR2/AF 不能用作过滤
输出 INFO 中的 AF 与 DR2 是按**目标样本**（仅目标样本）估计的：纯合参考位点 DR2 ≈ 0，杂合/纯合突变位点 DR2 高。若按 DR2 ≥ 0.8 过滤，会系统性保留携带突变的位点、丢弃参考纯合位点，使所有评分虚高。因此：
- 变异筛选**只用参考人群指标**：1000G 东亚 504 人 MAF ≥ 0.01（`FV/kg_all.afreq`）；
- 打分用剂量（DS）；
- 另做 样本 最大 GP ≥ 0.9 的敏感性分析（GP 过滤对纯合/杂合比例影响很小：杂合占比 0.33 → 0.32）。

### 10.4 PRS（`scripts/09_prs.py` 子集版；`scripts/15_prs_imputed.py` 填补版）
- 评分选择：从 PGS Catalog 全目录元数据中筛选开发人群为东亚或含东亚的多祖源评分（台湾人体生物资料库 Chen 2025 系列、GIGASTROKE 东亚、Yengo 东亚身高、Zhang 2024 多祖源血脂等），19 个疾病 + 12 个血液指标/精神类。
- 位点匹配：harmonized GRCh37 坐标 + 效应/其他等位基因与 1000G REF/ALT 一致；ID 统一为 `chr:pos`。
- 参考分布：1000G 东亚 504 人（`work/audit/prs_full/kg_all`，真实测序基因型）在**与 样本 相同的位点集**上打分；百分位与 z 值相对这 504 人（另给 CHB+CHS 208 人）。
- 子集版（第一轮）：仅直接分型位点，覆盖率 7–35%。
- 子集误差实测（`audit_prs_validation.py`）：对 504 人分别用完整评分与子集评分，Spearman 0.63–0.96，百分位绝对偏差中位数 4–16、第 90 分位 14–42；按 pgsc_calc 默认最小重叠 0.75，子集版不作为正式 PRS。
- 填补版：覆盖率 70–97%，权重保留 50–97%（哮喘评分缺少高权重位点，权重保留 49.5%）。
- 口径：只有填补版、子集版、GP 敏感性三者方向一致的评分才在正文提及；百分位只相对 504 人参考。

---

## 11. 报告生成

| 版本 | 脚本 | 说明 |
|---|---|---|
| v1（原始，保留） | `12_html_report.py` + `report_i18n.py` | 中英切换；数字部分硬编码（审计指出的缺陷） |
| v2 | `17_build_v2_data.py` → `work/results_v2/report_data_v2.json`；`18_html_report_v2.py` + `report_v2_template.html` + `report_v2_plain.py` | 全部数字由 JSON 注入（占位符替换）；通俗文案；三层置信度摘要；证据等级标签 |

视觉规范：Lieflat Charts（`larashero3-dotcom/lieflat-charts`）的报告模板 R04 版式与 Basics/Lupi 图库语法：结论式标题 + 副标题图例 + 图 + 大写来源行；可数量用 tick rows（一刻度 = 固定单位）；彩色模式下线宽 ×1.8、透明度地板 0.85；图内数字 Inter 800；条形 80–130 ms/根、点阵 12 ms/个的入场错峰；每图底部大写编码注记；色系锁定“青瓷蓝”单色相蓝阶（明度即数据）。深浅色主题通过 CSS token 同时支持。

校验：`node --check` 语法；jsdom 渲染两种语言均无运行时错误。

---

## 12. 可复现性

- 从 `work/data/input.txt`（配置项 input） 与 §2 的参考数据出发，按脚本编号顺序运行即可重建全部结果；`README.md` 列出每个脚本的输入输出。
- 随机性：ROH、PCA、plink2 打分均为确定性；填补自检遮蔽用固定种子 20260908；Beagle 使用默认固定种子。
- 已知的工具坑：
  - plink2 ARM64 需源码编译并提供 cblas；
  - plink2 `--pmerge` 不支持不同样本集合并 → 用投影；
  - plink2 PCA 仅常染色体 → X 需改名；
  - plink2 导入带剂量的多等位 VCF 报错 → `--import-max-alleles 2`；
  - 1000G phase 3 VCF 有重复标记 → 去重；
  - 单样本 Beagle 的 DR2/AF 依赖目标基因型 → 不能用于过滤；
  - yhaplo 在本数据上停在 BT → 自写下行器。

---

## 13. 局限（写进任何结论之前都要记得）

1. 芯片只覆盖约 60 万个预选常见位点：罕见致病变异、indel、CNV、SV、重复扩增、完整 HLA 与 CYP2D6 分型基本不在范围内；“未检出”≠“没有”。
2. 单倍群终端的支持位点少（Y 终端 1 个），下游无法分辨；mtDNA 仅覆盖约 27%。
3. 祖源结论是参考面板内的相对位置，不是籍贯或民族身份；HGDP 子群样本量小（NorthernHan n=10）。
4. PRS：参考仅 504 人；参考为测序基因型而 样本 为填补基因型；评分本身在东亚人群中的预测力有限；百分位不是疾病概率。
5. 单位点性状的 D 级结果在个人层面没有预测价值。
6. 任何用药相关结论在临床使用前需临床级检测确认。

---

## 14. 建议的下一步（详见 `work/audit/06_reanalysis_plan.md`）

1. 抽血复核血脂、尿酸、胆红素、血压、血糖与血型；
2. 向 WeGene 索取原始测序/强度数据；
3. 用填补结果做 HLA 分型（HIBAG，东亚模型）；
4. 家系数据的 IBD 分析；
5. 不再增加行为、智力、运动天赋类单位点。
