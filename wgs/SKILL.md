---
name: genome-atlas-wgs
description: Audited, reproducible analysis of a personal whole-genome sequencing delivery (FASTQ / BAM / CRAM, GRCh37) into a bilingual (中文/English) single-file HTML report. Covers QC and a callable mask, Y and mtDNA haplogroups with heteroplasmy, 1000G ancestry PCA, local ancestry with a calibrated reference, ancient-DNA projection, ClinVar and ACMG screening, star-allele pharmacogenomics with copy number, HLA typing with disease and drug associations, polygenic scores, structural variants, repeat expansions, archaic introgression, blood groups, KIR, phasing and somatic signals. Use when a user hands you their own raw sequencing reads and asks to analyse, explore or report on their genome. Never for diagnosis or dosing.
---

# Genome Atlas WGS — 个人全基因组测序分析 skill

你（AI agent）的任务是**运行一条已经审计过的流水线，然后诚实地写结论**。脚本负责数字，你负责把数字变成分层可信度的、大众能读懂的句子。不要为了"让报告更有趣"绕过审计规则。

这是 `genome-atlas` 的姊妹 skill。那一个处理消费级芯片导出的 15 MB 基因型表；这一个处理原始测序数据。**两者不要混用**：芯片流水线的很多约定（"没测到就是没有"、子集 PRS 的误差修正）在 WGS 上是错的，反过来也一样。

## 适用与不适用

- 适用：个人全基因组测序交付件。CRAM 或 BAM（已排序、已建索引），或一对 FASTQ（用 `setup/01_align.sh` 先比对）。30× 左右，GRCh37/hg19。
- 也适用：厂商同时给了 VCF 的情况（填 `vendor_vcf`，流水线会用它并自己补齐参考纯合）。
- 不适用：芯片导出（用 `genome-atlas`）；GRCh38 坐标（先转到 GRCh37 或换参考）；低于 10× 的数据（拷贝数、STR、体细胞几步不可靠）；医学诊断、用药剂量决定。
- 输出：`work/report/report.html`（单文件、中英切换）+ `work/wgs/**` 的全部中间表。

## 与芯片版最重要的差别

**厂商交付的 VCF 只列变异位点。"不在 VCF 里"不等于"参考纯合"，也可能是根本没测到。**
所有与参考人群的比较（PCA、PRS、局部祖源、ROH）都必须先用可调用掩膜（默认 MAPQ ≥ 20、深度 ≥ 8）把参考纯合**显式写出来**，掩膜外的位点记为缺失。`03_complete_set.py` 做这件事，后面每一步都吃它的输出。跳过这一步，PCA 会被拉向中心，PRS 会系统性偏低。

其余三条同样重要：

- **X 染色体常被整条按单倍体调用**，包括本该二倍体的拟常染色体区。`02_recall_x_mt.sh` 用正确的倍性重新调用，X 上任何杂合率统计都要用重调结果。
- **PharmCAT 只接受 GRCh38**。必须先在 GRCh37 上填好参考纯合再转换坐标，否则参考等位基因在两版本间相反的位点会被误判（CYP3A5 会从 \*3/\*3 变成 \*1/\*1）。`10_pharmcat.sh` 已经这样做。
- **只用 PASS 变异下结论**。NO_PASS 集合单独查一遍是为了确认过滤器没有丢掉真变异，不是结果。

## 第一步：环境、参考数据、配置

```bash
cp config.example.yaml config.yaml        # 填 sample_id / name / sex / reads 路径 / 参考人群
bash setup/00_install_tools.sh            # conda 环境 + PharmCAT/Picard/Cyrius/FLARE/SMN/haplogrep3
bash setup/02_download_references.sh      # 参考基因组、1000G、ClinVar、Ensembl、YFull、AADR、Sprime ≈ 60 GB
bash setup/03_build_prs_reference.sh      # 用配置里的参考超级人群建评分参照集
python3 setup/04_fetch_scores.py          # 下载 panel/pgs_scores.tsv 里列出的 PGS 评分文件
# FASTQ 起步时先跑：bash setup/01_align.sh
```

`config.yaml` 与 `work/` 已在 `.gitignore`。**个人数据只能出现在 work/ 里。**

## 第二步：按顺序运行（任何一步失败就停下来修，不要绕过）

一键：`bash run_all.sh`（第 30 步会停下等你写文案）。逐步跑时的顺序与产物：

| # | 命令 | 做什么 | 关键输出 |
|---|---|---|---|
| 01 | `bash scripts/01_normalize.sh` | 左对齐、拆多等位、取 PASS、建可调用掩膜 | `wgs/00_input/callable.bed` |
| 02 | `bash scripts/02_recall_x_mt.sh` | 按正确倍性重调 X，MT 全深度 pileup | `wgs/00_input/X.recall.vcf.gz` |
| 03 | `python3 scripts/03_complete_set.py` | **补全基因型集**（参考纯合显式填充） | `wgs/02_complete/*.1kg.pgen` |
| 04 | `bash scripts/04_ancestry_pca.sh` + `04b` + `04c` | 全球与人群内 PCA、逐染色体轴指数 | `wgs/04_ancestry/summary.txt` |
| 05 | `python3 scripts/05_y_haplogroup.py` | Y 单倍群（YFull 树 + CRAM pileup） | `wgs/03_haplo/y_haplogroup_yfull.txt` |
| 06 | `python3 scripts/06_mtdna.py` + `06b` | mtDNA 单倍群、异质性、疾病位点筛查 | `wgs/03_haplo/mt_*.tsv` |
| 07 | `bash scripts/07_annotate.sh` + `07b` + `07c` | ClinVar、后果注释、频率、LoF 表 | `wgs/05_clinvar/clinvar_PLP.tsv` |
| 08–09 | `python3 scripts/08_aadr_extract.py` → `bash 09_aadr_pca.sh` → `09b` | 古代 DNA 面板与投影 | `wgs/11_aadr/near_ancient.tsv` |
| 10–11 | `bash scripts/10_pharmcat.sh`，`python3 scripts/11_pgx_extra.py` | 星号等位基因、CYP2D6 拷贝数、HLA、扩展位点 | `wgs/06_pgx/pharmcat/pharmcat_summary.tsv` |
| 12 | `python3 scripts/12_prs.py` | 多基因评分（覆盖率 96–100%） | `wgs/07_prs/prs_wgs.tsv` |
| 13 | `python3 scripts/13_sv_filter.py` | 结构变异 + 深度重定基因型 | `wgs/08_sv/sv_filtered.tsv` |
| 14–15 | `bash scripts/14_phase_reads.sh` + `14b`，`bash 15_phase_statistical.sh` | 读段相位、统计相位 | `wgs/10_phase/summary.txt` |
| 16–17 | `bash scripts/16_local_ancestry.sh` + `16b`，`17`，`17b` | 局部祖源与**留出校准** | `wgs/12_localanc/calibration.tsv` |
| 18–19 | `python3 scripts/18_archaic.py`，`19_archaic_genes.py` | 古老人类片段与所覆盖基因 | `wgs/15_archaic/gene_families.tsv` |
| 20–22 | `20_mutation_spectrum.py`，`21_somatic.py`，`bash 22_telomere.sh` | 96 类突变谱、克隆性造血、端粒读段 | `wgs/17_mutspec/spectrum96.tsv` |
| 23–27 | `23_bloodgroups.py`，`24_kir.py`，`25_hla_disease.py`，`26_behaviour_prs.py`，`27_candidate_genes.py` | 血型、KIR、HLA 疾病关联、行为评分 | `wgs/16_panels/`, `wgs/19_hla_disease/` |
| 30 | `python3 scripts/30_build_report_data.py` | 汇总为报告数据 | `wgs/report_data.json` |
| 31 | `python3 scripts/31_html_report.py` | 渲染报告 | `work/report/report.html` |

需要多久：01–03 约 1 小时；10、12、16 各半小时到一小时；22 全 CRAM 扫描约 25 分钟；其余分钟级。总计半天到一天。

## 第三步：写结论（脚本不替你写）

复制 `example/report_text.yaml` 到 `work/report_text.yaml`，**逐句重写**。示例是一位东亚样本的真实文案，只用来看格式和语气，不能沿用其结论。每个值是 `[中文, English]`。

必须遵守的写法（来自审计，违反即返工）：

1. **分三层写摘要**：很可靠（读段里直接看得到的）/ 比较可靠（多位点或工具推断的）/ 仅供参考（多基因评分、行为位点、局部祖源百分比）。
2. **"没测到" ≠ "没有"**。凡写"未检出"，必须说明是在哪个范围内没检出，并注明可调用区的覆盖率。PharmCAT 在坐标转换中丢失的位点要单独列出。
3. **PGx 写双倍型，但注明工具与限制**。CYP2D6 的拷贝数来自 Cyrius，HLA 来自短读长分型——准确但不是临床级。凡涉及用药的每一句都带"用药前经临床级检测确认"。
4. **祖源只写参考面板内的相对位置**。局部祖源的百分比**必须与留出校准一起出现**：单独一个"88% 北方"没有意义，"88%，而参考人群里的北京汉平均 70%"才有。同时要给出欧洲/南亚对照面板拿到的比例作为噪声底。
5. **古代 DNA 只写"与哪些古代人群最接近"**，不写族群起源、不写"你是 X 人的后代"。距离随年代的变化趋势可以写，因果解释要留余地。
6. **古老人类片段**：先说总量与随机期望的比较，再说落在哪些基因家族上。片段跨度不是祖源比例，必须写明。
7. **多基因评分**：百分位写明"在参考人群的 N 个人中"。90 以上或 10 以下才值得提。**行为与精神类评分必须标注训练人群**；欧洲训练的评分用在非欧洲样本上时，只能作为对照出现，并解释为什么不能直接用。
8. **行为、认知、性格位点**统一为"部分研究提过，效应小，个人看不出"；不写"你更擅长""你的性格是"。候选基因位点要同时写出流行说法和证据实际支持什么。删除智力、人格、政治倾向、性取向、犯罪相关的推断。
9. **HLA 疾病关联**：阴性结论（乳糜泻的 DQ2/DQ8、发作性睡病的 DQB1\*06:02 等）与阳性同等重要，都要写。阳性写比值比区间与"绝大多数携带者不发病"。
10. **APOE** 写基因型与"人群水平风险因素"，不写致病、不写固定倍数。**ROH** 写"未见提示近亲婚配的长片段"，不写"父母无亲缘"。**体细胞**写检出下限（30× 约对应 10% 变异等位基因比例）。
11. **通俗**：目标读者是没学过遗传学的人。先说这意味着什么，再说可信到什么程度，一到两句。

## 第四步：自检后再交付

- 无头浏览器渲染 `report.html`，两种语言都不报错；`Uncaught` 为零。
- 页面里每个数字都能在 `work/wgs/**` 找到来源；不要在 yaml 里手打数字，用占位符。
- 对照 `docs/AUDIT_CHECKLIST.md` 逐条回答，写进 `work/audit_summary.md`。
- 报告末尾保留免责声明。

## 已知陷阱（都踩过）

- **VCF 只含变异位点**：见上文，这是最容易毁掉整份报告的一条。
- **PharmCAT 坐标转换**：33 个位点的参考等位基因在 GRCh37 与 GRCh38 之间相反。必须先填后转。94 个含插入缺失的位点无法双向转换，要在报告里列出（多在 RYR1、G6PD、CYP2D6）。
- **单样本 Beagle 的 DR2/AF 是目标基因型的函数**，绝不能用它过滤位点；只按参考人群 MAF 过滤。
- **1000G phase3 VCF 含重复标记**，Beagle 会中止；先去重。这些 VCF 是普通 gzip 不是 bgzip，`bcftools` 无法直接切片，但 Beagle 与 FLARE 可以直接读。
- **plink2 `--pca allele-wts` 的输出列数随输入格式变化**（bed 输入多一列 `PROVISIONAL_REF?`），`--score` 的列号要相应调整，否则报"mismatching allele codes"。
- **Delly 不接受路径中的空格**，先建软链接。**T1K 不读 CRAM**，要先抽出目标区与未比对读段成 BAM。
- **ybrowse 的祖先/衍生标注对深层主干 SNP 不可靠**。Y 下行从一个已知可靠的支系起算，在该支系内用"非参考 = 衍生"定向。
- **AADR 是 packedancestrymap 的转置格式（TGENO）**：48 字节头 + 每个个体一行、每 SNP 2 比特。解码后务必与 1000G 频率对照验证（相关系数应 > 0.98）。
- **端粒长度可能测不了**。有些交付流程会滤掉低复杂度重复读段。若含 ≥7 个 TTAGGG 重复的读段比例远低于 1e-4，如实写"无法估计"，不要给一个明显不可能的数字。
- **WhatsHap 需要单独的 conda 环境**，且在混合倍性上会失败，只跑常染色体、逐条染色体跑。
- **报告脚本里 `renderAll()` 必须在文件最末尾**，否则在它之后注册的图不会绘制。

## 可编辑的面板

`panel/` 下全是制表符分隔的表，改它们就能改分析范围，不用动脚本：

| 文件 | 内容 |
|---|---|
| `pgs_scores.tsv` | 要计算的多基因评分（疾病 + 行为），含训练人群标注 |
| `hla_disease.tsv` | HLA 疾病与药物关联，`rule` 列是可执行的判定表达式 |
| `pgx_extra.tsv` | PharmCAT 之外的药物基因组位点 |
| `blood_groups.tsv` | 血型系统的定义位点 |
| `candidate_genes.tsv` | 候选行为基因，含"流行说法"与"证据支持什么"两列 |
| `chip_hotspots.tsv` / `chip_genes.tsv` | 克隆性造血热点与基因区间 |
| `mt_disease.tsv` | 线粒体疾病位点 |
| `str_thresholds.tsv` | 重复扩增病的致病阈值 |
| `gene_families.tsv` | 古老人类片段注释用的基因家族正则 |
| `chrom_grch37.tsv` | 染色体长度与着丝粒位置 |

## 视觉规范

报告模板沿用 Lieflat Charts 的报告版式：结论式标题 + 副标题图例 + 图 + 大写来源行；青瓷蓝单色相，祖源与突变谱另用两个经过色觉障碍校验的强调色（`--ancn` / `--ancs` / `--arch`）。改文案改 yaml，改数字改脚本，不要手改 HTML。
