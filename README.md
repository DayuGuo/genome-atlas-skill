# Genome Atlas

把一个人的基因数据变成一份**经过审计的、可复现的**中英双语单文件 HTML 报告。仓库里有两条流水线：

| | 输入 | 目录 | 示例报告 |
|---|---|---|---|
| **芯片版** | WeGene / 23andMe 风格导出，约 60 万位点、15 MB | [`array/`](array) | [看报告](https://dayuguo.github.io/genome-atlas-skill/array/example/report.html) |
| **全基因组版** | FASTQ / BAM / CRAM，30× | [`wgs/`](wgs) | [看报告](https://dayuguo.github.io/genome-atlas-skill/wgs/example/report.html) |

落地页：<https://dayuguo.github.io/genome-atlas-skill/>

## 两者的关系

同一套设计思路，不同的输入，**约定不通用**。芯片版的"没测到不等于没有"、子集评分的误差修正，在全基因组上是错的；
全基因组版的可调用掩膜、坐标转换前先填参考纯合，在芯片上没有意义。各自的 `SKILL.md` 写清了适用与不适用。

## 它们做什么

**两者都有**：质控 · Y 与 mtDNA 单倍群 · 1000 Genomes 祖源主成分与投影 · ClinVar 交叉 · 纯合片段 ·
多基因评分（同位点集给参考人群打分，百分位才有意义）· 中英双语单文件报告

**只有全基因组版**：可调用掩膜 · 局部祖源与留出校准 · 古代 DNA 投影 · 星号等位基因药物基因组（含拷贝数）·
HLA 分型与疾病/药物关联 · 结构变异与拷贝数 · 重复扩增 · 古老人类渗入片段 · 读段与统计相位 ·
十四个血型系统 · KIR · 血液体细胞信号 · 96 类突变谱

## 它们不做什么

不做医学诊断，不给用药剂量，不推断性格、智力或族群归属。报告里每条读数都带证据等级（A/B/C/D）和限制说明；
涉及用药的每一句都带"用药前经临床级检测确认"。

## 快速开始

```bash
git clone https://github.com/DayuGuo/genome-atlas-skill.git
cd genome-atlas-skill/array     # 或 cd genome-atlas-skill/wgs
cp config.example.yaml config.yaml   # 填样本、性别、输入路径、参考人群
# 按各自 README 与 SKILL.md 装工具、下参考数据，然后
bash run_all.sh
```

`config.yaml` 与 `work/` 都在 `.gitignore` 里。**个人数据只出现在 work/ 下，永远不进仓库。**

## 给 AI agent 用

每个目录下的 `SKILL.md` 是给 agent 的完整契约：跑什么、按什么顺序、怎么写结论、哪些话不能写、踩过哪些坑。
先读它，再动手。写结论的规则不是建议，是审计出来的硬约束。

## 示例报告

两份示例来自同一位东亚男性的真实数据，本人同意公开。**只用来看格式和写法**，其中每条结论都属于那个样本，
不要沿用。仓库不含任何原始数据。

## 许可

代码 MIT。两个 `templates/` 沿用 [Lieflat Charts](https://github.com/larashero3-dotcom/lieflat-charts)
的报告版式，按 PolyForm Noncommercial 授权。参考数据各有引用要求，见 [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md)。
