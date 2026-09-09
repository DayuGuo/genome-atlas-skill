# Genome Atlas · 个人基因组芯片数据分析 skill

[English](README.en.md)

把一份消费级基因检测导出文件（WeGene / 23andMe 格式，GRCh37）变成一份**经审计、可复现、中英双语**的单文件 HTML 报告。为 Claude Code、Codex 等 AI 编程工具设计：`SKILL.md` 告诉工具怎么跑、怎么写结论、哪些话不能说。

示例报告：**https://dayuguo.github.io/genome-atlas-skill/**（在线预览；源文件 [`example/report.html`](example/report.html)，一位东亚男性样本，经本人同意公开）。

## 它做什么

| 模块 | 方法 | 输出 |
|---|---|---|
| QC 与对齐 | 检出率、杂合率、性别一致性、与 1000G 逐位点比对链方向与等位基因 | 零翻转零冲突才继续 |
| 父系 / 母系 | ISOGG 树下行 + YFull 核对；haplogrep3（带覆盖范围） | Y / mtDNA 单倍群及支持位点数 |
| 祖源 | 1000G 参考建空间、样本投影；HGDP 区域精细 PCA（自动 liftover） | 最近参考人群、kNN |
| 性状 / 药物基因组 | 76 个位点面板，Ensembl 链校验，证据分级 A–D，CPIC 口径 | 表格 + 通俗解读 |
| ClinVar | 坐标 + 等位基因匹配 | 覆盖范围内的记录 |
| ROH | 纯合片段扫描 | 近亲提示 |
| 多基因评分 | Beagle 填补 → 参考 MAF 过滤 → 与 504 名参考同位点集打分；子集误差实测 | 百分位 + 误差 |
| 报告 | Lieflat 风格模板，所有数字来自结果文件，文案由分析者按规则撰写 | `work/report/report.html` |

## 快速开始

```bash
git clone <this repo> && cd genome-atlas-skill
cp config.example.yaml config.yaml   # 填写样本 ID、显示名、性别、输入文件路径
pip install -r requirements.txt
bash setup/download_references.sh    # ≈45 GB；plink2 x86_64 自行下载，ARM64 见 setup/build_plink2_arm64.sh
bash run_all.sh                      # 跑到第 19 步停下，等你写 work/report_text.yaml
python3 scripts/18_html_report_v2.py
```

用 AI 工具时：把本仓库作为 skill 加载（Claude Code：放入 `.claude/skills/genome-atlas/`；Codex：按其 skills 目录约定），然后说"用 genome-atlas 分析我的 xxx.txt"。工具会按 `SKILL.md` 执行并撰写文案。

## 必读的规则

- **个人数据永远只在 `work/`**（已 gitignore）。不要把 `config.yaml`、`work/` 或含基因型的文件提交。
- "没测到" ≠ "没有"：芯片只覆盖约 60 万个常见位点。
- 药物基因组只写"检出/未检出的等位基因"；任何用药决定需临床级检测确认。
- PRS 只报填补后结果并附实测误差；百分位相对于 1000G 的参考人群。
- 行为、认知、运动、长寿类位点只能写"部分研究提过，个人看不出"。

细则见 `SKILL.md`、`docs/EVIDENCE_LEVELS.md`、`docs/AUDIT_CHECKLIST.md`；方法学见 `docs/METHODS.md`。

## 目录

```
SKILL.md                 给 AI 工具的操作规范（先读这个）
config.example.yaml      配置模板
run_all.sh               一键流水线
scripts/                 01–18 分析脚本 + audit_* 审计脚本 + query.py 交互查询
panel/                   位点面板（性状、证据等级、PGx、PGS 评分列表、英文名）
templates/               报告模板与通用文案（Lieflat 风格，PolyForm Noncommercial）
docs/                    方法学、证据等级、审计清单
setup/                   参考数据下载、ARM64 plink2 编译
example/                 示例报告与示例文案
```

## 许可

代码 MIT。`templates/` 改编自 [Lieflat Charts](https://github.com/larashero3-dotcom/lieflat-charts)，沿用其 PolyForm Noncommercial 1.0.0 许可（见 `THIRD_PARTY_NOTICES.md`）。参考数据与工具各有许可，不在本仓库内分发。

本项目仅供研究、教育与个人探索，不构成医学诊断或用药建议。
