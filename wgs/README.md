# Genome Atlas WGS

把一份个人全基因组测序数据（FASTQ / BAM / CRAM，GRCh37）变成一份经过审计的、中英双语的单文件 HTML 报告。

这是 [`genome-atlas`](../genome-atlas-skill) 的姊妹项目。那一个处理消费级芯片导出的 15 MB 基因型表，
这一个处理原始测序数据。**两者的约定不通用**，不要混着用。

## 它做什么

质控与可调用掩膜 · Y 与 mtDNA 单倍群（含异质性与线粒体疾病位点）· 1000G 祖源 PCA ·
局部祖源与留出校准 · 古代 DNA 投影 · ClinVar 与 ACMG 筛查 · 星号等位基因药物基因组（含拷贝数）·
HLA 分型与疾病/药物关联 · 多基因评分 · 结构变异与拷贝数 · 重复扩增 · 古老人类渗入片段 ·
十四个血型系统 · KIR · 读段与统计相位 · 血液体细胞信号 · 96 类突变谱

## 它不做什么

不做医学诊断，不给用药剂量，不推断性格、智力或族群归属。报告里每条结论都带证据等级和限制说明。

## 快速开始

```bash
cp config.example.yaml config.yaml     # 填 sample_id / sex / reads 路径 / 参考人群
bash setup/00_install_tools.sh
bash setup/02_download_references.sh   # 约 60 GB
bash setup/03_build_prs_reference.sh
python3 setup/04_fetch_scores.py
bash run_all.sh                        # 半天到一天
# 然后按 SKILL.md 第三步写 work/report_text.yaml，再渲染：
python3 scripts/31_html_report.py
```

`config.yaml` 与 `work/` 都在 `.gitignore` 里。**个人数据只出现在 work/。**

## 给 AI agent 用

`SKILL.md` 是给 agent 的完整契约：跑什么、按什么顺序、怎么写结论、哪些话不能写、踩过哪些坑。
先读它，再动手。

## 目录

```
SKILL.md              agent 契约（先读这个）
config.example.yaml   配置模板
setup/                工具安装、参考数据下载、评分参照集
scripts/              编号脚本，按顺序运行
panel/                可编辑的位点面板（改这里就能改分析范围）
templates/            报告模板
docs/                 方法学、审计清单、证据等级
example/              示例报告与文案（只作格式参考，不要沿用结论）
```

## 许可

代码 MIT。`templates/` 沿用 Lieflat Charts 的报告版式，按 PolyForm Noncommercial 授权。
参考数据各有各的引用要求，见 `THIRD_PARTY_NOTICES.md`。
