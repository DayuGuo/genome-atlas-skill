"""Shared configuration for the genome-atlas pipeline.
Reads config.yaml at the repository root (copy config.example.yaml). Nothing personal lives in the repo:
raw data, references, tools and every output go under WORK (default ./work, git-ignored)."""
import os, pathlib
try:
    import yaml
except ImportError:  # minimal fallback: key: value lines
    yaml = None
ROOT = pathlib.Path(__file__).resolve().parents[1]
def _load():
    f = ROOT / "config.yaml"
    if not f.exists(): return {}
    if yaml: return yaml.safe_load(open(f)) or {}
    out = {}
    for line in open(f):
        line = line.split("#")[0].strip()
        if ":" in line:
            k, v = line.split(":", 1); out[k.strip()] = v.strip().strip('"').strip("'")
    return out
CFG = _load()
SAMPLE = str(CFG.get("sample_id", os.environ.get("SAMPLE_ID", "SAMPLE")))           # ID used inside VCF/plink files
NAME_ZH = str(CFG.get("name_zh", SAMPLE)); NAME_EN = str(CFG.get("name_en", SAMPLE))  # display names in the report
SEX = str(CFG.get("sex", "male")).lower()                                             # 'male' | 'female'
WORK = pathlib.Path(CFG.get("work_dir", ROOT / "work")).expanduser().resolve()
RAW = pathlib.Path(CFG.get("input", WORK / "data" / "input.txt")).expanduser().resolve()
TOOLS = pathlib.Path(CFG.get("tools_dir", WORK / "tools")).expanduser().resolve()
REF = pathlib.Path(CFG.get("ref_dir", WORK / "data" / "ref")).expanduser().resolve()
PLINK2 = str(CFG.get("plink2", TOOLS / "plink2"))
HAPLOGREP3 = str(CFG.get("haplogrep3", TOOLS / "haplogrep3"))
THREADS = str(CFG.get("threads", 8))
DATA = WORK / "data"; RESULTS = WORK / "results"; RESULTS2 = WORK / "results_v2"; AUDIT = WORK / "audit"; REPORT = WORK / "report"; PGS = WORK / "data" / "pgs"
for d in (DATA, RESULTS, RESULTS2, AUDIT, REPORT, PGS): d.mkdir(parents=True, exist_ok=True)
PANEL = ROOT / "panel"; TEMPLATES = ROOT / "templates"
# ---- reference populations (1000 Genomes labels). Defaults are East Asian; change for other ancestries.
REF_SUPERPOP = str(CFG.get("ref_superpop", "EAS"))                       # super-population used for the within-group PCA and PRS percentiles
SUBPOPS = list(CFG.get("ref_subpops", ["CHB", "CHS"]))                   # sub-populations reported as a second percentile column
AXIS = list(CFG.get("axis_pops", ["CHS", "CHB"]))                        # [pop at 0, pop at 1] for the per-chromosome axis index
HGDP_REGION = str(CFG.get("hgdp_region", "EAST_ASIA"))                   # HGDP 'region' column value
HGDP_FOCUS = list(CFG.get("hgdp_focus", ["Han", "NorthernHan"]))         # HGDP groups drawn in the dark colour
HGDP_SECOND = list(CFG.get("hgdp_second", ["Japanese"]))                 # HGDP groups drawn in the mid colour
