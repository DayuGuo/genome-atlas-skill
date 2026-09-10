# source this from bash scripts: exports SAMPLE, WORK, RAW, TOOLS, REF, PLINK2, THREADS
eval "$(python3 - <<'PY'
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent if '__file__' in dir() else '.'))
PY
)"
_HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
eval "$(python3 -c "import sys; sys.path.insert(0,'$_HERE'); import config as c; print(f'export SAMPLE=\"{c.SAMPLE}\" WORK=\"{c.WORK}\" RAW=\"{c.RAW}\" TOOLS=\"{c.TOOLS}\" REF=\"{c.REF}\" PLINK2=\"{c.PLINK2}\" THREADS=\"{c.THREADS}\" RESULTS=\"{c.RESULTS}\" RESULTS2=\"{c.RESULTS2}\" AUDIT=\"{c.AUDIT}\"')")"
