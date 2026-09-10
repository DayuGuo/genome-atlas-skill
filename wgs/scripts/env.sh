# source this before any shell step: reads config.yaml through wgsconfig.py so shell and Python agree.
# usage:  source scripts/env.sh
_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
eval "$(python3 - "$_here" <<'PY'
import sys, pathlib, shlex
sys.path.insert(0, sys.argv[1])
import wgsconfig as c
out = {
  "PROJ": str(c.P), "WGS": str(c.W), "REF_DIR": str(c.REF), "TOOLS": str(c.TOOLS), "PGSDIR": str(c.PGS),
  "SAMPLE": c.SAMPLE, "SEX": c.SEX, "THREADS": str(c.THREADS), "MEM_GB": str(c.MEM_GB),
  "CRAM": c.READS, "YBAM": c.Y_READS, "VENDOR_VCF": c.VENDOR_VCF,
  "REF": c.FASTA, "REF_PATH": c.FASTA, "REF38": c.FASTA38,
  "CHAIN": c.CHAIN_19_38, "CHAIN_BACK": c.CHAIN_38_19,
  "KG_PFILE": c.KG_PFILE, "KG_VCF_DIR": c.KG_VCF_DIR, "CLINVAR": c.CLINVAR, "GFF3": c.GFF3,
  "PRS_REF": c.PRS_REF, "AADR": c.AADR, "YTREE": c.YTREE,
  "PLINK2": c.PLINK2, "PICARD": c.PICARD, "PHARMCAT_DIR": c.PHARMCAT_DIR, "FLARE": c.FLARE,
  "JAVA": c.JAVA, "WHATSHAP": c.WHATSHAP,
  "MIN_DP": str(c.MIN_DP), "MIN_MQ": str(c.MIN_MQ),
  "LA_A": ",".join(c.LA_NORTH), "LA_B": ",".join(c.LA_SOUTH), "LA_CTRL": ",".join(c.LA_CONTROL),
  "LA_LAB_A": c.LA_LABELS[0], "LA_LAB_B": c.LA_LABELS[1],
  "SUPERPOP": c.SUPERPOP, "SUBPOPS": ",".join(c.SUBPOPS),
}
for k, v in out.items():
    print(f"export {k}={shlex.quote(v)}")
PY
)"
unset _here
# tools installed by setup/00_install_tools.sh live in a conda env; put it first if it exists
[ -d "$TOOLS/env/bin" ] && export PATH="$TOOLS/env/bin:$PATH"
