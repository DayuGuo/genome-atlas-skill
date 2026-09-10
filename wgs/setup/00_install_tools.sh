#!/usr/bin/env bash
# Install every command-line tool the pipeline needs into a self-contained conda environment under work/tools/env,
# plus the Java and Python tools that are not on conda. Works on x86_64 and aarch64.
set -euo pipefail
source "$(dirname "$0")/../scripts/env.sh"
MAMBA=${MAMBA:-$(command -v mamba || command -v micromamba || command -v conda)}
[ -n "$MAMBA" ] || { echo "need mamba/micromamba/conda on PATH"; exit 1; }

# 1. core toolchain
"$MAMBA" create -y -p "$TOOLS/env" -c conda-forge -c bioconda \
  "samtools>=1.19" "bcftools>=1.19" htslib mosdepth delly bedtools expansionhunter t1k \
  "openjdk>=17" pysam cyvcf2 pandas numpy scipy pyarrow requests pyyaml pyliftover matplotlib pip
# WhatsHap pins an older Python, so it gets its own environment
"$MAMBA" create -y -p "$TOOLS/whatshap" -c conda-forge -c bioconda whatshap "python=3.10"

# 2. plink2
mkdir -p "$TOOLS"
if [ ! -x "$TOOLS/plink2" ]; then
  ARCH=$(uname -m)
  if [ "$ARCH" = "x86_64" ]; then
    curl -sSL -o /tmp/p2.zip https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241222.zip
    unzip -o -q /tmp/p2.zip -d "$TOOLS" && chmod +x "$TOOLS/plink2"
  else
    echo "plink2 has no aarch64 build: run setup/00b_build_plink2_arm64.sh"
  fi
fi

# 3. Java and Python tools that are not packaged
mkdir -p "$TOOLS/pharmcat" "$TOOLS/flare"
PC=${PHARMCAT_VERSION:-3.4.0}
B=https://github.com/PharmGKB/PharmCAT/releases/download/v$PC
for f in pharmcat-$PC-all.jar pharmcat-preprocessor-$PC.tar.gz pharmcat_positions_$PC.vcf.bgz pharmcat_positions_$PC.vcf.bgz.csi; do
  [ -f "$TOOLS/pharmcat/$f" ] || curl -sSL -o "$TOOLS/pharmcat/$f" "$B/$f"
done
tar xzf "$TOOLS/pharmcat/pharmcat-preprocessor-$PC.tar.gz" -C "$TOOLS/pharmcat"
[ -f "$TOOLS/picard.jar" ] || curl -sSL -o "$TOOLS/picard.jar" \
  "$(curl -s https://api.github.com/repos/broadinstitute/picard/releases/latest | grep -o 'https[^"]*picard.jar')"
[ -f "$TOOLS/flare/flare.jar" ] || curl -sSL -o "$TOOLS/flare/flare.jar" https://faculty.washington.edu/browning/flare.jar
[ -d "$TOOLS/Cyrius" ] || git clone -q https://github.com/Illumina/Cyrius.git "$TOOLS/Cyrius"
[ -d "$TOOLS/SMNCopyNumberCaller" ] || git clone -q https://github.com/Illumina/SMNCopyNumberCaller.git "$TOOLS/SMNCopyNumberCaller"
"$TOOLS/env/bin/pip" install -q -r "$TOOLS/Cyrius/requirements.txt" -r "$TOOLS/pharmcat/preprocessor/requirements.txt"
[ -x "$TOOLS/haplogrep3" ] || { curl -sSL -o /tmp/h3.zip https://github.com/genepi/haplogrep3/releases/latest/download/haplogrep3-3.2.2-linux.zip; unzip -o -q /tmp/h3.zip -d "$TOOLS"; chmod +x "$TOOLS/haplogrep3"; }
echo "tools installed under $TOOLS"
