#!/usr/bin/env bash
# plink2 has no ARM64 (aarch64) binary; build from source against a local OpenBLAS (conda-forge) prefix.
# x86_64 users: download a binary from https://www.cog-genomics.org/plink/2.0/ into $TOOLS/plink2 instead.
set -euo pipefail
source "$(dirname "$0")/../scripts/env.sh"; mkdir -p "$TOOLS"; cd "$TOOLS"
conda create -y -q -p ./blas -c conda-forge openblas liblapack
[ -d plink-ng ] || git clone -q --depth 1 https://github.com/chrchang/plink-ng.git
cd plink-ng/2.0 && B="$TOOLS/blas"
make -j"$(nproc)" CXXFLAGS="-O2 -std=c++14 -I$B/include -DNDEBUG" CFLAGS="-O2 -I$B/include -DNDEBUG" BLASFLAGS64="-L$B/lib -lopenblas -llapack -Wl,-rpath,$B/lib"
cp bin/plink2 "$TOOLS/plink2" && "$TOOLS/plink2" --version
