#!/bin/bash
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
set -ex
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate chipseq
cd $ROOT/4FP_blood/_setup/ref/index
ln -sf $ROOT/mouse_ref/build/genome.fa GRCm39.fa
samtools faidx GRCm39.fa
bowtie2-build --threads 64 GRCm39.fa GRCm39
bowtie2-inspect -n GRCm39 | head -25
echo BUILD_DONE
