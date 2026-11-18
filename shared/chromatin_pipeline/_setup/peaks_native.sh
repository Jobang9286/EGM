#!/bin/bash
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
# Peak calling at native (unmatched) depth, for both marks.
set -x
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate chipseq
B=$ROOT/4FP_blood/chip_heatmap/bin
for s in 4FP_H3K56ac_ps1 4FP_H3K56ac_ps2 CP_H3K56ac_ps1 CP_H3K56ac_ps2; do
  $B/03_peaks.sh h3k56ac_noinput $ROOT/4FP_blood/h3k56ac "$s"
done
for s in 4FP_H3K27me3_ps1 4FP_H3K27me3_ps2 CP_H3K27me3_ps1 CP_H3K27me3_ps2; do
  $B/03_peaks.sh h3k27me3_noinput $ROOT/4FP_blood/h3k27me3 "$s"
done
echo PEAKS_NATIVE_DONE
