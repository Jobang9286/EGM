#!/bin/bash
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
# Per-promoter CPM signal for each arm, from the arm-mean bigwigs.
# A continuous score avoids the threshold amplification the presence/absence
# peak classes introduce, and it is the closer analogue of a per-gene
# correlation coefficient used as GSEA input.
set -x
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate chipseq
R=$ROOT/4FP_blood/_setup/ref
B=$ROOT/4FP_blood
for m in h3k27me3 h3k27ac atac; do
  multiBigwigSummary BED-file \
    -b $B/$m/heatmap/4FP.bw $B/$m/heatmap/CP.bw \
    --BED $R/pc_promoters_2kb.bed \
    --labels 4FP CP \
    -p 32 --outRawCounts $B/rank/${m}_promoter_signal.tsv \
    -o $B/rank/${m}_promoter.npz
  echo "DONE_$m"
done
echo SIGNAL_DONE
