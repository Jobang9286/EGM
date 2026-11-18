#!/bin/bash
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
# Alignment + filtering for all four libraries. Independent of peak calling, so
# it can run while the MACS2 rebuild is pending.
set -x
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate chipseq
B=$ROOT/4FP_blood/chip_heatmap/bin
run_mark() {
  local preset=$1 wd=$2; shift 2
  for s in "$@"; do
    $B/01_align.sh  "$preset" "$wd" "$s" || return 1
    $B/02_filter.sh "$preset" "$wd" "$s" || return 1
    $B/02b_pseudorep.sh "$preset" "$wd" "$s" || return 1
  done
}
run_mark h3k56ac_noinput  $ROOT/4FP_blood/h3k56ac  4FP_H3K56ac  CP_H3K56ac  > $ROOT/4FP_blood/_setup/align_k56.log 2>&1 &
P1=$!
run_mark h3k27me3_noinput $ROOT/4FP_blood/h3k27me3 4FP_H3K27me3 CP_H3K27me3 > $ROOT/4FP_blood/_setup/align_k27.log 2>&1 &
P2=$!
wait $P1; echo "ALIGN_K56_RC=$?"
wait $P2; echo "ALIGN_K27_RC=$?"
echo ALIGN_DONE at $(date)
