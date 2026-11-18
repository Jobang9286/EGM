#!/bin/bash
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
# Wait for native-depth peak calling, then build both heatmaps.
# run_all_pr.sh re-runs the whole chain, but every earlier step skips work that
# already exists, so in practice this is just 04_heatmap.
set -x
L=$ROOT/4FP_blood/_setup/peaks_native.log
until grep -q PEAKS_NATIVE_DONE "$L" 2>/dev/null; do sleep 20; done
R=$ROOT/4FP_blood/_setup/run.sh
$R h3k56ac_noinput  $ROOT/4FP_blood/h3k56ac  4FP 4FP_H3K56ac  CP CP_H3K56ac  > $ROOT/4FP_blood/_setup/hm_k56.log 2>&1
echo "HM_K56_RC=$?"
$R h3k27me3_noinput $ROOT/4FP_blood/h3k27me3 4FP 4FP_H3K27me3 CP CP_H3K27me3 > $ROOT/4FP_blood/_setup/hm_k27.log 2>&1
echo "HM_K27_RC=$?"
echo HEATMAP_DONE at $(date)
