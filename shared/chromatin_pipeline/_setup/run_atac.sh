#!/bin/bash
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
# ATAC end to end. The stock atac preset already fits: Tn5 shift on, no control.
set -x
until grep -q K27AC_DONE $ROOT/4FP_blood/_setup/run_k27ac.log 2>/dev/null; do sleep 20; done
$ROOT/4FP_blood/_setup/run.sh atac \
  $ROOT/4FP_blood/atac 4FP 4FP CP CP
echo "ATAC_RC=$?"
echo ATAC_DONE at $(date)
