#!/bin/bash
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
# H3K27ac end to end, same chain as the other two marks.
set -x
$ROOT/4FP_blood/_setup/run.sh h3k27ac_noinput \
  $ROOT/4FP_blood/h3k27ac 4FP 4FP_H3K27ac CP CP_H3K27ac
echo "K27AC_RC=$?"
echo K27AC_DONE at $(date)
