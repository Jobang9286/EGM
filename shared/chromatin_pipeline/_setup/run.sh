#!/bin/bash
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
# run.sh PRESET WORKDIR LABEL_A SAMPLE_A LABEL_B SAMPLE_B
# Activates the chipseq env, then hands off to run_all_pr.sh.
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate chipseq
exec $ROOT/4FP_blood/chip_heatmap/bin/run_all_pr.sh "$@"
