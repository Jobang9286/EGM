#!/bin/bash
# run_all_pr.sh PRESET WORKDIR LABEL_A SAMPLE_A LABEL_B SAMPLE_B [CTRL_A] [CTRL_B]
#
# FASTQ -> heatmap for arms that have ONE library each. The two "replicates"
# 04_heatmap.sh needs are pseudo-replicates made by 02b_pseudorep.sh. Every step
# is resumable: rerunning skips finished samples.
#
#   run_all_pr.sh h3k56ac_noinput /data/k56 4FP 4FP_H3K56ac CP CP_H3K56ac
#
# The pseudo-replicates also carry the coverage: each half is CPM-normalised on
# its own, so their bigwigAverage is the full-depth CPM profile.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/_common.sh"
PRESET="$1"; WD="$2"; LA="$3"; SA="$4"; LB="$5"; SB="$6"; CA="${7:-}"; CB="${8:-}"
hm_load_config "$PRESET"
SUF="$(hm_peak_suffix)"

ALL=("$SA" "$SB")
[ -n "$CA" ] && ALL+=("$CA")
[ -n "$CB" ] && ALL+=("$CB")
for s in "${ALL[@]}"; do
  "$HERE/01_align.sh"  "$PRESET" "$WD" "$s"
  "$HERE/02_filter.sh" "$PRESET" "$WD" "$s"
done

for s in "$SA" "$SB"; do "$HERE/02b_pseudorep.sh" "$PRESET" "$WD" "$s"; done
for s in "${SA}_ps1" "${SA}_ps2"; do "$HERE/03_peaks.sh" "$PRESET" "$WD" "$s" "$CA"; done
for s in "${SB}_ps1" "${SB}_ps2"; do "$HERE/03_peaks.sh" "$PRESET" "$WD" "$s" "$CB"; done

pk(){ echo "$WD/${1}_peaks.${SUF}.filt"; }
bam(){ echo "$WD/alignments/${1}.dedup.bam"; }
"$HERE/04_heatmap.sh" "$PRESET" "$WD/heatmap" \
  "$LA" "$(pk "${SA}_ps1"),$(pk "${SA}_ps2")" "$(bam "${SA}_ps1"),$(bam "${SA}_ps2")" \
  "$LB" "$(pk "${SB}_ps1"),$(pk "${SB}_ps2")" "$(bam "${SB}_ps1"),$(bam "${SB}_ps2")"
