#!/bin/bash
# 03_peaks.sh PRESET WORKDIR SAMPLE [CONTROL_SAMPLE]
#   MACS2 (narrow or broad per preset) + blacklist removal.
#   Produces WORKDIR/<SAMPLE>_peaks.<narrow|broad>Peak.filt
set -euo pipefail
source "$(dirname "$0")/_common.sh"
hm_load_config "$1"; WD="$2"; S="$3"; CTRL="${4:-}"
hm_need macs2 bedtools
ALN="$WD/alignments"; cd "$WD"
SUF="$(hm_peak_suffix)"
OUTP="$WD/${S}_peaks.${SUF}.filt"

if [ -s "$OUTP" ]; then hm_log "$S: peaks exist, skipping"; exit 0; fi
# A preset that wants a control but was given none is a config error, so report
# it before touching any BAM.
if [ "$USE_CONTROL" -eq 1 ] && [ -z "$CTRL" ]; then
  echo "preset sets USE_CONTROL=1 but no control sample was given" >&2
  echo "usage: 03_peaks.sh PRESET WORKDIR SAMPLE CONTROL_SAMPLE" >&2
  exit 1
fi
samtools quickcheck "$ALN/$S.dedup.bam"

# BED input keeps the MACS2 call identical for paired and single-end data.
[ -s "$ALN/$S.bed" ] || bedtools bamtobed -i "$ALN/$S.dedup.bam" > "$ALN/$S.bed"

CTRL_ARG=()
if [ "$USE_CONTROL" -eq 1 ]; then
  [ -s "$ALN/$CTRL.bed" ] || bedtools bamtobed -i "$ALN/$CTRL.dedup.bam" > "$ALN/$CTRL.bed"
  CTRL_ARG=(-c "$ALN/$CTRL.bed")
  hm_log "$S: using control $CTRL"
fi

MODE_ARG=()
[ "$PEAK_MODE" = "broad" ] && MODE_ARG=(--broad --broad-cutoff "$MACS_BROAD_CUTOFF")

rm -rf "macs2_$S"
hm_log "$S: macs2 $PEAK_MODE"
macs2 callpeak -t "$ALN/$S.bed" "${CTRL_ARG[@]}" -n "$S" -g "$MACS_GENOME" \
  -B -q "$MACS_QVALUE" --outdir "macs2_$S" "${MODE_ARG[@]}" $MACS_EXTRA

bedtools intersect -v -a "macs2_$S/${S}_peaks.${SUF}" -b "$BLACKLIST" > "$OUTP"
hm_log "$S: $(wc -l < "$OUTP") peaks -> $OUTP"
