#!/bin/bash
# 04_heatmap.sh PRESET OUTDIR LABEL_A "pkA1,pkA2" "bamA1,bamA2" \
#                             LABEL_B "pkB1,pkB2" "bamB1,bamB2"
#
# Three region classes (A-specific / Shared / B-specific) built from within-arm
# replicate-reproducible peaks, then group-mean CPM coverage over each class.
set -euo pipefail
source "$(dirname "$0")/_common.sh"
hm_load_config "$1"
OUT="$2"; LA="$3"; PKA="$4"; BAMA="$5"; LB="$6"; PKB="$7"; BAMB="$8"
hm_need bedtools computeMatrix
# The vendored deeptools 3.3.2 tree stopped importing under matplotlib >= 3.9
# (deeptools/cm.py calls matplotlib.cm.register_cmap), so a site config may
# point PLOTHEATMAP_CMD at its own launcher. See lib/deeptools356_port/.
DTP="$HM_ROOT/lib/deeptools_patched"
PHEAT="${PLOTHEATMAP_CMD:-env PYTHONPATH=$DTP /usr/bin/python3 $DTP/plotheatmap_patched.py}"
mkdir -p "$OUT"; cd "$OUT"

# Fonts are set through a local matplotlibrc because neither installed deeptools
# version exposes --fontSize.
cat > mplrc <<RC
font.size: ${HEATMAP_FONT_SIZE:-10}
axes.titlesize: ${HEATMAP_TITLE_SIZE:-11}
xtick.labelsize: ${HEATMAP_XTICK_SIZE:-8}
ytick.labelsize: ${HEATMAP_YTICK_SIZE:-9.5}
legend.fontsize: ${HEATMAP_LEGEND_SIZE:-7.5}
figure.subplot.hspace: 0.12
RC
export MATPLOTLIBRC="$OUT/mplrc"

# ---- regions: reproducible peaks per arm, then the 3-way split ----
pk2bed(){ cut -f1-3 "$1" | sort -k1,1 -k2,2n > "$2"; }
IFS="," read -ra A <<< "$PKA"; IFS="," read -ra B <<< "$PKB"
# Intermediates are keyed by label, never by position: positional names made a
# re-run after an arm-order change silently reuse the other arm cached bigwig.
pk2bed "${A[0]}" "${LA}_r1.bed"; pk2bed "${A[1]}" "${LA}_r2.bed"
pk2bed "${B[0]}" "${LB}_r1.bed"; pk2bed "${B[1]}" "${LB}_r2.bed"
bedtools intersect -a "${LA}_r1.bed" -b "${LA}_r2.bed" > "${LA}_group.bed"
bedtools intersect -a "${LB}_r1.bed" -b "${LB}_r2.bed" > "${LB}_group.bed"
bedtools intersect -a "${LA}_group.bed" -b "${LB}_group.bed" -v > "${LA}_specific.bed"
bedtools intersect -a "${LB}_group.bed" -b "${LA}_group.bed" -v > "${LB}_specific.bed"
bedtools intersect -a "${LA}_group.bed" -b "${LB}_group.bed" > shared.bed
for f in "${LA}_specific.bed" shared.bed "${LB}_specific.bed"; do
  hm_log "$(printf "%-28s %s regions" "$f" "$(wc -l < "$f")")"
  [ -s "$f" ] || { echo "empty region set: $f" >&2; exit 1; }
done

# ---- coverage: CPM bigwig per replicate, then the arm mean ----
IFS="," read -ra BA <<< "$BAMA"; IFS="," read -ra BB <<< "$BAMB"
cpm_bw(){ local bam="$1" out="$2"
  [ -s "$out" ] && return 0
  hm_log "bamCoverage $(basename "$bam")"
  "$CONDA_BIN/bamCoverage" --bam "$bam" --normalizeUsing CPM -o "$out" \
    --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE" --exactScaling \
    --binSize "$BW_BINSIZE" --extendReads "$BW_EXTEND" -p "$THREADS" >/dev/null 2>&1; }
cpm_bw "${BA[0]}" "${LA}_r1.cpm.bw"; cpm_bw "${BA[1]}" "${LA}_r2.cpm.bw"
cpm_bw "${BB[0]}" "${LB}_r1.cpm.bw"; cpm_bw "${BB[1]}" "${LB}_r2.cpm.bw"
[ -s "$LA.bw" ] || "$CONDA_BIN/bigwigAverage" -b "${LA}_r1.cpm.bw" "${LA}_r2.cpm.bw" \
  -o "$LA.bw" -p "$THREADS" >/dev/null 2>&1
[ -s "$LB.bw" ] || "$CONDA_BIN/bigwigAverage" -b "${LB}_r1.cpm.bw" "${LB}_r2.cpm.bw" \
  -o "$LB.bw" -p "$THREADS" >/dev/null 2>&1

# ---- matrix ----
# A point-centred window suits punctate marks. Broad domains vary in width by
# orders of magnitude, so those are scaled to a common body length instead.
# deeptools labels a scaled body TSS/TES, which is wrong here: the regions are
# peak domains, not genes, and the two labels also collide with the flank tick.
if [ "$MATRIX_MODE" = "scale-regions" ]; then
  MODE_ARGS=(scale-regions --regionBodyLength "$MATRIX_BODY")
  XLAB="domain"
  EDGE_ARGS=(--startLabel "$MATRIX_START_LABEL" --endLabel "$MATRIX_END_LABEL")
else
  MODE_ARGS=(reference-point --referencePoint center)
  XLAB=""
  EDGE_ARGS=()
fi
if [ ! -s matrix_group.gz ]; then
  hm_log "computeMatrix ($MATRIX_MODE)"
  computeMatrix "${MODE_ARGS[@]}" \
    -R "${LA}_specific.bed" shared.bed "${LB}_specific.bed" \
    -S "$LA.bw" "$LB.bw" -b "$MATRIX_BEFORE" -a "$MATRIX_AFTER" \
    --skipZeros --missingDataAsZero -p "$THREADS" -o matrix_group.gz >/dev/null 2>&1
fi

# ---- heatmap ----
hm_log "plotHeatmap"
# Region groups to keep as heatmap rows but leave out of the summary profile.
# Only lib/deeptools356_port understands the flag, so it is skipped when a site is
# still on the vendored launcher.
HIDE_ARG=()
if [ -n "${HEATMAP_HIDE_PROFILE_GROUPS:-}" ] && [ -n "${PLOTHEATMAP_CMD:-}" ]; then
  HIDE_ARG=(--hideProfileGroups "$HEATMAP_HIDE_PROFILE_GROUPS")
fi
# Every deeptools --legendLocation choice sits inside the axes, on top of the
# curves. This moves it out to the right instead.
if [ "${HEATMAP_LEGEND_OUTSIDE:-0}" -eq 1 ] && [ -n "${PLOTHEATMAP_CMD:-}" ]; then
  HIDE_ARG+=(--legendOutside)
fi
# deeptools sizes the summary plot from heatmapwidth, so widening the panels
# also inflates it. Set it directly instead.
if [ -n "${HEATMAP_PROFILE_HEIGHT:-}" ] && [ -n "${PLOTHEATMAP_CMD:-}" ]; then
  HIDE_ARG+=(--profileHeight "$HEATMAP_PROFILE_HEIGHT")
fi
if [ -n "${HEATMAP_PROFILE_SPACER:-}" ] && [ -n "${PLOTHEATMAP_CMD:-}" ]; then
  HIDE_ARG+=(--profileSpacer "$HEATMAP_PROFILE_SPACER")
fi
# deeptools hardcodes the heatmap x tick size, so the matplotlibrc above only
# reaches the profile ticks.
if [ -n "${HEATMAP_XTICK_SIZE:-}" ] && [ -n "${PLOTHEATMAP_CMD:-}" ]; then
  HIDE_ARG+=(--heatmapTickSize "$HEATMAP_XTICK_SIZE")
fi
# The gap between the sample columns is set by constrained_layout, not by the
# GridSpec width ratios.
if [ -n "${HEATMAP_PANEL_WSPACE:-}" ] && [ -n "${PLOTHEATMAP_CMD:-}" ]; then
  HIDE_ARG+=(--panelWspace "$HEATMAP_PANEL_WSPACE")
fi
if [ -n "${HEATMAP_PANEL_WPAD:-}" ] && [ -n "${PLOTHEATMAP_CMD:-}" ]; then
  HIDE_ARG+=(--panelWpad "$HEATMAP_PANEL_WPAD")
fi
$PHEAT -m matrix_group.gz -out peakheatmap.png "${HIDE_ARG[@]}" \
  --regionsLabel "${LA}_specific" Shared "${LB}_specific" --xAxisLabel "$XLAB" \
  "${EDGE_ARGS[@]}" \
  --samplesLabel "$LA" "$LB" \
  --colorList "$HEATMAP_COLOR" "$HEATMAP_COLOR" \
  --heatmapHeight "$HEATMAP_HEIGHT" --heatmapWidth "$HEATMAP_WIDTH" \
  --missingDataColor white --dpi 200 >/dev/null 2>&1
hm_log "OK -> $OUT/peakheatmap.png"
