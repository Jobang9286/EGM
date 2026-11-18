#!/bin/bash
# 02_filter.sh PRESET WORKDIR SAMPLE
#   MAPQ filter, chrM removal, duplicate removal.
#   Produces WORKDIR/alignments/<SAMPLE>.dedup.bam
set -euo pipefail
source "$(dirname "$0")/_common.sh"
hm_load_config "$1"; WD="$2"; S="$3"
hm_need samtools bamtools
ALN="$WD/alignments"; TMP="$WD/picard_tmp"; mkdir -p "$TMP"
cd "$WD"

if [ -s "$ALN/$S.dedup.bam" ]; then hm_log "$S: dedup exists, skipping"; exit 0; fi
samtools quickcheck "$ALN/$S.sorted.bam"
rm -f "$ALN/$S.filtered.bam" "$ALN/$S.noChrM.bam" "$ALN/$S.dedup.sam"

# -isProperPair only exists for paired libraries; forcing it on single-end data
# silently empties the BAM.
if [ "$PAIRED" -eq 1 ]; then
  bamtools filter -in "$ALN/$S.sorted.bam" -out "$ALN/$S.filtered.bam" \
    -mapQuality ">=$MAPQ" -isProperPair true
else
  bamtools filter -in "$ALN/$S.sorted.bam" -out "$ALN/$S.filtered.bam" \
    -mapQuality ">=$MAPQ"
fi

if [ "$REMOVE_CHRM" -eq 1 ]; then
  samtools view -h "$ALN/$S.filtered.bam" | grep -v -E "chrM" | \
    samtools view -b -o "$ALN/$S.noChrM.bam" -
else
  mv "$ALN/$S.filtered.bam" "$ALN/$S.noChrM.bam"
fi

hm_log "$S: MarkDuplicates"
java -Xmx${JAVA_XMX}g -jar "$PICARD_JAR" MarkDuplicates \
  TMP_DIR="$TMP" VALIDATION_STRINGENCY=LENIENT \
  I="$ALN/$S.noChrM.bam" O="$ALN/$S.dedup.sam" \
  M="$ALN/$S.duplicate_metrics" REMOVE_DUPLICATES=true AS=true
java -Xmx${JAVA_XMX}g -jar "$PICARD_JAR" SortSam \
  TMP_DIR="$TMP" VALIDATION_STRINGENCY=LENIENT SO=coordinate \
  I="$ALN/$S.dedup.sam" O="$ALN/$S.dedup.bam" CREATE_INDEX=true
rm -f "$ALN/$S.dedup.sam" "$ALN/$S.filtered.bam" "$ALN/$S.noChrM.bam"
samtools quickcheck "$ALN/$S.dedup.bam"
[ -s "$ALN/$S.dedup.bam.bai" ] || samtools index "$ALN/$S.dedup.bam"
hm_log "$S: done -> $ALN/$S.dedup.bam"
