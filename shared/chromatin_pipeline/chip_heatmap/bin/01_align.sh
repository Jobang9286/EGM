#!/bin/bash
# 01_align.sh PRESET WORKDIR SAMPLE
#   bowtie2 emits no @RG; Picard 3.5 MarkDuplicates then dies on a null read
#   group, so one is stamped per sample here.
#   Expects WORKDIR/<SAMPLE>_R1.fastq.gz (+ _R2 when PAIRED=1).
#   Produces WORKDIR/alignments/<SAMPLE>.sorted.bam
set -euo pipefail
source "$(dirname "$0")/_common.sh"
hm_load_config "$1"; WD="$2"; S="$3"
hm_need bowtie2 samtools fastqc
ALN="$WD/alignments"; mkdir -p "$ALN" "$WD/trim_galore" "$WD/fastqc"
cd "$WD"

if [ -s "$ALN/$S.sorted.bam" ]; then hm_log "$S: alignment exists, skipping"; exit 0; fi

if [ "$PAIRED" -eq 1 ]; then
  [ -s "trim_galore/${S}_R1_val_1.fq.gz" ] || \
    "$TRIM_GALORE" ${TRIM_GALORE_EXTRA:-} --paired --gzip -o trim_galore/ "${S}_R1.fastq.gz" "${S}_R2.fastq.gz"
  fastqc -t "$THREADS" --nogroup -o fastqc/ "trim_galore/${S}_R1_val_1.fq.gz" || true
  hm_log "$S: bowtie2 (paired)"
  bowtie2 -p "$THREADS" -x "$REF_DIR/$GENOME" \
    --rg-id "$S" --rg "SM:$S" --rg "LB:$S" --rg "PL:ILLUMINA" \
    -1 "trim_galore/${S}_R1_val_1.fq.gz" -2 "trim_galore/${S}_R2_val_2.fq.gz" \
    -S "$ALN/$S.sam" 2> "$ALN/$S.bowtie2Log"
else
  [ -s "trim_galore/${S}_trimmed.fq.gz" ] || \
    "$TRIM_GALORE" ${TRIM_GALORE_EXTRA:-} --gzip -o trim_galore/ "${S}.fastq.gz"
  fastqc -t "$THREADS" --nogroup -o fastqc/ "trim_galore/${S}_trimmed.fq.gz" || true
  hm_log "$S: bowtie2 (single-end)"
  bowtie2 -p "$THREADS" -x "$REF_DIR/$GENOME" \
    --rg-id "$S" --rg "SM:$S" --rg "LB:$S" --rg "PL:ILLUMINA" \
    -U "trim_galore/${S}_trimmed.fq.gz" \
    -S "$ALN/$S.sam" 2> "$ALN/$S.bowtie2Log"
fi

samtools view -bS "$ALN/$S.sam" > "$ALN/$S.bam"
samtools sort -@ "$THREADS" "$ALN/$S.bam" > "$ALN/$S.sorted.bam"
samtools index "$ALN/$S.sorted.bam"
rm -f "$ALN/$S.sam" "$ALN/$S.bam"
samtools quickcheck "$ALN/$S.sorted.bam"
hm_log "$S: done -> $ALN/$S.sorted.bam"
