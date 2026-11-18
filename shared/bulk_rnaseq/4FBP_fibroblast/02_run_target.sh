#!/bin/bash
# 4FP/CP target quantification: trim_galore -> STAR(TranscriptomeSAM) -> RSEM
# Mirrors the project's existing pipeline, with three fixes over autopipe.sh:
#   (1) trim_galore --cores (autopipe.sh omits it -> single-threaded cutadapt)
#   (2) pigz instead of zcat for STAR --readFilesCommand
#   (3) --forward-prob 0  (library measured as dUTP reverse-stranded, ~99% purity)
set -uo pipefail
export PATH=$HOME/.conda/envs/rnaseq/bin:$PATH
PROJ=$PWD/fib_target
STAR_IDX=$PWD/mouse_ref/star_index
RSEM_REF=$PWD/mouse_ref/rsem/GRCm39
TRIM_CORES=8
STAR_T=40
RSEM_T=32
FP=0                      # forward-prob: 0 = reverse-stranded (dUTP)
cd "$PROJ"; mkdir -p logs trim_galore fastqc
export PATH PROJ STAR_IDX RSEM_REF TRIM_CORES STAR_T RSEM_T FP
ts(){ date '+%F %T'; }

run_sample() {
  local s="$1"; local log="$PROJ/logs/${s}.log"
  {
    echo "=== $s START $(ts) ==="
    # step1 trim_galore + fastqc
    if [ ! -s trim_galore/${s}_1_val_1.fq.gz ] || [ ! -s trim_galore/${s}_2_val_2.fq.gz ]; then
      trim_galore --paired --gzip --cores $TRIM_CORES --poly_g -o trim_galore/ \
        ${s}_1.fastq.gz ${s}_2.fastq.gz || { echo "STEP1 FAIL $s"; return 1; }
    fi
    echo "[step1 done] $s $(ts)"
    fastqc -t 8 --nogroup -o fastqc/ trim_galore/${s}_1_val_1.fq.gz trim_galore/${s}_2_val_2.fq.gz || true
    # step2 STAR
    if [ ! -s ${s}Aligned.toTranscriptome.out.bam ] || [ ! -s ${s}Log.final.out ]; then
      rm -rf TEMP_${s}
      rm -f ${s}Aligned.out.sam ${s}Aligned.toTranscriptome.out.bam ${s}Log.final.out
      STAR --runThreadN $STAR_T --genomeDir "$STAR_IDX" \
        --readFilesIn trim_galore/${s}_1_val_1.fq.gz trim_galore/${s}_2_val_2.fq.gz \
        --readFilesCommand "pigz -dc" --outTmpDir TEMP_${s} \
        --outFileNamePrefix ${s} --quantMode TranscriptomeSAM || { echo "STEP2 FAIL $s"; return 1; }
      rm -f ${s}Aligned.out.sam
    fi
    samtools quickcheck ${s}Aligned.toTranscriptome.out.bam || { echo "STEP2 BAM BAD $s"; return 1; }
    echo "[step2 done] $s $(ts)"
    # step3 RSEM
    if [ ! -s ${s}.genes.results ]; then
      rm -rf ${s}.temp ${s}.stat
      rm -f ${s}.genes.results ${s}.isoforms.results ${s}.transcript.bam
      rsem-calculate-expression --paired-end --bam -p $RSEM_T --forward-prob $FP \
        ${s}Aligned.toTranscriptome.out.bam "$RSEM_REF" ${s} || { echo "STEP3 FAIL $s"; return 1; }
    fi
    # reject near-empty RSEM output (RSEM exits 0 on a truncated BAM)
    tot=$(awk 'NR>1{s+=$5} END{printf "%.0f", s}' ${s}.genes.results)
    [ "${tot:-0}" -lt 1000 ] && { echo "STEP3 EMPTY $s total=$tot"; return 1; }
    echo "[step3 done] $s total_expected_count=$tot $(ts)"
    echo "=== $s DONE $(ts) ==="
  } >> "$log" 2>&1
}
export -f run_sample ts

SAMPLES="4FP-1 4FP-2 4FP-3 CP-1 CP-2 CP-3"
echo "TARGET PIPELINE START $(ts)" > logs/_driver.log
printf '%s\n' $SAMPLES | xargs -P 3 -I{} bash -c 'run_sample "$@"' _ {} >> logs/_driver.log 2>&1
echo "TARGET PIPELINE END $(ts)" >> logs/_driver.log
echo "results: $(ls *.genes.results 2>/dev/null | wc -l) / 6" >> logs/_driver.log
# NOTE on --poly_g: this env's trim_galore is a Rust reimplementation (v2.3.0) whose
# poly-G trimming is AUTO-DETECTED from the data by default. Auto-detect would make the
# behaviour vary per sample. Forced ON here because the target is NovaSeq (2-colour, where
# no-signal bases are called as high-quality G — measured up to 971 polyG reads per 200k in
# R2). The PRJNA316110 reference runs are HiSeq 2000/2500 (4-colour) and must be trimmed
# with --no_poly_g so each platform gets the processing its chemistry requires, recorded
# rather than inferred.
