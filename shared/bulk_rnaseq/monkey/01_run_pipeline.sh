#!/bin/bash
set -uo pipefail
export PATH=$HOME/.conda/envs/rnaseq/bin:$PATH
PROJ=$PWD/monkey/project
STAR_IDX=$PWD/monkey/organized/reference/star_index
RSEM_REF=$PWD/monkey/organized/reference/rsem/macFas
T=30
cd "$PROJ"
mkdir -p logs trim_galore fastqc
CP="$PROJ/.rna_checkpoint.json"
LK="$PROJ/.rna_checkpoint.lock"
[ -f "$CP" ] || echo '{"samples":{},"completed_steps":{}}' > "$CP"

mark() {  # sample step
  python3 - "$CP" "$LK" "$1" "$2" <<'PY'
import json,sys,fcntl
cp,lk,s,step=sys.argv[1:5]
with open(lk,"a") as lf:
    fcntl.flock(lf,fcntl.LOCK_EX)
    d=json.load(open(cp))
    d.setdefault("completed_steps",{}).setdefault(s,[])
    if step not in d["completed_steps"][s]:
        d["completed_steps"][s].append(step)
    json.dump(d,open(cp,"w"),indent=2,ensure_ascii=False)
    fcntl.flock(lf,fcntl.LOCK_UN)
PY
}
export -f mark
export PATH PROJ STAR_IDX RSEM_REF T CP LK

run_sample() {
  local s="$1"; local log="$PROJ/logs/${s}.log"
  {
    echo "=== $s START $(date) ==="
    # step1: trim_galore + fastqc
    if [ ! -s trim_galore/${s}_1_val_1.fq.gz ] || [ ! -s trim_galore/${s}_2_val_2.fq.gz ]; then
      trim_galore --paired --gzip -o trim_galore/ ${s}_1.fastq.gz ${s}_2.fastq.gz || { echo "STEP1 FAIL $s"; return 1; }
      fastqc -t $T --nogroup -o fastqc/ trim_galore/${s}_1_val_1.fq.gz || true
      fastqc -t $T --nogroup -o fastqc/ trim_galore/${s}_2_val_2.fq.gz || true
    fi
    mark "$s" step1; echo "[step1 done] $s $(date)"
    # step2: STAR
    if [ ! -s ${s}Aligned.toTranscriptome.out.bam ] || [ ! -s ${s}Log.final.out ]; then
      rm -rf TEMP_${s}; rm -f ${s}Aligned.out.sam ${s}Aligned.toTranscriptome.out.bam ${s}Log.final.out ${s}Log.out ${s}Log.progress.out ${s}SJ.out.tab
      STAR --runThreadN $T --genomeDir "$STAR_IDX" \
        --readFilesIn trim_galore/${s}_1_val_1.fq.gz trim_galore/${s}_2_val_2.fq.gz \
        --readFilesCommand zcat --outTmpDir TEMP_${s} \
        --outFileNamePrefix ${s} --quantMode TranscriptomeSAM || { echo "STEP2 FAIL $s"; return 1; }
      rm -f ${s}Aligned.out.sam
    fi
    samtools quickcheck ${s}Aligned.toTranscriptome.out.bam || { echo "STEP2 BAM BAD $s"; return 1; }
    mark "$s" step2; echo "[step2 done] $s $(date)"
    # step3: RSEM
    if [ ! -s ${s}.genes.results ]; then
      rm -rf ${s}.temp ${s}.stat; rm -f ${s}.genes.results ${s}.isoforms.results ${s}.transcript.bam
      rsem-calculate-expression --paired-end --bam -p $T \
        ${s}Aligned.toTranscriptome.out.bam "$RSEM_REF" ${s} || { echo "STEP3 FAIL $s"; return 1; }
    fi
    mark "$s" step3; echo "[step3 done] $s $(date)"
    echo "=== $s DONE $(date) ==="
  } >> "$log" 2>&1
}
export -f run_sample

SAMPLES="EGM_Brain_1 EGM_Brain_2 EGM_Brain_3 EGM_Kidney_1 EGM_Kidney_2 EGM_Kidney_3 EGM_Liver_1 EGM_Liver_2 EGM_Liver_3 EGM_Spleen_1 EGM_Spleen_2 EGM_Spleen_3 control_Brain_1 control_Brain_2 control_Brain_3 control_Kidney_1 control_Kidney_2 control_Kidney_3 control_Liver_1 control_Liver_2 control_Liver_3 control_Spleen_1 control_Spleen_2 control_Spleen_3"
echo "PIPELINE START $(date)" > "$PROJ/logs/_driver.log"
printf '%s\n' $SAMPLES | xargs -P 4 -I{} bash -c 'run_sample "$@"' _ {} >> "$PROJ/logs/_driver.log" 2>&1
echo "PIPELINE END $(date)" >> "$PROJ/logs/_driver.log"
