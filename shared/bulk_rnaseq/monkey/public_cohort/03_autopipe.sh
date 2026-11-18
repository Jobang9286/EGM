#!/bin/bash
# Autonomous server-side aligner for GSA CRA012195 monkey RNA-seq.
# Stages completed download pairs, then aligns pending samples per organ using the
# SAME commands as rna_mcp (trim_galore -> STAR -> RSEM). Runs one sample at a time
# (RAM-safe). Survives SSH disconnect (launched via nohup setsid).
export PATH=$HOME/.conda/envs/rnaseq/bin:$PATH
ROOT=$PWD/public_monkey
REF=$ROOT/reference
GENOME=macFas
LOG=$ROOT/autopipe.log
STAT=$ROOT/autopipe.status
STAR_T=96
RSEM_T=64
ts(){ date '+%F %T'; }

align_one(){
  local pd=$1 s=$2
  cd "$pd" || return 1
  mkdir -p trim_galore fastqc
  [ -f ${s}_1.fastq.gz ] || ln -sf ${s}_1.fq.gz ${s}_1.fastq.gz
  [ -f ${s}_2.fastq.gz ] || ln -sf ${s}_2.fq.gz ${s}_2.fastq.gz
  if [ ! -f trim_galore/${s}_1_val_1.fq.gz ] || [ ! -f trim_galore/${s}_2_val_2.fq.gz ]; then
    echo "[$(ts)] $s step1 trim_galore" >>"$LOG"
    trim_galore --paired --gzip -o trim_galore/ ${s}_1.fastq.gz ${s}_2.fastq.gz >>"$LOG" 2>&1 || { echo "[$(ts)] $s step1 FAIL" >>"$LOG"; return 1; }
  fi
  if [ ! -f ${s}Log.final.out ] || [ ! -f ${s}Aligned.toTranscriptome.out.bam ]; then
    echo "[$(ts)] $s step2 STAR" >>"$LOG"
    rm -rf TEMP_${s}; rm -f ${s}Aligned.out.sam
    STAR --runThreadN $STAR_T --genomeDir "$REF" \
      --readFilesIn trim_galore/${s}_1_val_1.fq.gz trim_galore/${s}_2_val_2.fq.gz \
      --readFilesCommand zcat --outTmpDir TEMP_${s} \
      --outFileNamePrefix ${s} --quantMode TranscriptomeSAM >>"$LOG" 2>&1 || { echo "[$(ts)] $s step2 FAIL" >>"$LOG"; return 1; }
    rm -f ${s}Aligned.out.sam
  fi
  if [ ! -f ${s}.genes.results ]; then
    echo "[$(ts)] $s step3 RSEM" >>"$LOG"
    rm -rf ${s}.temp ${s}.stat
    rsem-calculate-expression --paired-end --bam -p $RSEM_T \
      ${s}Aligned.toTranscriptome.out.bam "$REF/$GENOME" ${s} >>"$LOG" 2>&1 || { echo "[$(ts)] $s step3 FAIL" >>"$LOG"; return 1; }
  fi
  echo "[$(ts)] $s DONE" >>"$LOG"
}

echo "[$(ts)] autopipe started (pid $$)" >>"$LOG"
while true; do
  (cd "$ROOT/.." && bash "$ROOT/stage.sh") >>"$LOG" 2>&1
  for o in Brain Liver Kidney Spleen; do
    pd="$ROOT/proj/$o"; [ -d "$pd" ] || continue
    exec 9>"$pd/.rna_pipeline_run.lock"
    if flock -n 9; then
      for f1 in "$pd"/*_1.fq.gz; do
        [ -e "$f1" ] || continue
        s=$(basename "$f1" _1.fq.gz)
        [ -f "$pd/${s}.genes.results" ] && continue
        align_one "$pd" "$s"
        cd "$ROOT"
      done
      flock -u 9
    fi
    exec 9>&-
  done
  # status snapshot
  { echo "[$(ts)] dl.sh=$(ps aux|grep -c '[d]l.sh')";
    for o in Brain Liver Kidney Spleen; do
      echo "  $o staged=$(ls "$ROOT/proj/$o"/*_1.fq.gz 2>/dev/null|wc -l) done=$(ls "$ROOT/proj/$o"/*.genes.results 2>/dev/null|wc -l)";
    done; } > "$STAT"
  # stop condition
  dl=$(ps aux|grep -c '[d]l.sh')
  allok=1; anyzero=0
  for o in Brain Liver Kidney Spleen; do
    st=$(ls "$ROOT/proj/$o"/*_1.fq.gz 2>/dev/null|wc -l); dn=$(ls "$ROOT/proj/$o"/*.genes.results 2>/dev/null|wc -l)
    [ "$st" -eq 0 ] && anyzero=1
    [ "$st" -ne "$dn" ] && allok=0
  done
  if [ "$dl" -le 1 ] && [ "$allok" -eq 1 ] && [ "$anyzero" -eq 0 ]; then
    echo "[$(ts)] ALL COMPLETE — autopipe exiting" >>"$LOG"; break
  fi
  sleep 180
done
