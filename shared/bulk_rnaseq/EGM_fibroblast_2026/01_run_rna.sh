#!/usr/bin/env bash
# EGM blood RNA-seq: trim_galore -> STAR -> RSEM -> count matrix.
# Commands are the RNA_mcp step1/2/3 commands; this driver adds resource
# partitioning, per-sample checkpoints and logging so the whole set can run
# unattended.
set -uo pipefail

PROJ=$PWD/egm_blood
REF_STAR=$PWD/mouse_ref/star_index
REF_RSEM=$PWD/mouse_ref/rsem/GRCm39
CUT=$HOME/.conda/envs/chipseq/bin/cutadapt   # rnaseq env has no cutadapt
export PATH=$HOME/.conda/envs/rnaseq/bin:$PATH  # also fixes fastqc's perl

LOG=$PROJ/_logs; CK=$PROJ/_ckpt; mkdir -p "$LOG" "$CK"
SAMPLES=(GFP-1 GFP-2 GFP-3 GME-1 GME-2 GME-3 Eno1-1 Gclc-1 Map1s-1)

# 128 cores / 187 GB, machine otherwise idle.
P1_JOBS=9;  P1_CORES=4     # trim_galore --cores 4 costs ~12 cores with pigz -> ~108
P2_JOBS=4;  P2_THREADS=30  # STAR holds ~32 GB each; 4 x 32 = 128 GB < 187 GB
P3_JOBS=9;  P3_THREADS=14  # RSEM is CPU-bound, memory-light

done_p()  { [ -f "$CK/$1.$2.done" ]; }
mark()    { touch "$CK/$1.$2.done"; }
stamp()   { echo "[$(date +%H:%M:%S)] $*"; }

step1() {  # trim_galore + fastqc
  local s=$1; done_p "$s" step1 && { stamp "SKIP step1 $s"; return 0; }
  cd "$PROJ" || return 1
  mkdir -p trim_galore fastqc
  trim_galore --path_to_cutadapt "$CUT" --cores $P1_CORES --paired --gzip \
    -o trim_galore/ "${s}_1.fastq.gz" "${s}_2.fastq.gz" >"$LOG/$s.step1.log" 2>&1 || return 1
  fastqc -t 2 --nogroup -o fastqc/ "trim_galore/${s}_1_val_1.fq.gz" >>"$LOG/$s.step1.log" 2>&1 || return 1
  fastqc -t 2 --nogroup -o fastqc/ "trim_galore/${s}_2_val_2.fq.gz" >>"$LOG/$s.step1.log" 2>&1 || return 1
  mark "$s" step1; stamp "OK   step1 $s"
}

step2() {  # STAR -> TranscriptomeSAM
  local s=$1; done_p "$s" step2 && { stamp "SKIP step2 $s"; return 0; }
  cd "$PROJ" || return 1
  rm -rf "TEMP_$s"
  rm -f "${s}Aligned.out.sam" "${s}Aligned.toTranscriptome.out.bam" \
        "${s}Log.final.out" "${s}Log.out" "${s}Log.progress.out" "${s}SJ.out.tab"
  STAR --runThreadN $P2_THREADS --genomeDir "$REF_STAR" \
    --readFilesIn "trim_galore/${s}_1_val_1.fq.gz" "trim_galore/${s}_2_val_2.fq.gz" \
    --readFilesCommand zcat --outTmpDir "TEMP_$s" \
    --outFileNamePrefix "$s" --quantMode TranscriptomeSAM >"$LOG/$s.step2.log" 2>&1 || return 1
  rm -f "${s}Aligned.out.sam"
  [ -s "${s}Log.final.out" ] || { stamp "FAIL step2 $s (no Log.final.out)"; return 1; }
  samtools quickcheck "${s}Aligned.toTranscriptome.out.bam" || { stamp "FAIL step2 $s (bad BAM)"; return 1; }
  mark "$s" step2
  stamp "OK   step2 $s  unique=$(awk -F'\t' '/Uniquely mapped reads %/{gsub(/ /,"",$2);print $2}' "${s}Log.final.out")"
}

step3() {  # RSEM
  local s=$1; done_p "$s" step3 && { stamp "SKIP step3 $s"; return 0; }
  cd "$PROJ" || return 1
  rm -rf "$s.temp" "$s.stat"
  rm -f "$s.genes.results" "$s.isoforms.results" "$s.transcript.bam"
  rsem-calculate-expression --paired-end --bam -p $P3_THREADS \
    "${s}Aligned.toTranscriptome.out.bam" "$REF_RSEM" "$s" >"$LOG/$s.step3.log" 2>&1 || return 1
  local tot; tot=$(awk 'NR>1{s+=$5}END{printf "%.0f", s}' "$s.genes.results" 2>/dev/null)
  [ "${tot:-0}" -ge 1000 ] 2>/dev/null || { stamp "FAIL step3 $s (near-empty: $tot)"; return 1; }
  mark "$s" step3; stamp "OK   step3 $s  total_expected_count=$tot"
}

# run a function over all samples, at most N at a time
phase() {
  local fn=$1 jobs=$2; shift 2
  local pids=() rc=0
  for s in "${SAMPLES[@]}"; do
    while [ "$(jobs -rp | wc -l)" -ge "$jobs" ]; do sleep 5; done
    $fn "$s" & pids+=($!)
  done
  for p in "${pids[@]}"; do wait "$p" || rc=1; done
  return $rc
}

stamp "=== phase 1: trim_galore + fastqc  (${P1_JOBS} x --cores ${P1_CORES}) ==="
phase step1 $P1_JOBS || { stamp "PHASE 1 had failures"; }
stamp "=== phase 2: STAR  (${P2_JOBS} x ${P2_THREADS} threads) ==="
phase step2 $P2_JOBS || { stamp "PHASE 2 had failures"; }
stamp "=== phase 3: RSEM  (${P3_JOBS} x ${P3_THREADS} threads) ==="
phase step3 $P3_JOBS || { stamp "PHASE 3 had failures"; }

stamp "=== phase 4: count matrix ==="
python3 - <<'PY'
import glob, os, pandas as pd
PROJ = "egm_blood"
SAMPLES = ["GFP-1","GFP-2","GFP-3","GME-1","GME-2","GME-3","Eno1-1","Gclc-1","Map1s-1"]
# FASTQs were delivered as GME-*; the condition is EGM, so outputs are labelled EGM-*.
RENAME  = {"GME-1":"EGM-1","GME-2":"EGM-2","GME-3":"EGM-3"}
cnt, tpm = {}, {}
for s in SAMPLES:
    f = f"{PROJ}/{s}.genes.results"
    if not os.path.exists(f):
        print(f"  [warn] missing {f}"); continue
    d = pd.read_csv(f, sep="\t", index_col=0)
    cnt[s], tpm[s] = d["expected_count"], d["TPM"]
if cnt:
    C, T = pd.DataFrame(cnt), pd.DataFrame(tpm)
    C = C.rename(columns=RENAME); T = T.rename(columns=RENAME)
    C.to_csv(f"{PROJ}/txi.rsem.counts.txt", sep="\t")
    T.to_csv(f"{PROJ}/txi.rsem.abundance.TPM.txt", sep="\t")
    C.round().astype(int).to_csv(f"{PROJ}/counts_symbol.csv")
    print(f"  genes={C.shape[0]}  samples={C.shape[1]}  ({', '.join(C.columns)})")
    print("  library sizes (M):",
          ", ".join(f"{s}={C[s].sum()/1e6:.1f}" for s in C.columns))
PY
stamp "=== DONE ==="
