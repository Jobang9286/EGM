#!/bin/bash
# Wait for the STAR index build, verify STAR<->RSEM consistency, then run the target pipeline.
export PATH=$HOME/.conda/envs/rnaseq/bin:$PATH
REF=mouse_ref
L=fib_target/logs/chain.log
ts(){ date '+%F %T'; }
echo "[$(ts)] waiting for STAR index build" > "$L"
while pgrep -f "STAR-avx2 --runThreadN 48 --runMode genomeGenerate" > /dev/null; do sleep 20; done
echo "[$(ts)] STAR process exited" >> "$L"

if ! grep -q "finished successfully" "$REF/logs/star_build.log"; then
  echo "[$(ts)] ABORT: STAR did not report success" >> "$L"; tail -5 "$REF/logs/star_build.log" >> "$L"; exit 1
fi
for f in Genome SA SAindex chrName.txt transcriptInfo.tab exonGeTrInfo.tab; do
  [ -s "$REF/star_index/$f" ] || { echo "[$(ts)] ABORT: missing $f" >> "$L"; exit 1; }
done
echo "[$(ts)] STAR index OK ($(du -sh $REF/star_index | cut -f1))" >> "$L"

# STAR and RSEM must agree on the transcript set, or rsem-calculate-expression
# fails on the transcriptome BAM.
ST=$(head -1 "$REF/star_index/transcriptInfo.tab" | awk '{print $1}')
RS=$(grep -c '^>' "$REF/rsem/GRCm39.transcripts.fa")
echo "[$(ts)] transcripts  STAR=$ST  RSEM=$RS" >> "$L"
[ "$ST" = "$RS" ] || { echo "[$(ts)] ABORT: transcript count mismatch" >> "$L"; exit 1; }

echo "[$(ts)] launching target pipeline" >> "$L"
bash fib_target/run_target.sh >> "$L" 2>&1
echo "[$(ts)] chain done" >> "$L"
