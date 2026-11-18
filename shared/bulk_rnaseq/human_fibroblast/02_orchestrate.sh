#!/bin/bash
# Chain the remaining RNA_mcp phases once trimming and the STAR index build are done.
BASE=4_human_egm
PY=$HOME/RNA_mcp/.venv/bin/python
LOG=$BASE/_logs
STARLOG=human_ref/logs/star_build.log
stamp() { echo "[$(date '+%F %T')] $*"; }
stamp "orchestrator start"
while pgrep -f 'run_rna_mcp.py step1' >/dev/null; do sleep 20; done
stamp "step1 driver finished"
until grep -q 'finished successfully' "$STARLOG"; do
  pgrep -f 'runMode genomeGenerate' >/dev/null || { stamp "STAR index build exited without success"; exit 1; }
  sleep 20
done
[ -s human_ref/RNA_hg38/SAindex ] && [ -s human_ref/RNA_hg38/GRCh38.grp ] || { stamp "reference files missing"; exit 1; }
stamp "reference ready; launching align phase (run_full_pipeline)"
$PY -I $BASE/_scripts/run_rna_mcp.py align > $LOG/driver_align.log 2>&1
stamp "align phase exited rc=$?"
stamp "launching de + qc phases"
$PY -I $BASE/_scripts/run_rna_mcp.py de qc > $LOG/driver_de_qc.log 2>&1
stamp "de+qc phases exited rc=$?"
stamp "ORCHESTRATION_DONE"
