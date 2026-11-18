#!/bin/bash
BASE=4_human_egm
PY=$HOME/RNA_mcp/.venv/bin/python
LOG=$BASE/_logs
stamp() { echo "[$(date '+%F %T')] $*"; }
stamp "orchestrator3 start: waiting for all run_rna_mcp drivers (align + sample:Control-2)"
sleep 30
while pgrep -f 'run_rna_mcp.py' >/dev/null; do sleep 20; done
stamp "drivers finished; safety align pass (skips completed samples)"
$PY -I $BASE/_scripts/run_rna_mcp.py align > $LOG/driver_align2.log 2>&1
stamp "safety align pass exited rc=$?"
$PY -I $BASE/_scripts/run_rna_mcp.py de qc > $LOG/driver_de_qc.log 2>&1
stamp "de+qc phases exited rc=$?"
stamp "ORCHESTRATION_DONE"
