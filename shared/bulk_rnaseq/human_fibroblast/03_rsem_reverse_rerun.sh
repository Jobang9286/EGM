#!/bin/bash
ROOT=$PWD

set -u
BASE=$ROOT/4_human_egm; PROJ=$BASE/project; PY=$HOME/RNA_mcp/.venv/bin/python
# 1) wrapper: real rsem-calculate-expression with --strandedness reverse
cat > $BASE/_bin/rsem-calculate-expression-reverse <<EOF
#!/bin/bash
exec $HOME/.conda/envs/rnaseq/bin/rsem-calculate-expression --strandedness reverse "\$@"
EOF
chmod +x $BASE/_bin/rsem-calculate-expression-reverse
# 2) preserve the unstranded run
mkdir -p $PROJ/unstranded_run
cd $PROJ
mv *.genes.results *.isoforms.results *.transcript.bam unstranded_run/ 2>&1
for s in Control-1 Control-2 Control-3 EGM-1 EGM-2 EGM-3; do mv $s.stat unstranded_run/; done
mv txi.rsem.counts.txt txi.rsem.abundance.TPM.txt DESeq2_results.csv DESeq2_results_annotated.csv DESeq2_normalized_counts.tsv MA.png Volcano.png Heatmap.png unstranded_run/
cp .rna_checkpoint.json unstranded_run/rna_checkpoint.unstranded.json
# 3) drop step3 from the checkpoint so step3_rsem re-runs
$PY -I - <<EOF
import json
p = "$ROOT/4_human_egm/project/.rna_checkpoint.json"
d = json.load(open(p))
for s, steps in d["completed_steps"].items():
    if s != "_global":
        d["completed_steps"][s] = [x for x in steps if x != "step3"]
d["completed_steps"]["_global"] = [x for x in d["completed_steps"].get("_global", []) if x not in ("step4", "step5")]
json.dump(d, open(p, "w"), indent=2)
print("checkpoint:", {k: v for k, v in d["completed_steps"].items()})
EOF
# 4) phase_rsem_reverse in 01_run_rna_mcp.py passes this wrapper as rsem_bin
# 5) launch: RSEM(reverse) x6 -> step4-8 -> QC
( cd $ROOT && setsid nohup $PY -I $BASE/_scripts/run_rna_mcp.py rsem_reverse de qc > $BASE/_logs/driver_reverse.log 2>&1 < /dev/null & )
