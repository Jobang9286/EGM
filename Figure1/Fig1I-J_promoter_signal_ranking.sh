#!/bin/bash
ROOT=$PWD

## ---- protein-coding promoters, TSS +/- 2 kb ---------------------------------
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate chipseq
R=$ROOT/4FP_blood/_setup/ref
cut -f1,2 $R/index/GRCm39.fa.fai > $R/genome_sizes.txt
awk -F'\t' 'BEGIN{OFS="\t"} $5=="protein_coding"{print $1,$2,$3,$4,0,$6}' $R/gencode_vM38_tss.bed \
 | sort -k1,1 -k2,2n \
 | bedtools slop -b 2000 -g $R/genome_sizes.txt -i - > $R/pc_promoters_2kb.bed
wc -l $R/pc_promoters_2kb.bed
echo "=== 경계 초과 검증 ==="
awk 'NR==FNR{L[$1]=$2; next} $3>L[$1] || $2<0 {n++} END{print (n?n:0)" invalid intervals"}' $R/genome_sizes.txt $R/pc_promoters_2kb.bed
head -2 $R/pc_promoters_2kb.bed

## ---- promoter signal for both arms (runs _setup/signal.sh) ------------------
## Started in the background and polled after 55 s; it can be
## run in the foreground instead.
cd $ROOT/4FP_blood
setsid nohup _setup/signal.sh > _setup/signal.log 2>&1 < /dev/null &
sleep 55; grep -E 'DONE_|SIGNAL_DONE|Error|error' _setup/signal.log | tail -5
ls -la rank/*.tsv 2>/dev/null | awk '{print $5, $9}'

## ---- GO BP gene sets (MGI GAF propagated up the ontology) -------------------
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate chipseq
cd $ROOT/4FP_blood
cat > /tmp/mkgmt.py <<PY
import os, sys
from collections import defaultdict
from importlib.machinery import SourceFileLoader
_e = SourceFileLoader("e","$ROOT/4FP_blood/chip_heatmap/bin/07_go_enrich.py").load_module()
name, ns, parents = _e.parse_obo(_e.OBO); sys.setrecursionlimit(10000)
g2g = _e.load_annotations(_e.GAF, "P", parents)
t2g = defaultdict(set)
for g, ts in g2g.items():
    for t in ts:
        if ns.get(t) == "P":
            t2g[t].add(g)
out = "$ROOT/4FP_blood/rank/go_bp_mouse.gmt"
n = 0
with open(out, "w") as fh:
    for t, gs in sorted(t2g.items()):
        if 10 <= len(gs) <= 500:
            fh.write(f"{t}\t{name.get(t,t)}\t" + "\t".join(sorted(gs)) + "\n")
            n += 1
print(f"{n} GO BP sets (size 10-500) -> {out}")
PY
python /tmp/mkgmt.py
head -1 rank/go_bp_mouse.gmt | cut -c1-140
echo; echo "=== signal file head ==="; head -3 rank/h3k27me3_promoter_signal.tsv

## ---- collapse promoters to genes, log2(P-BP/C-BP) --------------------------
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate chipseq
cd $ROOT/4FP_blood/rank
R=$ROOT/4FP_blood/_setup/ref
for m in h3k27me3 h3k27ac atac; do
  awk -F'\t' 'BEGIN{OFS="\t"}
    NR==FNR { key=$1":"$2"-"$3; g[key]=$4; next }
    FNR==1  { next }
    { key=$1":"$2"-"$3
      if (key in g) { s4[g[key]]+=$4; sc[g[key]]+=$5; n[g[key]]++ } else { miss++ } }
    END {
      print "gene","sig_4FP","sig_CP","log2FC" > "'"$m"'_promoter_gene.tsv"
      for (k in n) {
        a=s4[k]/n[k]; b=sc[k]/n[k]
        printf "%s\t%.6f\t%.6f\t%.6f\n", k, a, b, log((a+0.01)/(b+0.01))/log(2) >> "'"$m"'_promoter_gene.tsv"
      }
      printf "%-10s genes %d   unmatched rows %d\n", "'"$m"'", length(n), (miss?miss:0)
    }' $R/pc_promoters_2kb.bed ${m}_promoter_signal.tsv
done
echo; head -3 h3k27me3_promoter_gene.tsv
echo "=== log2FC 분포 (H3K27me3) ==="
awk -F'\t' 'NR>1{v[NR]=$4; s+=$4; n++} END{
  asort(v); printf "n=%d  mean=%.3f  median=%.3f  p5=%.3f  p95=%.3f\n", n, s/n, v[int(n/2)], v[int(n*.05)], v[int(n*.95)]}' h3k27me3_promoter_gene.tsv
