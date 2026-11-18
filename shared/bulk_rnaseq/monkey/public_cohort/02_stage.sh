#!/bin/bash
# Hardlink COMPLETED _f1/_r2 pairs into proj/<organ>/ as _1/_2 (skips files still downloading)
cd public_monkey || exit 1
active=$(ps aux | grep curl | grep -v grep | grep -oE "CRR[0-9]+" | sort -u)
staged=0
for organ in Brain Liver Kidney Spleen; do
  d="fastq/$organ"; [ -d "$d" ] || continue
  mkdir -p "proj/$organ"
  for f1 in "$d"/*_f1.fq.gz; do
    [ -e "$f1" ] || continue
    crr=$(basename "$f1" _f1.fq.gz)
    r2="$d/${crr}_r2.fq.gz"
    [ -e "$r2" ] || continue
    echo "$active" | grep -q "$crr" && continue   # still downloading
    if [ ! -e "proj/$organ/${crr}_1.fq.gz" ]; then
      ln -f "$f1" "proj/$organ/${crr}_1.fq.gz"
      ln -f "$r2" "proj/$organ/${crr}_2.fq.gz"
      staged=$((staged+1)); echo "staged $organ/$crr"
    fi
  done
done
echo "newly staged: $staged"
