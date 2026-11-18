#!/bin/bash
cd public_monkey
tail -n +2 download_list.csv | while IFS=, read organ age crr title u1 u2; do
  d="fastq/$organ"; mkdir -p "$d"
  for u in "$u1" "$u2"; do
    f="$d/$(basename $u)"
    echo "$f"
  done
done > allfiles.txt
tail -n +2 download_list.csv | awk -F, '{print $3","$1","$5"\n"$3","$1","$6}' | \
  xargs -P4 -I{} bash -c '
    IFS=, read crr organ url <<< "{}"
    d="fastq/$organ"; mkdir -p "$d"
    out="$d/$(basename $url)"
    curl -sS -C - --retry 5 --retry-delay 10 "$url" -o "$out" && echo "OK $out" || echo "FAIL $out"
  ' >> dl_progress.log 2>&1
echo "ALL_DONE $(date)" >> dl_progress.log
