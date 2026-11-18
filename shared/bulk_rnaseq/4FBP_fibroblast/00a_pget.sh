#!/bin/bash
# Parallel ranged GET — beats per-connection throttling.  pget.sh <url> <out> <nchunks>
URL="$1"; OUT="$2"; N="${3:-16}"
LEN=$(curl -sSI --max-time 30 "$URL" | awk -F': ' 'tolower($1)=="content-length"{gsub(/\r/,"",$2); print $2}')
[ -z "$LEN" ] && { echo "no content-length for $URL"; exit 1; }
CH=$(( (LEN + N - 1) / N ))
rm -f "$OUT".part.*
export URL OUT CH LEN
seq 0 $((N-1)) | xargs -P"$N" -I{} bash -c '
  i=$1; s=$(( i * CH )); e=$(( s + CH - 1 )); [ $e -ge $LEN ] && e=$((LEN-1))
  curl -sS --retry 6 --retry-delay 5 --retry-all-errors -r ${s}-${e} \
       -o "${OUT}.part.$(printf %03d $i)" "$URL"' _ {}
cat "$OUT".part.* > "$OUT"
GOT=$(stat -c%s "$OUT")
if [ "$GOT" = "$LEN" ]; then rm -f "$OUT".part.*; echo "OK $OUT $GOT bytes"
else echo "SIZE MISMATCH $OUT got=$GOT want=$LEN"; exit 1; fi
