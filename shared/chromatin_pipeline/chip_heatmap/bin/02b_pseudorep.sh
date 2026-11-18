#!/bin/bash
# 02b_pseudorep.sh PRESET WORKDIR SAMPLE
#   Split a dedup BAM into two disjoint pseudo-replicates.
#   Produces WORKDIR/alignments/<SAMPLE>_ps1.dedup.bam and _ps2.dedup.bam
#
# 04_heatmap.sh defines an arm's regions as the peaks reproducible between two
# replicates. This experiment has one library per arm, so the reproducibility
# filter is built from pseudo-replicates instead (the ENCODE fallback for
# unreplicated ChIP): each half is peak-called on its own and only peaks found
# in both survive. It is a self-consistency filter, not a biological replicate,
# and it cannot tell you anything about biological variability.
#
# The split is keyed on read name, so mates stay together, and it is
# deterministic, so a rerun reproduces the same halves. Illumina names end in
# :x:y cluster coordinates, which are independent of genomic position, so their
# parity is a fair coin with respect to peaks; a char-sum covers other names.
set -euo pipefail
source "$(dirname "$0")/_common.sh"
hm_load_config "$1"; WD="$2"; S="$3"
hm_need samtools
ALN="$WD/alignments"
O1="$ALN/${S}_ps1.dedup.bam"; O2="$ALN/${S}_ps2.dedup.bam"

if [ -s "$O1" ] && [ -s "$O2" ]; then hm_log "$S: pseudo-replicates exist, skipping"; exit 0; fi
samtools quickcheck "$ALN/$S.dedup.bam"
hm_log "$S: splitting into pseudo-replicates"

samtools view -h -@ 8 "$ALN/$S.dedup.bam" \
| awk -v C1="samtools view -b -@ 4 -o $O1.tmp -" \
      -v C2="samtools view -b -@ 4 -o $O2.tmp -" \
      -v REF="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789:._-" '
    /^@/ { print | C1; print | C2; next }
    {
      n = split($1, a, ":")
      if (n >= 3 && a[n] ~ /^[0-9]+$/ && a[n-1] ~ /^[0-9]+$/) {
        h = (a[n] + a[n-1]) % 2
      } else {
        s = 0
        for (i = 1; i <= length($1); i++) s += i * index(REF, substr($1, i, 1))
        h = s % 2
      }
      if (h == 0) print | C1; else print | C2
    }
    END { close(C1); close(C2) }'

for o in "$O1" "$O2"; do
  [ -s "$o.tmp" ] || { echo "pseudo-replicate split produced no output: $o" >&2; exit 1; }
  mv "$o.tmp" "$o"
  samtools index -@ 8 "$o"
  samtools quickcheck "$o"
done
n1=$(samtools view -c -@ 8 "$O1"); n2=$(samtools view -c -@ 8 "$O2")
hm_log "$S: ps1 $n1 reads / ps2 $n2 reads"
