#!/bin/bash
# 05_annotate.sh PRESET HEATMAP_DIR [REGION_BED ...]
#   Annotate region sets against GENCODE vM38 genes.
#   Defaults to the three sets 04_heatmap.sh writes.
#   Produces <region>.annotated.tsv, <region>.genes.txt and a printed summary.
#
# Nearest-TSS assignment, not a causal claim: a peak is reported with the gene
# whose TSS is closest, plus the signed distance, so a distal peak sitting
# 300 kb from anything is visible as such rather than silently credited to a
# gene. Classes are promoter (|d| <= 2 kb), genic (inside the gene body),
# proximal (|d| <= 50 kb) and distal.
#
# GENCODE vM38 is the same annotation the RNA-seq used, so gene symbols join
# straight onto the DE table with no ortholog or build mapping.
set -euo pipefail
source "$(dirname "$0")/_common.sh"
hm_load_config "$1"; HM="$2"; shift 2
hm_need bedtools sort
TSS="$REF_DIR/gencode_vM38_tss.bed"
GENES="$REF_DIR/gencode_vM38_genes.bed"
for f in "$TSS" "$GENES"; do
  [ -s "$f" ] || { echo "missing annotation reference: $f" >&2; exit 1; }
done
PROMOTER=2000
PROXIMAL=50000

REGIONS=("$@")
if [ "${#REGIONS[@]}" -eq 0 ]; then
  for f in "$HM"/*_specific.bed "$HM"/shared.bed; do [ -s "$f" ] && REGIONS+=("$f"); done
fi

for BED in "${REGIONS[@]}"; do
  base="${BED%.bed}"
  hm_log "annotating $(basename "$BED") ($(wc -l < "$BED") regions)"
  sort -k1,1 -k2,2n "$BED" | cut -f1-3 > "$base.sorted.tmp"

  # nearest TSS, ties broken to the first, signed distance (-D b: upstream of a
  # + strand gene is negative)
  # bedtools warns once per file about mixed chr/scaffold naming; it is not an
  # error here, both files use the same GENCODE names.
  bedtools closest -a "$base.sorted.tmp" -b "$TSS" -D b -t first \
    2> >(grep -v 'inconsistent naming convention' >&2) > "$base.closest.tmp"
  # does the peak sit inside any gene body
  bedtools intersect -a "$base.sorted.tmp" -b "$GENES" -c \
    2> >(grep -v 'inconsistent naming convention' >&2) > "$base.genic.tmp"

  { printf "chrom\tstart\tend\tgene\tgene_type\tstrand\tdist_to_tss\tclass\n"
    paste "$base.closest.tmp" <(cut -f4 "$base.genic.tmp") \
    | awk -F'\t' -v OFS='\t' -v P="$PROMOTER" -v X="$PROXIMAL" '{
        d=$10; ad=(d<0?-d:d); ingene=$11
        # bedtools closest reports "." and distance -1 when a contig carries
        # no gene at all; left alone that reads as |d|=1, i.e. a promoter.
        if ($7 == ".") { print $1,$2,$3,".","NA",".","NA","unassigned"; next }
        cls = (ad<=P) ? "promoter" : (ingene>0 ? "genic" : (ad<=X ? "proximal" : "distal"))
        print $1,$2,$3,$7,$8,$9,d,cls
      }'
  } > "$base.annotated.tsv"

  awk -F'\t' 'NR>1 && $4!="." {print $4}' "$base.annotated.tsv" | sort -u > "$base.genes.txt"
  awk -F'\t' 'NR>1 && $4!="." && $5=="protein_coding" {print $4}' "$base.annotated.tsv" | sort -u > "$base.genes_pc.txt"

  echo "  class distribution:"
  awk -F'\t' 'NR>1{c[$8]++} END{for(k in c) printf "    %-10s %6d  (%.1f%%)\n", k, c[k], 100*c[k]/(NR-1)}' \
    "$base.annotated.tsv" | sort -k2 -rn
  echo "  genes: $(wc -l < "$base.genes.txt") total, $(wc -l < "$base.genes_pc.txt") protein-coding"
  rm -f "$base.sorted.tmp" "$base.closest.tmp" "$base.genic.tmp"
done
