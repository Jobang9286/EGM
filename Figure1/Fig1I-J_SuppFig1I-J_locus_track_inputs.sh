#!/bin/bash
ROOT=$PWD

## ---- Cdkn2a / HoxA, 200-bp bins, three marks --------------------------------
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate chipseq
B=$ROOT/4FP_blood; R=$B/_setup/ref
mkdir -p $B/tracks && cd $B/tracks
# two loci, generous flanks so the domain edges are visible
printf "chr4\t89170000\t89250000\tCdkn2a_locus\nchr6\t52120000\t52250000\tHoxA_cluster\n" > loci.bed
bedtools makewindows -b loci.bed -w 200 -i src | sort -k1,1 -k2,2n > bins.bed
wc -l bins.bed
for m in atac h3k27ac h3k27me3; do
  multiBigwigSummary BED-file -b $B/$m/heatmap/4FP.bw $B/$m/heatmap/CP.bw \
    --BED bins.bed --labels 4FP CP -p 16 \
    --outRawCounts ${m}_bins.tsv -o ${m}_bins.npz 2>/dev/null
done
# stitch into one tidy table, re-attaching the locus name
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{k[$1":"$2"-"$3]=$4; next}
  FNR==1{next}{key=$1":"$2"-"$3; print MARK, k[key], $1, $2, $3, $4, $5}' MARK=x bins.bed atac_bins.tsv > /dev/null
{ echo -e "mark\tlocus\tchrom\tstart\tend\tsig_4FP\tsig_CP"
  for m in ATAC:atac H3K27ac:h3k27ac H3K27me3:h3k27me3; do
    lab=${m%%:*}; f=${m##*:}
    awk -F'\t' -v M="$lab" 'BEGIN{OFS="\t"} NR==FNR{k[$1":"$2"-"$3]=$4; next}
      FNR==1{next}{key=$1":"$2"-"$3; if(key in k) print M, k[key], $1, $2, $3, $4, $5}' bins.bed ${f}_bins.tsv
  done
} > locus_signal.tsv
wc -l locus_signal.tsv; head -3 locus_signal.tsv
awk -F'\t' 'NR>1{n[$1"|"$2]++} END{for(k in n) print k, n[k]}' locus_signal.tsv

## ---- Ccnb1 (+/- 15 kb) and Itga5, ATAC, 100-bp bins -------------------------
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate chipseq
B=$ROOT/4FP_blood; R=$B/_setup/ref; cd $B/tracks
grep -P '\tCcnb1\t' $R/gencode_vM38_genes.bed
printf "chr13\t100420000\t100460000\tCcnb1\nchr15\t103240000\t103290000\tItga5\n" > loci3.bed
# use the real Ccnb1 coordinates
awk -F'\t' '$4=="Ccnb1"{s=$2-15000; e=$3+15000; print $1"\t"s"\t"e"\tCcnb1"}' $R/gencode_vM38_genes.bed > loci3.bed
printf "chr15\t103240000\t103290000\tItga5\n" >> loci3.bed
cat loci3.bed
bedtools makewindows -b loci3.bed -w 100 -i src | sort -k1,1 -k2,2n > bins3.bed
multiBigwigSummary BED-file -b $B/atac/heatmap/4FP.bw $B/atac/heatmap/CP.bw \
  --BED bins3.bed --labels 4FP CP -p 16 --outRawCounts atac_bins3.tsv -o atac_bins3.npz 2>/dev/null
{ echo -e "locus\tchrom\tstart\tend\tsig_4FP\tsig_CP"
  awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{k[$1":"$2"-"$3]=$4; next}
    FNR==1{next}{key=$1":"$2"-"$3; if(key in k) print k[key], $1, $2, $3, $4, $5}' bins3.bed atac_bins3.tsv
} > atac_locus_signal2.tsv
wc -l atac_locus_signal2.tsv
echo "--- protein-coding genes in the Ccnb1 view ---"
read c s e n < <(head -1 loci3.bed)
awk -F'\t' -v c=$c -v s=$s -v e=$e '$1==c && $3>s && $2<e && $5=="protein_coding"' $R/gencode_vM38_genes.bed

## ---- local, after copying tracks/locus_signal.tsv from the HPC server -------
## H3K27me3 rows only, for the reproduce_H3K27me3 package
C=$ROOT/chromatin; RP=$C/reproduce_H3K27me3
awk -F'\t' 'NR==1 || $1=="H3K27me3"' "$C/tracks/locus_signal.tsv" > "$RP/data/locus_signal_H3K27me3.tsv"
