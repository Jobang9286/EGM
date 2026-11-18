#!/bin/bash
ROOT=$PWD
# H3K27me3, FASTQ to the three-class peak heatmap. Runs on the HPC server, not locally.
# Everything downstream in this directory starts from the outputs of this step.
#
# Design note: one library per arm, so the reproducibility filter inside each
# arm is a pseudo-replicate split (02b_pseudorep.sh), not biological replicates.
# No input/IgG was available, so the preset sets USE_CONTROL=0 and MACS2 falls
# back to local lambda from the ChIP track itself.
set -euo pipefail
B=$ROOT/4FP_blood
P=$B/chip_heatmap

# 01 align -> 02 filter -> 02b pseudo-replicates -> 03 MACS2 broad -> 04 heatmap
$B/_setup/run.sh h3k27me3_noinput $B/h3k27me3 4FP 4FP_H3K27me3 CP CP_H3K27me3
#   -> $B/h3k27me3/heatmap/{4FP_specific,shared,CP_specific}.bed
#   -> $B/h3k27me3/heatmap/peakheatmap.png          (scale-regions, the preset default)
#   -> $B/h3k27me3/heatmap/{4FP,CP}.bw              arm-mean CPM coverage

# The preset displays a scaled domain body. That assumes broad domains; MACS2
# --broad here returned peaks of median width 406 bp covering 0.57% of the
# genome, so the point-centred variant is the faithful one for this data.
$P/bin/04_heatmap.sh h3k27me3_noinput_center $B/h3k27me3/heatmap_center \
  4FP "$B/h3k27me3/4FP_H3K27me3_ps1_peaks.broadPeak.filt,$B/h3k27me3/4FP_H3K27me3_ps2_peaks.broadPeak.filt" \
      "$B/h3k27me3/alignments/4FP_H3K27me3_ps1.dedup.bam,$B/h3k27me3/alignments/4FP_H3K27me3_ps2.dedup.bam" \
  CP  "$B/h3k27me3/CP_H3K27me3_ps1_peaks.broadPeak.filt,$B/h3k27me3/CP_H3K27me3_ps2_peaks.broadPeak.filt" \
      "$B/h3k27me3/alignments/CP_H3K27me3_ps1.dedup.bam,$B/h3k27me3/alignments/CP_H3K27me3_ps2.dedup.bam"

# --- inputs for the local scripts -------------------------------------------
R=$B/_setup/ref

# (a) promoter signal, the GSEA ranking metric.  A continuous score avoids the
#     threshold amplification the presence/absence peak classes introduce.
multiBigwigSummary BED-file \
  -b $B/h3k27me3/heatmap/4FP.bw $B/h3k27me3/heatmap/CP.bw \
  --BED $R/pc_promoters_2kb.bed --labels 4FP CP -p 24 \
  --outRawCounts $B/rank/h3k27me3_promoter_signal.tsv -o /dev/null
# multiBigwigSummary drops the BED name column, so gene names are joined back by
# coordinate; log2FC uses a 0.01 pseudocount.
awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{k[$1":"$2"-"$3]=$4; next}
  FNR==1{next}{key=$1":"$2"-"$3; if(key in k){g=k[key]; a[g]+=$4; c[g]+=$5; n[g]++}}
  END{print "gene","sig_4FP","sig_CP","log2FC"
      for(x in n){p=a[x]/n[x]; q=c[x]/n[x];
        printf "%s\t%.6f\t%.6f\t%.6f\n", x,p,q,log((p+0.01)/(q+0.01))/log(2)}}' \
  $R/pc_promoters_2kb.bed $B/rank/h3k27me3_promoter_signal.tsv \
  > $B/rank/h3k27me3_promoter_gene.tsv

# (b) binned coverage over the two loci shown as tracks
printf "chr4\t89170000\t89250000\tCdkn2a_locus\nchr6\t52120000\t52250000\tHoxA_cluster\n" > /tmp/loci.bed
bedtools makewindows -b /tmp/loci.bed -w 200 -i src | sort -k1,1 -k2,2n > /tmp/bins.bed
multiBigwigSummary BED-file -b $B/h3k27me3/heatmap/4FP.bw $B/h3k27me3/heatmap/CP.bw \
  --BED /tmp/bins.bed --labels 4FP CP -p 16 --outRawCounts /tmp/k27me3_bins.tsv -o /dev/null
{ echo -e "mark\tlocus\tchrom\tstart\tend\tsig_4FP\tsig_CP"
  awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{k[$1":"$2"-"$3]=$4; next}
    FNR==1{next}{key=$1":"$2"-"$3; if(key in k) print "H3K27me3", k[key], $1,$2,$3,$4,$5}' \
    /tmp/bins.bed /tmp/k27me3_bins.tsv
} > $B/tracks/locus_signal_H3K27me3.tsv

# (c) GO biological process sets, from the MGI GAF propagated up the ontology
#     (see chip_heatmap/bin/07_go_enrich.py for the parser these reuse)
