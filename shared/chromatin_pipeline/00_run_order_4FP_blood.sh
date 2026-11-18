#!/bin/bash
ROOT=$PWD

## ============================================================================
## 1. Environment
## ============================================================================
cd $ROOT/4FP_blood/_setup
nohup ./mkenv.sh > mkenv.log 2>&1 &

## the bioconda MACS2 build does not import on RHEL 9
cd $ROOT/4FP_blood/_setup
setsid nohup ./fixmacs2.sh > fixmacs2.log 2>&1 < /dev/null &

## ============================================================================
## 2. Reference
## ============================================================================
cd $ROOT/4FP_blood/_setup/ref
nohup ./getref.sh > getref.log 2>&1 &

## mm10 ENCODE blacklist v2 lifted to mm39
bash $ROOT/4FP_blood/_setup/ref/blacklist_liftover_mm39.sh

## drop the RefSeq-named prebuilt index, build from the GENCODE FASTA
cd $ROOT/4FP_blood/_setup/ref
rm -f GRCm39.*.bt2 GRCm39.zip GRCm39.dict GRCm39.manifest.json mm10-bl3.bed mm10-blacklist.v2.bed
mkdir -p index
nohup ./build_index.sh > build_index.log 2>&1 &

## GENCODE TSS/gene BEDs, MGI GAF, TF list and GO ontology for steps 7-10
## (independent of steps 3-6)
bash $ROOT/4FP_blood/_setup/ref/annotation_refs.sh

## ============================================================================
## 3. FASTQ names
## ============================================================================
## the pipeline expects <sample>_R1/_R2.fastq.gz
cd $ROOT/4FP_blood
for m in h3k56ac:H3K56ac h3k27me3:H3K27me3; do
  d=${m%%:*}; M=${m##*:}
  cd $ROOT/4FP_blood/$d
  for arm in 4FP CP; do
    ln -sf ${arm}_${M}_1.fastq.gz ${arm}_${M}_R1.fastq.gz
    ln -sf ${arm}_${M}_2.fastq.gz ${arm}_${M}_R2.fastq.gz
  done
  echo "=== $d ==="; ls -l | grep _R
done

## read counts (diagnostic only)
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate chipseq
cd $ROOT/4FP_blood
nohup bash -c 'for f in h3k56ac/4FP_H3K56ac_1.fastq.gz h3k56ac/CP_H3K56ac_1.fastq.gz h3k27me3/4FP_H3K27me3_1.fastq.gz h3k27me3/CP_H3K27me3_1.fastq.gz; do
  echo "$f $(( $(unpigz -c -p 8 "$f" | wc -l) / 4 )) pairs"; done' > _setup/readcounts.txt 2>&1 &

## ============================================================================
## 4. H3K27me3 and H3K56ac: align, filter, pseudo-replicates, peaks,
##    region classes and bigWigs
## ============================================================================
## align, filter, pseudo-replicates (both marks)
cd $ROOT/4FP_blood
setsid nohup _setup/align_only.sh > _setup/align.log 2>&1 < /dev/null &

## MACS2 on each pseudo-replicate at native depth
cd $ROOT/4FP_blood/_setup
setsid nohup ./peaks_native.sh > peaks_native.log 2>&1 < /dev/null &

## run.sh -> run_all_pr.sh; finished steps are skipped, so this runs
## 04_heatmap.sh (h3k56ac_noinput, then h3k27me3_noinput: regions, bigWigs, images)
cd $ROOT/4FP_blood
setsid nohup _setup/heatmap.sh > _setup/heatmap.log 2>&1 < /dev/null &

## ============================================================================
## 5. H3K27me3 heatmap: final render
##    heatmap/ = scale-regions, heatmap_center/ = point-centred +/-3 kb;
##    same regions and bigWigs.
## ============================================================================
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate chipseq
cd $ROOT/4FP_blood/chip_heatmap
W=$ROOT/4FP_blood/h3k27me3
pk(){ echo "$W/${1}_peaks.broadPeak.filt"; }; bam(){ echo "$W/alignments/${1}.dedup.bam"; }
for co in "h3k27me3_noinput:heatmap" "h3k27me3_noinput_center:heatmap_center"; do
  cfg=${co%%:*}; out=${co##*:}
  bin/04_heatmap.sh "$cfg" "$W/$out" \
    4FP "$(pk 4FP_H3K27me3_ps1),$(pk 4FP_H3K27me3_ps2)" "$(bam 4FP_H3K27me3_ps1),$(bam 4FP_H3K27me3_ps2)" \
    CP  "$(pk CP_H3K27me3_ps1),$(pk CP_H3K27me3_ps2)"   "$(bam CP_H3K27me3_ps1),$(bam CP_H3K27me3_ps2)" 2>&1 | tail -1
done

## ============================================================================
## 6. H3K27ac, then ATAC
## ============================================================================
## FASTQ names for H3K27ac
cd $ROOT/4FP_blood/h3k27ac
for a in 4FP CP; do
  ln -sf ${a}_H3K27ac_1.fastq.gz ${a}_H3K27ac_R1.fastq.gz
  ln -sf ${a}_H3K27ac_2.fastq.gz ${a}_H3K27ac_R2.fastq.gz
done

## H3K27ac end to end
cd $ROOT/4FP_blood/_setup
setsid nohup ./run_k27ac.sh > run_k27ac.log 2>&1 < /dev/null &

## ATAC end to end (run_atac.sh waits for K27AC_DONE)
cd $ROOT/4FP_blood/atac
for a in 4FP CP; do ln -sf ${a}_1.fastq.gz ${a}_R1.fastq.gz; ln -sf ${a}_2.fastq.gz ${a}_R2.fastq.gz; done
cd $ROOT/4FP_blood/_setup
setsid nohup ./run_atac.sh > run_atac.log 2>&1 < /dev/null &

## ============================================================================
## 7. Nearest-TSS annotation of the region classes
## ============================================================================
## H3K56ac (final) and H3K27me3 (superseded below)
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate chipseq
cd $ROOT/4FP_blood/chip_heatmap
echo "############ H3K56ac ############"
bin/05_annotate.sh h3k56ac_noinput $ROOT/4FP_blood/h3k56ac/heatmap 2>/dev/null
echo "############ H3K27me3 (center regions = same as scale-regions) ############"
bin/05_annotate.sh h3k27me3_noinput $ROOT/4FP_blood/h3k27me3/heatmap 2>/dev/null

## ATAC, H3K27ac and H3K27me3 with the final 05_annotate.sh
## (gene-less contigs -> "unassigned"); these are the annotation files used downstream.
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate chipseq
cd $ROOT/4FP_blood/chip_heatmap
for p in "atac:atac" "h3k27ac_noinput:h3k27ac" "h3k27me3_noinput:h3k27me3"; do
  bin/05_annotate.sh "${p%%:*}" $ROOT/4FP_blood/${p##*:}/heatmap 2>/dev/null | grep -E 'annotating|unassigned'
done

## ============================================================================
## 8. TF and GO over-representation (exploratory, console output only)
## ============================================================================
## TF over-representation
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate chipseq
cd $ROOT/4FP_blood
chmod +x chip_heatmap/bin/06_tf_enrich.py
python chip_heatmap/bin/06_tf_enrich.py _setup/ref/mouse_TF_genes.txt \
  atac/heatmap h3k27ac/heatmap h3k56ac/heatmap h3k27me3/heatmap 2>&1 | tail -24

## GO BP against a same-mark background
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate chipseq
cd $ROOT/4FP_blood
python chip_heatmap/bin/07_go_enrich.py atac/heatmap h3k27ac/heatmap --top 10 2>&1 | tail -50
python chip_heatmap/bin/07_go_enrich.py h3k27me3/heatmap --p 0.05 --top 12 2>&1 | tail -32

## 4FP-only vs CP-only genes directly
python chip_heatmap/bin/08_go_contrast.py h3k27ac/heatmap atac/heatmap --top 10 2>&1 | grep -v Warning
python chip_heatmap/bin/08_go_contrast.py h3k56ac/heatmap h3k27me3/heatmap --top 8 2>&1 | grep -v Warning

## ============================================================================
## 9. Tables for R
##     export_R/ is copied to the local project folder chromatin/ (regions_annotated,
##     genes, go_bp, summary, genes_de_long, gene_matrix .tsv and by_mark/)
## ============================================================================
## regions, gene annotation and GO BP tables
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate chipseq
cd $ROOT/4FP_blood
rm -rf export_R/by_mark
python chip_heatmap/bin/09_export.py $ROOT/4FP_blood/export_R \
  ATAC:atac/heatmap H3K27ac:h3k27ac/heatmap H3K27me3:h3k27me3/heatmap 2>&1 | grep -v Warning

## join to the RNA-seq DE table (genes_de_long.tsv is read by Figure1/)
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate chipseq
cd $ROOT/4FP_blood
python chip_heatmap/bin/10_join_de.py export_R $ROOT/4FP_export/4FP_CP_DE_limma_voom.csv \
  ATAC:atac/heatmap H3K27ac:h3k27ac/heatmap H3K27me3:h3k27me3/heatmap

## ============================================================================
## 10. Next: Figure1/Fig1I-J_promoter_signal_ranking.sh
##     builds _setup/ref/pc_promoters_2kb.bed, runs _setup/signal.sh and collapses the
##     promoter signal per gene for the fgsea ranking.
## ============================================================================
