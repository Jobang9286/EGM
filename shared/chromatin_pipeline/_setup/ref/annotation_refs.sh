#!/bin/bash
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"

## ---- GENCODE vM38 gene bodies and TSSs -----------------------------------------
cd $ROOT/4FP_blood/_setup/ref
awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene"{
  match($9, /gene_name "[^"]+"/); nm=substr($9, RSTART+11, RLENGTH-12);
  match($9, /gene_type "[^"]+"/); ty=substr($9, RSTART+11, RLENGTH-12);
  print $1, $4-1, $5, nm, ty, $7
}' $ROOT/mouse_ref/build/genome.gtf | sort -k1,1 -k2,2n > gencode_vM38_genes.bed
awk -F'\t' 'BEGIN{OFS="\t"} {t = ($6=="+") ? $2 : $3-1; print $1, t, t+1, $4, $5, $6}' gencode_vM38_genes.bed \
  | sort -k1,1 -k2,2n > gencode_vM38_tss.bed
echo "genes: $(wc -l < gencode_vM38_genes.bed)   TSS: $(wc -l < gencode_vM38_tss.bed)"

## ---- MGI GO annotation (GAF) ---------------------------------------------------
cd $ROOT/4FP_blood/_setup/ref
curl -sL -o mgi.gaf.gz "https://current.geneontology.org/annotations/mgi.gaf.gz"

## ---- mouse TF genes from the GAF -----------------------------------------------
cd $ROOT/4FP_blood/_setup/ref
# GO:0003700 DNA-binding transcription factor activity (and its RNA Pol II child)
# add the activator/repressor children, which are TF-specific terms
zcat mgi.gaf.gz | awk -F'\t' '
  $5=="GO:0003700" || $5=="GO:0000981" || $5=="GO:0001228" || $5=="GO:0001227" || $5=="GO:0140110" {print $3}' \
  | sort -u > mouse_TF_genes.txt
echo "TF genes: $(wc -l < mouse_TF_genes.txt)"

## ---- GO ontology for 07_go_enrich.py -------------------------------------------
cd $ROOT/4FP_blood/_setup/ref
curl -sL -o go-basic.obo https://current.geneontology.org/ontology/go-basic.obo
echo "terms: $(grep -c '^\[Term\]' go-basic.obo)"
