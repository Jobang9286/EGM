#!/bin/bash
export PATH=$HOME/.conda/envs/rnaseq/bin:$PATH
REF=$PWD/mouse_ref
cd "$REF"; mkdir -p build logs rsem star_index

## ---- download: GENCODE M38 GRCm39 primary-assembly genome FASTA and GENCODE M38 GTF
GENCODE_FA_URL="<GENCODE M38 GRCm39 primary assembly genome .fa.gz URL>"
GENCODE_GTF_URL="<GENCODE M38 annotation .gtf.gz URL>"
bash pget.sh "$GENCODE_GTF_URL" genome.gtf.gz 16 >> logs/ref_dl.log 2>&1
bash pget.sh "$GENCODE_FA_URL"  genome.fa.gz  16 >> logs/ref_dl.log 2>&1
echo REF_DL_DONE >> logs/ref_dl.log

## ---- decompress
zcat genome.fa.gz  > build/genome.fa
zcat genome.gtf.gz > build/genome.gtf
echo DECOMP_DONE >> logs/decomp.log

## ---- transcript -> gene-symbol map for RSEM: one line per transcript,
## "<gene_name><TAB><transcript_id with version>", in GTF order
awk -F'\t' '$3=="transcript"{match($9,/gene_name "[^"]+"/); g=substr($9,RSTART+11,RLENGTH-12);
  match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART+15,RLENGTH-16); print g"\t"t}' \
  build/genome.gtf > build/t2g_symbol.txt

## ---- RSEM reference
rsem-prepare-reference --gtf ${REF}/build/genome.gtf \
  --transcript-to-gene-map ${REF}/build/t2g_symbol.txt \
  ${REF}/build/genome.fa ${REF}/rsem/GRCm39 > logs/rsem_build.log 2>&1 && echo RSEM_REF_DONE >> logs/rsem_build.log

## ---- STAR index
STAR-avx2 --runThreadN 48 --runMode genomeGenerate --genomeDir ${REF}/star_index \
  --genomeFastaFiles ${REF}/build/genome.fa --sjdbGTFfile ${REF}/build/genome.gtf \
  --sjdbOverhang 100 --limitGenomeGenerateRAM 128000000000 \
  --outFileNamePrefix ${REF}/logs/star_ > logs/star_build.log 2>&1 && echo STAR_INDEX_DONE >> logs/star_build.log
