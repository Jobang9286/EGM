#!/bin/bash
ROOT=$PWD
export PATH=$HOME/.conda/envs/rnaseq/bin:$PATH
cd ${ROOT}/monkey/organized/reference
mkdir -p build rsem star_index

## ---- decompress -----------------------------------------------------------------
zcat genome.fa.gz  > build/genome.fa
zcat genome.gtf.gz > build/genome.gtf
echo DECOMP_DONE > ${ROOT}/monkey/decomp.log

## ---- RSEM reference -------------------------------------------------------------
rsem-prepare-reference --gtf build/genome.gtf build/genome.fa rsem/macFas \
  > ${ROOT}/monkey/rsem_ref_build.log 2>&1

## ---- STAR index -----------------------------------------------------------------
STAR-avx2 --runThreadN 32 --runMode genomeGenerate --genomeDir star_index \
  --genomeFastaFiles build/genome.fa --sjdbGTFfile build/genome.gtf --sjdbOverhang 149 \
  > ${ROOT}/monkey/star_index_build.log 2>&1
