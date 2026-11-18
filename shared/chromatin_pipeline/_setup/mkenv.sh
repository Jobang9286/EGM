#!/bin/bash
set -x
source "$(conda info --base)/etc/profile.d/conda.sh"
conda create -y -n chipseq -c conda-forge -c bioconda \
  python=3.10 bowtie2 samtools bamtools bedtools macs2 deeptools \
  fastqc trim-galore cutadapt picard ucsc-bedgraphtobigwig pigz
