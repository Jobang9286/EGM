#!/bin/bash
STAR --runThreadN 32 --genomeDir STAR_ref/mm10 --readFilesIn 22_JY_RNA_seq/raw_data/GME-1_1.fastq.gz 22_JY_RNA_seq/raw_data/GME-1_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix GME_1
STAR --runThreadN 32 --genomeDir STAR_ref/mm10 --readFilesIn 22_JY_RNA_seq/raw_data/GME-2_1.fastq.gz 22_JY_RNA_seq/raw_data/GME-2_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix GME_2
STAR --runThreadN 32 --genomeDir STAR_ref/mm10 --readFilesIn 22_JY_RNA_seq/raw_data/GME-3_1.fastq.gz 22_JY_RNA_seq/raw_data/GME-3_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix GME_3
STAR --runThreadN 32 --genomeDir STAR_ref/mm10 --readFilesIn 22_JY_RNA_seq/raw_data/GFP-1_1.fastq.gz 22_JY_RNA_seq/raw_data/GFP-1_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix negative_GME_1
STAR --runThreadN 32 --genomeDir STAR_ref/mm10 --readFilesIn 22_JY_RNA_seq/raw_data/GFP-2_1.fastq.gz 22_JY_RNA_seq/raw_data/GFP-2_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix negative_GME_2
STAR --runThreadN 32 --genomeDir STAR_ref/mm10 --readFilesIn 22_JY_RNA_seq/raw_data/GFP-3_1.fastq.gz 22_JY_RNA_seq/raw_data/GFP-3_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix negative_GME_3
