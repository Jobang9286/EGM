#!/bin/bash
GENOME_DIR="STAR_genome"

STAR --runThreadN 20 --outSAMtype BAM SortedByCoordinate --genomeDir ${GENOME_DIR} --quantMode GeneCounts --readFilesCommand zcat --readFilesIn RNA/4FP-1_1.fastq.gz RNA/4FP-1_2.fastq.gz --outFileNamePrefix 4FP-1_
STAR --runThreadN 10 --outSAMtype BAM SortedByCoordinate --genomeDir ${GENOME_DIR} --quantMode GeneCounts --readFilesCommand zcat --readFilesIn RNA/4FP-2_1.fastq.gz RNA/4FP-2_2.fastq.gz --outFileNamePrefix 4FP-2_
STAR --runThreadN 10 --outSAMtype BAM SortedByCoordinate --genomeDir ${GENOME_DIR} --quantMode GeneCounts --readFilesCommand zcat --readFilesIn RNA/4FP-3_1.fastq.gz RNA/4FP-3_2.fastq.gz --outFileNamePrefix 4FP-3_
STAR --runThreadN 10 --outSAMtype BAM SortedByCoordinate --genomeDir ${GENOME_DIR} --quantMode GeneCounts --readFilesCommand zcat --readFilesIn RNA/CP-1_1.fastq.gz RNA/CP-1_2.fastq.gz --outFileNamePrefix CP-1_
STAR --runThreadN 10 --outSAMtype BAM SortedByCoordinate --genomeDir ${GENOME_DIR} --quantMode GeneCounts --readFilesCommand zcat --readFilesIn RNA/CP-2_1.fastq.gz RNA/CP-2_2.fastq.gz --outFileNamePrefix CP-2_
STAR --runThreadN 10 --outSAMtype BAM SortedByCoordinate --genomeDir ${GENOME_DIR} --quantMode GeneCounts --readFilesCommand zcat --readFilesIn RNA/CP-3_1.fastq.gz RNA/CP-3_2.fastq.gz --outFileNamePrefix CP-3_
