#!/bin/bash
ROOT=$PWD
T=${ROOT}/39_fibroblasts/2_RNA_seq_tissue

## ---- brain (run in EGM/brain) ------------------------------------------------
cd ${T}/EGM/brain
STAR --runThreadN 32 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./Con_Brain__1_1.fastq.gz ./Con_Brain__1_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix brain_control_rep1 --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 16 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./Con_Brain__2_1.fastq.gz ./Con_Brain__2_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix brain_control_rep2 --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 32 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./Con_Brain__3_1.fastq.gz ./Con_Brain__3_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix brain_control_rep3 --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 16 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./EGM_Brain__1_1.fastq.gz ./EGM_Brain__1_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix brain_egm_rep1 --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 16 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./EGM_Brain__2_1.fastq.gz ./EGM_Brain__2_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix brain_egm_rep2 --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 16 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./EGM_Brain__3_1.fastq.gz ./EGM_Brain__3_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix brain_egm_rep3 --outSAMtype BAM SortedByCoordinate

## ---- kidney (run in EGM/kidney) ----------------------------------------------
cd ${T}/EGM/kidney
STAR --runThreadN 32 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./Con_Kidney-1_1.fastq.gz ./Con_Kidney-1_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix kidney_control_rep1 --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 32 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./Con_Kidney-2_1.fastq.gz ./Con_Kidney-2_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix kidney_control_rep2 --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 32 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./Con_Kidney-3_1.fastq.gz ./Con_Kidney-3_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix kidney_control_rep3 --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 16 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./EGM_Kidney-1_1.fastq.gz ./EGM_Kidney-1_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix kidney_egm_rep1 --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 16 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./EGM_Kidney-2_1.fastq.gz ./EGM_Kidney-2_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix kidney_egm_rep2 --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 16 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./EGM_Kidney-3_1.fastq.gz ./EGM_Kidney-3_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix kidney_egm_rep3 --outSAMtype BAM SortedByCoordinate

## ---- liver (run in EGM/liver) ------------------------------------------------
cd ${T}/EGM/liver
STAR --runThreadN 32 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./Con_Liver-1_1.fastq.gz ./Con_Liver-1_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix liver_control_rep1 --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 16 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./Con_Liver-2_1.fastq.gz ./Con_Liver-2_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix liver_control_rep2 --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 32 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./Con_Liver-3_1.fastq.gz ./Con_Liver-3_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix liver_control_rep3 --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 32 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./EGM_Liver-1_1.fastq.gz ./EGM_Liver-1_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix liver_egm_rep1 --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 16 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./EGM_Liver-2_1.fastq.gz ./EGM_Liver-2_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix liver_egm_rep2 --outSAMtype BAM SortedByCoordinate
STAR --runThreadN 32 --genomeDir ${ROOT}/STAR_ref/GRCm39 --readFilesIn ./EGM_Liver-3_1.fastq.gz ./EGM_Liver-3_2.fastq.gz --readFilesCommand zcat --outFileNamePrefix liver_egm_rep3 --outSAMtype BAM SortedByCoordinate
