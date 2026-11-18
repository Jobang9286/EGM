#!/bin/bash
GTF=$PWD/STAR_ref/mm39.gtf

## 1) only_egm_effect_in_tissue/counts_egm.txt -> counts_egm.csv,
##    counts_egm_kid_liv.csv, counts_egm_{kid,liv,brain}.csv
featureCounts -T 64 -a ${GTF} --ignoreDup -o counts_egm.txt Brain_HGPS_control_1.bam Brain_HGPS_control_2.bam Brain_HGPS_control_3.bam Brain_HGPS_EGM_1.bam Brain_HGPS_EGM_2.bam Brain_HGPS_EGM_3.bam Kidney_HGPS_control_1.bam Kidney_HGPS_control_2.bam Kidney_HGPS_control_3.bam Kidney_HGPS_EGM_1.bam Kidney_HGPS_EGM_2.bam Kidney_HGPS_EGM_3.bam Liver_HGPS_control_1.bam Liver_HGPS_control_2.bam Liver_HGPS_control_3.bam Liver_HGPS_EGM_1.bam Liver_HGPS_EGM_2.bam Liver_HGPS_EGM_3.bam

## 2) Liver_final/counts_liver_for_time.txt -> Liver_final/counts_liver_for_time.csv,
##    counts_liver_revised.csv, counts_liver_only_egm.csv
featureCounts -T 64 -a ${GTF} --ignoreDup -o counts.txt Liver_HGPS_control_1.bam Liver_HGPS_control_2.bam Liver_HGPS_control_3.bam Liver_HGPS_EGM_1.bam Liver_HGPS_EGM_2.bam Liver_HGPS_EGM_3.bam Liver_WT_12m_1.bam Liver_WT_12m_2.bam Liver_WT_12m_3.bam Liver_WT_12m_4.bam Liver_WT_18m_1.bam Liver_WT_18m_2.bam Liver_WT_18m_3.bam Liver_WT_1m_1.bam Liver_WT_1m_2.bam Liver_WT_1m_3.bam Liver_WT_1m_4.bam Liver_WT_1m_5.bam Liver_WT_1m_6.bam Liver_WT_21M_1.bam Liver_WT_21M_2.bam Liver_WT_21M_3.bam Liver_WT_21M_4.bam Liver_WT_21M_5.bam Liver_WT_21M_6.bam Liver_WT_24M_1.bam Liver_WT_24M_2.bam Liver_WT_24M_3.bam Liver_WT_27m_1.bam Liver_WT_27m_2.bam Liver_WT_27m_3.bam Liver_WT_27m_4.bam Liver_WT_3M_1.bam Liver_WT_3M_2.bam Liver_WT_3M_3.bam Liver_WT_3M_4.bam Liver_WT_3M_5.bam Liver_WT_3M_6.bam Liver_WT_6m_1.bam Liver_WT_6m_2.bam Liver_WT_6m_3.bam Liver_WT_6m_4.bam

## 3) Brain_final/bam/counts_brain.txt -> Brain_final/bam/counts_brain_time.csv,
##    counts_brain_revised.csv, counts_brain_only_egm.csv
featureCounts -T 64 -a ${GTF} --ignoreDup -o counts.txt Brain_HGPS_control_1.bam Brain_HGPS_control_2.bam Brain_HGPS_control_3.bam Brain_HGPS_EGM_1.bam Brain_HGPS_EGM_2.bam Brain_HGPS_EGM_3.bam Brain_WT_12M_1.bam Brain_WT_12M_2.bam Brain_WT_12M_3.bam Brain_WT_12M_4.bam Brain_WT_18M_1.bam Brain_WT_18M_2.bam Brain_WT_18M_3.bam Brain_WT_18M_4.bam Brain_WT_1M_1.bam Brain_WT_1M_2.bam Brain_WT_1M_3.bam Brain_WT_1M_4.bam Brain_WT_1M_5.bam Brain_WT_1M_6.bam Brain_WT_21M_1.bam Brain_WT_21M_2.bam Brain_WT_21M_3.bam Brain_WT_21M_4.bam Brain_WT_21M_5.bam Brain_WT_21M_6.bam Brain_WT_24M_1.bam Brain_WT_24M_2.bam Brain_WT_24M_3.bam Brain_WT_24M_4.bam Brain_WT_27M_1.bam Brain_WT_27M_2.bam Brain_WT_27M_3.bam Brain_WT_27M_4.bam Brain_WT_3M_1.bam Brain_WT_3M_2.bam Brain_WT_3M_3.bam Brain_WT_3M_4.bam Brain_WT_3M_5.bam Brain_WT_3M_6.bam Brain_WT_6M_1.bam Brain_WT_6M_2.bam Brain_WT_6M_3.bam Brain_WT_6M_4.bam

## 4) Kidney_final/total/counts_kidney_final.txt -> Kidney_final/total/counts_kidney_*.csv,
##    incl. counts_kidney_time.csv and counts_kidney_only_egm.csv
featureCounts -T 64 -a ${GTF} --ignoreDup -o counts.txt Kidney_HGPS_control_1.bam Kidney_HGPS_control_2.bam Kidney_HGPS_control_3.bam Kidney_HGPS_EGM_1.bam Kidney_HGPS_EGM_2.bam Kidney_HGPS_EGM_3.bam Kidney_WT_1M_1.bam Kidney_WT_1M_2.bam Kidney_WT_1M_3.bam Kidney_WT_1M_4.bam Kidney_WT_21M_1.bam Kidney_WT_21M_2.bam Kidney_WT_21M_3.bam Kidney_WT_21M_4.bam Kidney_WT_21M_5.bam Kidney_WT_21M_6.bam Kidney_WT_24M_1.bam Kidney_WT_24M_2.bam Kidney_WT_24M_3.bam Kidney_WT_27M_1.bam Kidney_WT_27M_2.bam Kidney_WT_27M_3.bam Kidney_WT_27M_4.bam Kidney_WT_3M_1.bam Kidney_WT_3M_2.bam Kidney_WT_3M_3.bam Kidney_WT_6M_1.bam Kidney_WT_6M_2.bam Kidney_WT_6M_3.bam Kidney_WT_6M_4.bam

## 5) public/kidney/bam/counts_kidney.txt -> public/kidney/bam/counts_kidney.csv
##    (earlier kidney set; superseded by 4 and 6)
featureCounts -T 64 -a ${GTF} --ignoreDup -o counts.txt Kidney_HGPS_control_1.bam Kidney_HGPS_control_2.bam Kidney_HGPS_control_3.bam Kidney_HGPS_EGM_1.bam Kidney_HGPS_EGM_2.bam Kidney_HGPS_EGM_3.bam Kidney_WT_1M_1.bam Kidney_WT_1M_2.bam Kidney_WT_1M_3.bam Kidney_WT_1M_4.bam Kidney_WT_1M_5.bam Kidney_WT_1M_6.bam Kidney_WT_21M_1.bam Kidney_WT_21M_2.bam Kidney_WT_21M_3.bam Kidney_WT_21M_4.bam Kidney_WT_21M_5.bam Kidney_WT_21M_6.bam Kidney_WT_24M_1.bam Kidney_WT_24M_2.bam Kidney_WT_24M_3.bam Kidney_WT_27M_1.bam Kidney_WT_27M_2.bam Kidney_WT_27M_3.bam Kidney_WT_27M_4.bam Kidney_WT_3M_1.bam Kidney_WT_3M_2.bam Kidney_WT_3M_3.bam Kidney_WT_3M_5.bam Kidney_WT_3M_6.bam Kidney_WT_6M_1.bam Kidney_WT_6M_2.bam Kidney_WT_6M_3.bam Kidney_WT_6M_4.bam

## 6) public/kidney/bam/bam/counts_kidney_time.txt -> public/kidney/bam/bam/counts_kidney_time.csv
##    (WT kidney aging series 1-27 months)
featureCounts -T 64 -a ${GTF} --ignoreDup -o counts.txt Kidney_WT_12M_1.bam Kidney_WT_12M_2.bam Kidney_WT_12M_3.bam Kidney_WT_12M_4.bam Kidney_WT_18M_1.bam Kidney_WT_18M_2.bam Kidney_WT_18M_3.bam Kidney_WT_18M_4.bam Kidney_WT_1M_1.bam Kidney_WT_1M_2.bam Kidney_WT_1M_3.bam Kidney_WT_1M_4.bam Kidney_WT_1M_5.bam Kidney_WT_1M_6.bam Kidney_WT_24M_1.bam Kidney_WT_24M_2.bam Kidney_WT_24M_3.bam Kidney_WT_27M_1.bam Kidney_WT_27M_2.bam Kidney_WT_27M_3.bam Kidney_WT_27M_4.bam Kidney_WT_3M_1.bam Kidney_WT_3M_2.bam Kidney_WT_3M_3.bam Kidney_WT_3M_5.bam Kidney_WT_3M_6.bam Kidney_WT_6M_1.bam Kidney_WT_6M_2.bam Kidney_WT_6M_3.bam Kidney_WT_6M_4.bam
