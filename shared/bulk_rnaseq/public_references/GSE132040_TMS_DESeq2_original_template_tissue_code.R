## ---- libraries -------------------------------------------------------------
library(tximport)
library(biomaRt)
library(Rsubread)
library(DESeq2)
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("edgeR")
library(gplots)
library(pheatmap)

## ---- count matrix and sample metadata ------------------------------------------
data(humanGender)
library(SummarizedExperiment)
library(ggplot2)
library(RRHO2)

setwd("41_RNA_seq_test/")

mat <- as.data.frame(read.csv("./matrix_all.csv"))
row.names(mat) <- mat[,"gene"]
mat <- as.data.frame(mat[,-1])

meta <- as.data.frame(read.csv("./meta_info.csv"))
row.names(meta) <- meta[,"Sample.name"]
meta <- as.data.frame(meta[,-1])
meta

dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = meta,
                              design = ~ range)

## ---- kidney: aaa vs ccc ----------------------------------------------------------
dds_sk = dds[, dds$tissue %in% c("Kidney")]
meta_sk = colData(dds_sk)

dds_sk  = DESeq(dds_sk)

res_sk <- results(dds_sk, contrast=c("range","aaa","ccc"))
res_sk

resOrdered_sk <- res_sk[order(res_sk$ pvalue),]
summary(resOrdered_sk)

resSig_sk <- subset(resOrdered_sk, pvalue < 0.05)
summary(resSig_sk)

resSig_sk

write.table(resSig_sk, file="kd_young_vs_old.csv", sep=',')
