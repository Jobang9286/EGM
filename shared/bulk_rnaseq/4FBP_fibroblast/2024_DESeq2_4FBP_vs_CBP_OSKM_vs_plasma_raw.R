## ---- read counts per gene ---------------------------------------------------
library(Rsubread)
library(DESeq2)
library("Rsamtools")
library("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("GenomicAlignments")
library("BiocParallel")
library("edgeR")


setwd("22_JY_RNA_seq")



##negative control(GFP)_vs_GME_overexpression_to_tailtip_mFb##




bam.files <- list.files(path = "bam_file_plasma/", pattern = ".bam$", full.names = TRUE)

bam.files

txdb <- makeTxDbFromGFF("../STAR_ref/mm10/mm10.gtf", format="gtf")

ebg <- exonsBy(txdb, by="gene")
ebg

register(MulticoreParam(64))



se <- summarizeOverlaps(features=ebg, 
                        reads=bam.files,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )

se
head( assay(se) )
rowRanges(se)
colnames(se)


sample <- c("GME","GME","GME","negative_GME", "negative_GME", "negative_GME", "P_Young", "P_Young", "P_Young", "P_Young")
sample <- c("hFib","hFib","hFib","iDN", "iDN", "iDN", "iDN")
sample <- c("4FP_P","4FP_P","4FP_P","Neg_P", "Neg_P", "Neg_P")


sample <- as.data.frame(sample)

colData(se) <- DataFrame(sample)

colData(se)

se$sample

countdata <- assay(se)


head(countdata, 3)

dds <- DESeqDataSet(se, design = ~ sample)
dds

coldata <- colData(se)
coldata

## ---- DESeq2 -------------------------------------------------------------------
countdata <- assay(se)

countdata
countdata_2 <- countdata
countdata_2
countdata_2 <- countdata_2[which(rowSums(countdata_2) > 0),]
Counts <- countdata_2

condition <- factor(c("GME","GME","GME","negative_GME","negative_GME","negative_GME"))
coldata <- data.frame(row.names = colnames(Counts), condition)


condition <- factor(c("4FP","4FP","4FP","negative_plasma","negative_plasma","negative_plasma"))
coldata <- data.frame(row.names = colnames(Counts), condition)


coldata

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~sample)
dds <- DESeq(dds)
res <- results(dds)
res <- na.omit(res)
res <- res[res$baseMean > 5,]
res

res <- res[order(-res$stat),]
res

res_raw <- res@listData
write.csv(res, "./fib_vs_og_cata.csv")

## ---- save results ------------------------------------------------------------
write.csv(res, "../39_fibroblasts/RNA_seq/final_4FP_Plasma_vs_neg_Plasma/OSKM_vs_plasma_raw.csv")
