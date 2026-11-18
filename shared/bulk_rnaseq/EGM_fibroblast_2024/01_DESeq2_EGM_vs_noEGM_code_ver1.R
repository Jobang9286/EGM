# Required packages
library('edgeR')    # v3.12.1
library("DESeq2")   # v1.20.0
library('limma')
library('pheatmap') # v3.26.9
sessionInfo()

library("RNAseqQC")
library("DESeq2")
library("ensembldb")
library("dplyr")
library("ggplot2")
library("purrr")
library("tidyr")
library("tibble")
library("magrittr")

library("GenomicFeatures")


setwd('39_fibroblasts/RNA_seq/final_EGM_vs_no_EGM/')
HTseqpath="."
list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=TRUE)
HTseqfiles <-list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=T)
HTseq <- do.call("cbind", lapply(HTseqfiles, read.csv, header = TRUE, row.names=1, sep='\t'))
HTseqfiles_names <- list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=F)
colnames(HTseq) <- HTseqfiles_names

smCountTable.1a <- HTseq[,grep("Fib", colnames(HTseq))]
colnames(smCountTable.1a)

smCountTable.1 <- cbind(smCountTable.1a[,grep("O_3F_Fib", colnames(smCountTable.1a))], smCountTable.1a[,grep("O_Fib", colnames(smCountTable.1a))])
colnames(smCountTable.1)

txdb <- makeTxDbFromGFF("../../../22_JY_RNA_seq/old_fib_0923/mef_2/genes.gtf", format="gtf")

ebg <- exonsBy(txdb, by="gene")
ebg

exonic.gene.sizes <- sum(width(reduce(ebg)))


ind <- match(rownames(smCountTable.1), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
length <- data.frame(gene.symbol=rownames(smCountTable.1),
                     exonic.size=exonic.gene.sizes.ord)
length1<- as.numeric(length[,2])

fpkm <- rpkm(smCountTable.1, gene.length=length1,normalized.lib.sizes=F, log=FALSE)
keep <- rowSums(fpkm > 0.2) >= 1
smCountTable.filt <- smCountTable.1[keep,] # maintain 14310 genes for further analysis

condition = c(rep('O_3F_Fib',3), rep('O_Fib',3))
condition

smCountTable.filt

exp.data <- data.frame(
  row.names=colnames(smCountTable.1),
  condition = condition)

exp.data

dds <- DESeqDataSetFromMatrix(countData = smCountTable.1,
                              colData = exp.data,
                              design = ~ condition)
dds
dds <- DESeq(dds)

resultsNames(dds)

plot_total_counts(dds)
plot_library_complexity(dds)
vsd <- vst(dds)
mean_sd_plot(vsd)
plot_sample_clustering(vsd, anno_vars = c("condition"), distance = "euclidean")

dds <- estimateSizeFactors(dds)
dds <- sizeFactors(dds)
nor_counts <- counts(dds, normalized=TRUE)

nor_counts

resultsNames(dds)

res <- results(dds)
res

hist(res$pvalue) # Looks good!

res <- results(dds, cooksCutoff=T, independentFiltering = T)
# order by FDR

dds
res <- results(dds, contrast=c("condition","O_3F_Fib","O_Fib"))
res

resOrdered <- res1[order(res1$ pvalue),]
summary(resOrdered)

resSig <- subset(resOrdered, pvalue < 0.05) # 2563 DEGs at FDR < 0.05
summary(resSig)

resSIG_FC <- subset(resSig, abs(log2FoldChange) > 0.5) # 885 DEGs at FDR < 0.05 and 1.5 FC
summary(resSIG_FC)

write.table(res1, file="DESeq2_all.csv", sep=',')
write.table(resSig, file="egm_vs_old_mef.csv", sep=',')

# save DEG table
save(resOrdered, file = "DESeq2_DEG_all_genes.Rda")
