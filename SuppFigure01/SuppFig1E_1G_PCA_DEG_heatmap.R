## ---------------------------------------------------------------------------
## Part A - read counting and DESeq2
## ---------------------------------------------------------------------------
library(Rsubread)
library(DESeq2)
library("Rsamtools")
library("GenomicFeatures")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("GenomicAlignments")
library("BiocParallel")
library("edgeR")
library(gplots)
library(pheatmap)


setwd("22_JY_RNA_seq")


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

sample <- c("Neg_P","Neg_P","Neg_P","4FP_P","4FP_P","4FP_P")

sample <- as.data.frame(sample)

colData(se) <- DataFrame(sample)

colData(se)

se$sample

countdata <- assay(se)


head(countdata, 3)

dds <- DESeqDataSet(se, design = ~ sample)

dds <- DESeq(dds)

res <- results(dds, contrast=c("sample","4FP_P","Neg_P"))
res


mcols(res, use.names=TRUE)
summary(res)

res.05 <- results(dds, alpha=.1)
table(res.05$padj < 0.05)
summary(res.05)

## ---------------------------------------------------------------------------
## Part B - packages
## ---------------------------------------------------------------------------
library('edgeR')    # v3.12.1
library("DESeq2")   # v1.20.0
library('limma')
library('pheatmap') # v3.26.9

## ---------------------------------------------------------------------------
## Part B - Supp Fig 1G, DEG heatmap
## ---------------------------------------------------------------------------
res.05

dds$sample

resxxx <- results(dds, contrast=c("sample","4FP_P","Neg_P"))

resxxx

dds$sample

vst <- varianceStabilizingTransformation(dds, blind =F)
vstMat <- assay(vst)
vstMat
# batch-correct using limma function
v.data <- removeBatchEffect(vstMat)
v.data
genes<- as.data.frame(resxxx)
# Make gene names in CAPTIAL letter
row.names(genes) <- toupper(row.names(genes))
# Extract all DEGs that are FDR <0.05 and >1.5 in Fold Change, higher in young
deg <- genes[which(genes$pvalue < 0.05 & abs(genes$log2FoldChange) > 0),] # 1022 DEGs

#### Extract expression values (vst) for all DEGs between young and old that goes up with age
# Upload vst values from DESeq2

# Extract vst values for DEGs
vst.deg <-as.data.frame(vstMat[which(toupper(row.names(vstMat)) %in% row.names(deg)),])

vst.deg

# Generate heatmap
pdf("Heatmap_DESeq2_BATCH_Corrected_vst_DEGs_fibs_in_vitro_080618.pdf", width =4.3, height =5, onefile=F)
par(mar=c(4.1,4.1,6,1), xpd=TRUE)
pheatmap(vst.deg,
         border_color = "black",
         cluster_col=F,
         cluster_row=T,
         scale="row",
         annotation_legend=TRUE,
         clustering_distance_cols="manhattan",
         clustering_method="complete",
         fontsize_col=8,
         fontsize_row=3.5,
         show_rownames=F,
         show_colnames=T,
         fontface="bold",
         treeheight_col=10,
         treeheight_row=0,
         cellwidth = NA,
         cellheight = NA,
         color = colorRampPalette(c("#114D8A","#DBEFF6", "#FEFCBA","#FDCB82","#D22B28"))(50))
dev.off()

## ---------------------------------------------------------------------------
## Part B - Supp Fig 1E, PCA
## ---------------------------------------------------------------------------
colnames(v.data)

# Change the names so that it reflects their corresponding RE propenisty
colnames(v.data) <- c(rep("4FP",3), rep("Plasma",3))
colnames(v.data)

v.data

# PCA using prcomp function
e.pca <- prcomp(t(v.data))

# get the precentage for each component
summary(e.pca)
e.pca

# Importance of components:
#                         PC1     PC2     PC3     PC4      PC5      PC6    PC7      PC8      PC9    PC10     PC11     PC12     PC13     PC14     PC15      PC16     PC17     PC18
# Standard deviation     55.341 41.6419 38.6432 34.5661 31.43161 29.29025 27.278 26.37877 24.96769 24.3679 23.05823 22.62491 21.29661 21.03624 19.56099 7.948e-13 4.84e-13 2.78e-13
# Proportion of Variance  0.214  0.1212  0.1043  0.0835  0.06904  0.05995  0.052  0.04863  0.04356  0.0415  0.03715  0.03577  0.03169  0.03092  0.02674 0.000e+00 0.00e+00 0.00e+00
# Cumulative Proportion   0.214  0.3352  0.4395  0.5231  0.59208  0.65204  0.704  0.75266  0.79622  0.8377  0.87487  0.91064  0.94234  0.97326  1.00000 1.000e+00 1.00e+00 1.00e+00

# find the color code for orange and darkblue
col2rgb("orange")   # 255, 165, 0
col2rgb("darkblue") # 0, 0, 139

# Make PCA plot based PC1 and PC2
pdf("PCA_DESeq2_BATCH_Corrected_vst_RNAseq_Fib_PC1-2.no.label.pdf", width =4.3, height =5 )
par(mar=c(4.1,4.1,6,1), xpd=TRUE)
plot(e.pca$x[,c(1,2)],
     pch = 15,
     xlab = 'PC1 (21%)',
     ylab = 'PC2 (12%)',
     cex.lab = 1.5,
     col = "white",
     cex=0.6,
     xlim = c(-10,10),
     ylim = c(-50,50))
points(e.pca$x[grep("4FP", colnames(v.data)) ,c(1,2)],  pch=5, cex=2,
       col=rgb(red=96,green=152,blue=255, alpha = 204, maxColorValue = 255))

points(e.pca$x[grep("Plasma", colnames(v.data)) ,c(1,2)],pch=5, cex=2,
       col=rgb(red=175,green=171,blue=171, alpha = 204, maxColorValue = 255))

e.pca$x

legend("topright",c("Old", "EGM"), inset=c(0,0),horiz=T,
       col=c("orange","Lightblue","darkblue"), text.col=c("orange","Lightblue","darkblue"), pch=c(17,17,17), cex = 0.3, pt.cex=1.5)
