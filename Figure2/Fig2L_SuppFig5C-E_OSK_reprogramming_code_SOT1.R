## ---- setup, normalisation, clustering, UMAP; DimPlot = Supp Fig 5D ----------
library(SOT)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)
library(SingleCellExperiment)
library(clusterProfiler)
library(ReactomePA)
library(destiny)
library(Matrix)
library(Seurat)
library(slingshot)


library(scater)
# install loomR from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/loomR', ref = 'develop')
library(loomR)
library(Seurat)
library(patchwork)

setwd("43_SOT/")

file_dir = system.file("extdata", "c1_osk.csv.gz", package = "SOT")
exprs = read.csv(file_dir, row.names = 1,check.names = FALSE)
exprs[1:5,1:5]
dim(exprs)


sce = BuildSCE(normcounts=exprs, trim = 5.1/ncol(exprs))

day = c("mef", "osk_d0", "osk_d1", "osk_d2", "osk_d3", "osk_d4", "osk_d5", "osk_d7", "osk_d8", "ips", "esc")
sce$Day = factor(gsub("(mef|osk_d0|osk_d1|osk_d2|osk_d3|osk_d4|osk_d5|osk_d7|osk_d8|ips|esc).*", "\\1", colnames(sce)),levels=day)
metadata(sce)$`Day color` = c("mef" = "#BEAED4",
                              "osk_d0" = "#fdac86",
                              "osk_d1" = "#FFFF99",
                              "osk_d2" = "#386CB0",
                              "osk_d3" = "#F0027F",
                              "osk_d5" = "#BF5B17",
                              "osk_d7" = "#666666",
                              "osk_d8" = "#ffd56f",
                              "ips" = "#7FBC41",
                              "esc" = "#A6CEE3")

sce <- runPCA(sce)

sce.seurat <- as.Seurat(sce, counts = "normcounts", data = "logcounts")
sce.seurat
Idents(sce.seurat) <- "Day"
p1 <- DimPlot(sce.seurat, reduction = "PCA", group.by = "Source") + NoLegend()
p2 <- RidgePlot(sce.seurat, features = "ACTB", group.by = "Source")
p1 + p2

DimPlot(sce.seurat)

sce.seurat$Day


sce.seurat <- NormalizeData(sce.seurat)
all.genes <- rownames(sce.seurat)
sce.seurat <- ScaleData(sce.seurat, features = all.genes)
sce.seurat <- FindVariableFeatures(object = sce.seurat)
sce.seurat <- RunPCA(sce.seurat, features = VariableFeatures(object = sce.seurat))
sce.seurat <- FindNeighbors(sce.seurat, dims = 1:5)
sce.seurat <- FindClusters(sce.seurat, resolution = 0.5)
sce.seurat <- RunUMAP(sce.seurat, dims = 1:5)
DimPlot(sce.seurat, reduction = "umap", label = F, cols = c("#ACA4C9", "#535B36", "#F7A58A", "#4E7BAB", "#EBEEA2", "#4E7BAB", "#B18454", "#E74A48", "#18B5B7", "#916367"), pt.size = 0.6)

## ---- iPSC + ESC markers (ipsc_esc.csv), the pool the oncogene list was taken from ----
Idents(sce.seurat) <- "Day"

ipsc_esc.markers <- FindMarkers(sce.seurat, ident.1 = c("ips", "esc"), logfc.threshold = 0.5, only.pos = T)


head(ipsc_esc.markers, n = 5)
write.csv(ipsc_esc.markers, "./ipsc_esc.csv")

## ---- EGM-Young (65-gene) module score; violin by time point (Fig 2M) --------
young_ear_signature <- as.matrix(read.csv("../41_RNA_seq_test/genelist_ear_uu.list.csv"))
young_ear_signature <- as.character(young_ear_signature[,2])

sce.seurat = AddModuleScore(
  object = sce.seurat,
  features = list(young_ear_signature),
  name = 'young_ear_signature', ctrl = 3)

VlnPlot(sce.seurat, features = "young_ear_signature1", pt.size = 0, split.by = "Day", cols = c("#ACA4C9", "#F7A58A", "#EBEEA2", "#4E7BAB", "#E74A48", "#916367", "#935DA7", "#F6C883", "#FF9F1C", "#7FBF70")) + geom_boxplot(outlier.size = 0, fill="white",width=0.1)

## ---- Fig 2L: module score on the UMAP ---------------------------------------
FeaturePlot(sce.seurat, features = "young_ear_signature1", cols = c("yellow", "red"), min.cutoff = 0, max.cutoff = 0.05, pt.size = 1)

## ---- Fig 2L: variant with max.cutoff 0.03 -----------------------------------
FeaturePlot(sce.seurat, features = "young_ear_signature1", cols = c("yellow", "red"), min.cutoff = 0, max.cutoff = 0.03, pt.size = 0.5)

## ---- iPSC/ESC-high oncogene list (Supp Fig 5E; same list in Supp 5F) --------
cancer_genes <- c("Kit", "Sf3b1", "Ccne1", "Mdm4", "Mycn", "Pim1", "Srsf2", "Xpo1", "Eif4a2", "Fgf4", "Etv5", "Hdac1", "Sgk1", "Cenpa", "Dek", "Hsp90aa1", "Ncoa3", "Nup98", "Set", "Cad", "Gab1", "Hsp90ab1", "Srsf3", "Bub1b", "Ccnb1ip1", "Clp1", "Cnbp", "Folr1", "Inpp5d", "Kat6b", "Kdm2b", "Mthfd2", "Naca", "Parp2", "Pum1", "Rnaseh2a", "Sall4", "Tcea1", "Tfrc", "Top2a")

## ---- Supp Fig 5E: oncogene module score by time point -----------------------
sce.seurat = AddModuleScore(
  object = sce.seurat,
  features = list(cancer_genes),
  name = 'cancer_genes', ctrl = 3)

VlnPlot(sce.seurat, features = "cancer_genes1", pt.size = 0, split.by = "Day", cols = c("#ACA4C9", "#F7A58A", "#EBEEA2", "#4E7BAB", "#E74A48", "#916367", "#935DA7", "#F6C883", "#FF9F1C", "#7FBF70"))+ geom_boxplot(outlier.size = 0, fill="white",width=0.1)

## ---- Slingshot trajectory; UMAP coloured by pseudotime = Supp Fig 5C --------
sort_NSC_1 <- sce.seurat

sce_1 <- as.SingleCellExperiment(sort_NSC_1)

shuffle <- sample(ncol(sce_1))

sce_1 <- slingshot(sce_1, reducedDim = 'UMAP', clusterLabels = sce_1$Day, start.clus = 'mef', approx_points = 200)
layout(matrix(c(1, 1, 2, 3), 2))
par(mar = c(4.5, 4, 1, 1))

plot(reducedDims(sce_1)$UMAP, col = cell_colors, pch=16, asp = 1)
lines(SlingshotDataSet(sce_1), lwd=2, col='black')

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

sce_1$Day

cell_colors <- cell_pal(sce_1$condition, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(sce_1$seurat_clusters, hue_pal())

plot(reducedDims(sce_1)$UMAP[shuffle, ], asp = 1, pch = 16, xlab = "UMAP_1", ylab = "UMAP_2",
     col = hcl.colors(100, alpha = .5)[cut(sce_1$slingPseudotime_1, breaks = 100)][shuffle])

## ---- pseudotimeED (stored in sce_seurat.RData; used by Fig 2N) --------------
library("tradeSeq")

sce_2 <- sce_1

slingsce<-SlingshotDataSet(sce_2)


pseudotimeED <- slingPseudotime(slingsce, na = FALSE)

cellWeightsED <- slingCurveWeights(slingsce)
Days <- colData(sce_2)$Day

