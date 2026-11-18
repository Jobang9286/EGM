# ---- load object, cluster (res 0.2) -------------------------------------------
seu2 <- readRDS("EGM_multiflex/seu_harmony_rmov12.rds")
library(Seurat)
library(Matrix)
library(reticulate)
seu2 <- FindClusters(seu2, resolution = 0.2)

# ---- merge clusters 3 and 4 (0-6 -> 0-5) ---------------------------------------
remap <- c("0"="0","1"="1","2"="2","3"="3","4"="3","5"="4","6"="5")
seu2$seurat_clusters <- factor(
unname(remap[as.character(seu2$RNA_snn_res.0.2)]),
levels = c("0","1","2","3","4","5")
)
Idents(seu2) <- "seurat_clusters"

library(RColorBrewer)

# 0~5 → 1~6 리넘버
seu2$seurat_clusters <- factor(
as.integer(as.character(seu2$seurat_clusters)) + 1,
levels = 1:6
)
Idents(seu2) <- "seurat_clusters"

# ---- Fig 3B --------------------------------------------------------------------
# 시각화
DimPlot(seu2, reduction = "umap_harmony", label = F,
cols = brewer.pal(6, "Set2"), pt.size = 1.2)

# ---- Fig 3C --------------------------------------------------------------------
DimPlot(seu2, reduction = "umap_harmony", label = F,
cols = brewer.pal(6, "Set2"), pt.size = 1.2, split.by = "group")
