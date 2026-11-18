suppressPackageStartupMessages({
  library(Seurat); library(RColorBrewer); library(dplyr)
})

# ---- object, res 0.2 clusters --------------------------------------------------
seu2 <- readRDS("EGM_multiflex/seu_harmony_rmov12.rds")
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

# ---- sample 라벨 + factor 순서 -----------------------------
seu2$sample <- factor(
  dplyr::recode(as.character(seu2$sample), "A1"="A", "Y1"="Y"),
  levels = c("A", paste0("R", 0:9), "Y")
)

# ---- UMAP split by sample -------------------------------------------------------
DimPlot(seu2, reduction = "umap_harmony", label = FALSE,
        cols = brewer.pal(6, "Set2"), pt.size = 1.8,
        split.by = "sample", ncol = 4)
