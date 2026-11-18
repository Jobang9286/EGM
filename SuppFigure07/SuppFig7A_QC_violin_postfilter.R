suppressPackageStartupMessages({
  library(Seurat); library(dplyr)
})

# ---- post-QC object ------------------------------------------------------------
seu2 <- readRDS("EGM_multiflex/seu_harmony_rmov12.rds")

# ---- sample 라벨 + factor 순서 -----------------------------
seu2$sample <- factor(
  dplyr::recode(as.character(seu2$sample), "A1"="A", "Y1"="Y"),
  levels = c("A", paste0("R", 0:9), "Y")
)

# ---- percent.mt / percent.rb ----------------------------------------------------
seu2[["percent.mt"]] <- PercentageFeatureSet(seu2, pattern = "^mt-")
seu2[["percent.rb"]] <- PercentageFeatureSet(seu2, pattern = "^RP[SL]")

# ---- 6. QC 필터링 (임계값은 분포 보고 조정) --------------
seu2 <- subset(seu2,
               subset = nFeature_RNA >= 500 &
                        nFeature_RNA <= 7500 &
                        percent.mt   <  10
)

# ---- QC violin plots with thresholds -------------------------------------------
library(patchwork)
library(ggplot2)
p <- VlnPlot(seu2,
             features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             group.by = "sample", pt.size = 0, ncol = 3,
             combine = FALSE)
# 각 panel에 임계선 추가
p[[1]] <- p[[1]] + geom_hline(yintercept = c(500, 7500),
                              color = "red", linetype = "dashed", linewidth = 0.6)
p[[3]] <- p[[3]] + geom_hline(yintercept = 10,
                              color = "red", linetype = "dashed", linewidth = 0.6)
wrap_plots(p, ncol = 1)
