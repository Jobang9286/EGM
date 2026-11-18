# ---- object, res 0.2 clusters --------------------------------------------------
seu2 <- readRDS("EGM_multiflex/seu_harmony_rmov12.rds")
library(Seurat)
seu2 <- FindClusters(seu2, resolution = 0.2)
# ---- merge clusters 3 and 4 (0-6 -> 0-5) ---------------------------------------
remap <- c("0"="0","1"="1","2"="2","3"="3","4"="3","5"="4","6"="5")
seu2$seurat_clusters <- factor(
  unname(remap[as.character(seu2$RNA_snn_res.0.2)]),
  levels = c("0","1","2","3","4","5")
)
Idents(seu2) <- "seurat_clusters"
# 0~5 → 1~6 리넘버
seu2$seurat_clusters <- factor(
  as.integer(as.character(seu2$seurat_clusters)) + 1,
  levels = 1:6
)
Idents(seu2) <- "seurat_clusters"

library(ggplot2)
library(dplyr)
library(viridisLite)
# ---- 0. (이미 했으면 skip) sample 라벨 + factor 순서 -----------
seu2$sample <- factor(
  dplyr::recode(as.character(seu2$sample), "A1"="A", "Y1"="Y"),
  levels = c("A", paste0("R", 0:9), "Y")
)
# ---- 1. 비율 데이터 ------------------------------------------
sample_order <- c("A", paste0("R", 0:9), "Y")
df <- as.data.frame(table(Cluster = seu2$seurat_clusters,
                          Sample  = seu2$sample)) %>%
  group_by(Cluster) %>%
  mutate(Fraction = Freq / sum(Freq) * 100) %>%
  ungroup()
df$Sample <- factor(df$Sample, levels = sample_order)
# ---- 2. 팔레트 -----------------------------------------------
pal <- c(A = "gray",
         setNames(viridis(10), paste0("R", 0:9)),
         Y = "#dee")

# ---- A, Y segments outlined in red dashed lines ----------
border_col <- setNames(ifelse(sample_order %in% c("A","Y"), "red", NA),
                       sample_order)
border_lty <- setNames(ifelse(sample_order %in% c("A","Y"), "dashed", "solid"),
                       sample_order)
ggplot(df, aes(x = Cluster, y = Fraction,
               fill = Sample, color = Sample, linetype = Sample)) +
  geom_col(width = 0.85, linewidth = 0.5) +
  scale_fill_manual(values = pal, breaks = sample_order) +
  scale_color_manual(values = border_col, breaks = sample_order, guide = "none") +
  scale_linetype_manual(values = border_lty, breaks = sample_order, guide = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.5)) +
  labs(x = "Clusters", y = "Sample Fraction (%)") +
  theme_classic(base_size = 13) +
  theme(legend.title    = element_blank(),
        legend.key.size = unit(0.4, "cm"))
