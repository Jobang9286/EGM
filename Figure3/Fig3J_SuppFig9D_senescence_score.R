PY_CELLORACLE <- path.expand("~/.conda/envs/celloracle_env/bin/python")   # env with anndata

# ---- object + diffusion-map clusters ------------------------------------------
suppressPackageStartupMessages({
library(Seurat); library(slingshot); library(SingleCellExperiment)
library(ggplot2); library(dplyr); library(RColorBrewer)
})
options(Seurat.object.assay.version = "v3")
base_dir <- "EGM_multiflex/oracle"
seu <- readRDS(file.path(base_dir, "oracle_input_hvg.rds"))
# DM 위에서 Leiden res 0.1 재계산
seu <- FindNeighbors(seu, reduction = "diffmap", dims = 1:5,
graph.name = c("dm_nn", "dm_snn"), verbose = FALSE)
seu <- FindClusters(seu, graph.name = "dm_snn",
resolution = 0.1, algorithm = 4,
leiden_method = "igraph", random.seed = 42,
cluster.name = "dm_leiden_0.1", verbose = FALSE)
seu$dm_cluster <- factor(seu$dm_leiden_0.1)
cat("dm_cluster regenerated, levels:", levels(seu$dm_cluster), "\n")

# ---- senescence (senepy hub0) score + Supp 9D ---------------------------------
suppressPackageStartupMessages({
library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
library(viridisLite); library(reticulate); library(scales)
})
options(Seurat.object.assay.version = "v3")
base_dir <- "EGM_multiflex/oracle"
# ---- 1. Senepy gene list 읽기 (weighted) -----------------------
sen_tbl <- read.table(
file.path(base_dir, "phase3_figs", "hub0_genelist.txt"),
sep = "\t", header = FALSE, col.names = c("gene", "score"),
stringsAsFactors = FALSE)
cat("Total senepy genes:", nrow(sen_tbl), "\n")
cat("Score range:", range(sen_tbl$score), "\n")
# 상위 100, 300, all 버전 비교
top100 <- head(sen_tbl, 100)
top300 <- head(sen_tbl, 300)
# ---- 2. Full-gene h5ad에서 expression 가져와서 weighted score 계산 ----
Sys.setenv(RETICULATE_PYTHON = PY_CELLORACLE)
use_python(PY_CELLORACLE, required = TRUE)
ad <- import("anndata", convert = FALSE)
adata_full <- ad$read_h5ad(file.path(base_dir, "seu_R0R9_for_oracle_with_raw.h5ad"))
all_genes <- py_to_r(adata_full$var_names$values)
# 사용 가능한 senepy gene
sen_in_data <- sen_tbl[sen_tbl$gene %in% all_genes, ]
cat("Senepy genes in adata:", nrow(sen_in_data), "/", nrow(sen_tbl), "\n")
# logcounts 가져오기
sub_ad <- adata_full[, sen_in_data$gene]
expr <- as.matrix(py_to_r(sub_ad$layers["logcounts"]$toarray()))
rownames(expr) <- py_to_r(sub_ad$obs_names$values)
colnames(expr) <- sen_in_data$gene
# 가중 평균 (hub score를 weight로)
w <- sen_in_data$score
w <- w / sum(w)
senepy_score <- expr[colnames(seu), ] %*% w
seu$senepy_aged <- as.numeric(senepy_score)
# Seurat AddModuleScore 버전도 동시 (control-corrected)
seu <- AddModuleScore(seu,
features = list(senepy_top300 = sen_in_data$gene[1:300]),
name = "senepy_top300_", nbin = 24, ctrl = 100)
cat("senepy_aged (weighted) range:", round(range(seu$senepy_aged), 3), "\n")
cat("senepy_top300_1 (Seurat) range:", round(range(seu$senepy_top300_1), 3), "\n")
# ---- 3. DM 위에 visualization ----------------------------------
dc <- Embeddings(seu, "diffmap")
df <- data.frame(
DC1 = dc[, 1], DC2 = dc[, 2],
sample = factor(seu$sample, levels = paste0("R", 0:9)),
senepy = seu$senepy_aged,
senepy_seu = seu$senepy_top300_1,
pt = seu$palantir_pseudotime,
cluster = seu$dm_cluster
)
# (a) DM with senepy score                                   <- Supp Fig 9D
p_dm <- ggplot(df[order(df$senepy), ], aes(DC1, DC2, color = senepy)) +
geom_point(size = 0.4, alpha = 0.85) +
scale_color_gradient2(low = "#2c5599", mid = "white", high = "#c0392b",
midpoint = median(df$senepy),
name = "senepy\naged") +
theme_classic() + coord_fixed() +
labs(title = "Senepy aged signature on DM")
# (b) per-sample density
library(ggridges)
p_ridge <- ggplot(df, aes(senepy, sample, fill = sample)) +
geom_density_ridges(scale = 1.5, alpha = 1) +
scale_fill_manual(values = viridis(10)) +
theme_bw() + theme(legend.position = "none") +
labs(x = "Senepy aged score", y = NULL,
title = "Senepy aged score per sample")
# (c) score vs pseudotime
p_pt <- ggplot(df, aes(pt, senepy)) +
geom_point(size = 0.3, alpha = 0.4, color = "#888") +
geom_smooth(method = "gam", formula = y ~ s(x, k = 8),
color = "#c0392b", linewidth = 1.2) +
theme_bw() +
labs(x = "Palantir pseudotime", y = "Senepy aged score",
title = "Senepy ~ pseudotime")
print((p_dm | p_ridge) / p_pt)

suppressPackageStartupMessages({
library(ggplot2); library(dplyr)
})
base_dir <- "EGM_multiflex/oracle"
df <- data.frame(
pt = seu$palantir_pseudotime,
senepy = seu$senepy_aged
)

# ---- Fig 3J --------------------------------------------------------------------
p_senepy_clean <- ggplot(df, aes(pt, senepy)) +
geom_smooth(method = "loess", span = 0.6,
color = "#2c5599", fill = "#2c5599",
alpha = 0.4, linewidth = 1.2, level = 0.99) +
scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
theme_classic(base_size = 12) +
theme(plot.title = element_text(face = "bold", hjust = 0,
size = 13, lineheight = 1.2),
plot.title.position = "plot",
axis.title.x = element_blank(),
axis.line = element_line(linewidth = 0.4),
axis.ticks = element_line(linewidth = 0.4),
panel.grid = element_blank()) +
labs(title = "Gene Signature of Senescent fibroblast\n(senepy hub0, ~800 genes)",
y = "Senepy aged score")
print(p_senepy_clean)
