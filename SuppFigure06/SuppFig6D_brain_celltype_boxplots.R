# ---- libraries, root ----
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ComplexHeatmap); library(circlize); library(ggplot2); library(patchwork)})
root    <- "EGM_project_for_deconvolution"
# ---- brain DESeq2 results ----
brain_dir <- file.path(root, "integrated_analysis_for_lizard_factors/brain")
bd <- readRDS(file.path(brain_dir, "DEG_results_InstaPrism/Brain_DEG_Complete_Results.rds"))
# ---- cell types ----
brain_ct <- c("Excitatory Neurons","Inhibitory Neurons","Astrocytes","Microglia","ODC")
# ---- per-cell-type DESeq2 tables ----
brain_deg <- setNames(lapply(brain_ct, function(ct) bd$deg_results[[ct]]$full_results %>% rename(gene = Gene)), brain_ct)
# ---- genes per cell type ----
need <- list(`Excitatory Neurons` = c("Igfbp4","Igfbp3","Nrxn1","Nrp1"), ODC = c("Cdkn1a","Cdkn1c","Cryab","Kl"), Astrocytes = c("Cd68","Cd74","Nfkbia","Gria1"))
# ---- dds objects ----
suppressPackageStartupMessages(library(DESeq2))
dds_list <- setNames(lapply(brain_ct, function(ct) bd$deg_results[[ct]]$dds), brain_ct)
# ---- DESeq2-normalised counts and groups ----
Zb_norm <- setNames(lapply(names(need), function(ct) log2(counts(dds_list[[ct]], normalized = TRUE) + 1)), names(need))
brain_meta <- data.frame(Sample = colnames(Zb_norm[[1]]), Group = as.character(colData(dds_list[[1]])$Group), stringsAsFactors = FALSE)
# ---- padj lookup, 3-row layout ----
padj_of <- function(ct, g) { d <- brain_deg[[ct]]; d$padj[match(g, d$gene)] }
design3 <- c(area(1,1,1,4), area(2,1), area(2,2), area(2,3), area(2,4),
             area(3,1,3,4), area(4,1), area(4,2), area(4,3), area(4,4),
             area(5,1,5,4), area(6,1), area(6,2), area(6,3), area(6,4))
# ---- condition colours ----
cond_colors <- c("HGPS_control" = "#D9D9D9", "HGPS_EGM" = "#7FA6FF")   # 참고 이미지: Control 연회색, EGM 연파랑
# ---- point size, panel margin ----
point_size <- 1.6 * 0.7    # 30% 축소 -> 1.12
panel_margin <- margin(2, 1, 2, 1)   # plot 간 간격 축소 (top, right, bottom, left; pt)
tr <- "transparent"
# ---- strip_text ----
strip_text <- function(p) p + theme(plot.title = element_text(colour = tr), plot.caption = element_text(colour = tr),
                                    axis.title.x = element_text(colour = tr), axis.title.y = element_text(colour = tr),
                                    axis.text.x = element_text(colour = tr), axis.text.y = element_text(colour = tr))
# ---- row titles and 2-row patchwork layout ----
library(grid)
row_title <- function(txt, col = "black") wrap_elements(full = textGrob(txt, gp = gpar(fontsize = 15, fontface = "bold", col = col)))
design <- c(area(1, 1, 1, 4), area(2, 1), area(2, 2), area(2, 3), area(2, 4),
            area(3, 1, 3, 4), area(4, 1), area(4, 2), area(4, 3), area(4, 4))
# ---- base theme ----
base_theme <- function() theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 16), plot.caption = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1), plot.margin = panel_margin)
# ---- gene_df with explicit metadata ----
gene_df <- function(Z_matrix, gene, meta) {
  gene_expr <- Z_matrix[gene, ]
  data.frame(Sample = names(gene_expr), Expression = gene_expr,
             Group = meta$Group[match(names(gene_expr), meta$Sample)]) %>%
    filter(Group %in% c("HGPS_control", "HGPS_EGM")) %>% mutate(Group = factor(Group, levels = c("HGPS_control", "HGPS_EGM")))
}
# ---- p label ----
fmt_p <- function(p) if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)
# ---- axis expansion constants ----
EXP_MULT <- 0.20     # 축 확장
SEARCH   <- 0.10     # 눈금 탐색 범위 (확장보다 작게 -> 계산된 눈금이 항상 패널 안에 들어옴)
# ---- y-axis rule ----
y_scale_half <- function(v) {
  r <- range(v); span <- diff(r); ex <- r + c(-1, 1) * span * SEARCH
  b <- NULL
  for (s in c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 15, 20)) {
    bb <- seq(ceiling(ex[1] / s) * s, floor(ex[2] / s) * s, by = s)
    if (length(bb) >= 2 && length(bb) <= 3) { b <- bb; break }
    if (length(bb) < 2) break
  }
  if (is.null(b)) b <- seq(floor(r[1] / 0.5) * 0.5, ceiling(r[2] / 0.5) * 0.5, by = 0.5)  # 범위가 좁을 때
  lim <- c(min(r[1], min(b)), max(r[2], max(b)))          # 눈금이 밖에 있으면 축을 그만큼 넓힘
  scale_y_continuous(breaks = b, limits = lim, expand = expansion(mult = 0.06))
}
# ---- boxplot function (uses y_scale_half) ----
plot_gene_box_clean <- function(Z_matrix, gene, meta, pval = NULL, ylab = "Expression (log2 TPM + 1)") {
  if (!gene %in% rownames(Z_matrix)) return(NULL)
  d <- gene_df(Z_matrix, gene, meta)
  ggplot(d, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(width = 0.6, colour = "black", linewidth = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = point_size, alpha = 0.8) +
    scale_fill_manual(values = cond_colors) + scale_x_discrete(labels = c("Control", "EGM")) +
    y_scale_half(d$Expression) +
    labs(title = gene, x = "", y = ylab, caption = if (!is.null(pval)) fmt_p(pval) else NULL) + base_theme()
}
# ---- row builder on normalised counts ----
make_row_b <- function(ct, genes, strip = FALSE) lapply(seq_along(genes), function(i) {
  p <- plot_gene_box_clean(Zb_norm[[ct]], genes[i], brain_meta, pval = padj_of(ct, genes[i]), ylab = "Expression (log2 normalized counts)")
  if (strip) p <- strip_text(p); if (i != 1) p <- p + theme(axis.title.y = element_blank()); p })
# ---- text version p_brain ----
build_brain_box <- function(strip = FALSE, row_gap = 0.1) {
  tcol <- if (strip) "transparent" else "black"
  wrap_plots(c(list(row_title("Excitatory neurons", tcol)), make_row_b("Excitatory Neurons", need$`Excitatory Neurons`, strip),
               list(row_title("Oligodendrocytes (ODC)", tcol)), make_row_b("ODC", need$ODC, strip),
               list(row_title("Astrocytes", tcol)), make_row_b("Astrocytes", need$Astrocytes, strip)),
             design = design3, heights = c(row_gap, 1, row_gap, 1, row_gap, 1))
}
p_brain    <- build_brain_box(FALSE, 0.1)
p_brain_nt <- build_brain_box(TRUE, 0.02)          # 텍스트 없는 버전: 행 제목 자리를 최소화
# ---- output folder ----
tmp_dir <- "output"
# ---- no-text version ----
strip_text_tight <- function(p) p + theme(
  plot.title   = element_text(colour = tr, size = 16, face = "bold.italic", hjust = 0.5, margin = margin(b = 0)),
  plot.caption = element_text(colour = tr, size = 12, hjust = 0.5, margin = margin(t = 0)),
  axis.title.x = element_text(colour = tr, margin = margin(t = 0)), axis.title.y = element_text(colour = tr, margin = margin(r = 0)),
  axis.text.x  = element_text(colour = tr, angle = 45, hjust = 1, size = 12, margin = margin(t = 0)),
  axis.text.y  = element_text(colour = tr, size = 11, margin = margin(r = 0)))
make_row_b2 <- function(ct, genes) lapply(seq_along(genes), function(i) {
  p <- strip_text_tight(plot_gene_box_clean(Zb_norm[[ct]], genes[i], brain_meta, pval = padj_of(ct, genes[i]), ylab = "Expression (log2 normalized counts)"))
  if (i != 1) p <- p + theme(axis.title.y = element_blank()); p })
p_brain_nt <- wrap_plots(c(list(row_title("", tr)), make_row_b2("Excitatory Neurons", need$`Excitatory Neurons`),
                           list(row_title("", tr)), make_row_b2("ODC", need$ODC),
                           list(row_title("", tr)), make_row_b2("Astrocytes", need$Astrocytes)),
                         design = design3, heights = c(0.01, 1, 0.01, 1, 0.01, 1))
ggsave(file.path(tmp_dir, "brain_box4_notext.png"), p_brain_nt, width = 12, height = 9, dpi = 110)
