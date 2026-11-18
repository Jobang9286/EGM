# ---- libraries, inputs ----
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ComplexHeatmap); library(circlize); library(ggplot2); library(patchwork)})
root    <- "EGM_project_for_deconvolution"
deg_dir <- file.path(root, "integrated_analysis_for_lizard_factors/kidney/deg_results/cell_type_specific")
cell_types <- c("PT","Endo","LOH","DCT","CD-Trans","Podo","CD-PC","Macro")
deg_list <- setNames(lapply(cell_types, function(ct)
  read.csv(file.path(deg_dir, paste0(ct, "_DEG_HGPS_EGM_vs_control.csv")))), cell_types)
load(file.path(root, "integrated_analysis_for_lizard_factors/kidney/deg_results_2/kidney_deconvolution_main_objects.RData"))
# ---- p-value label ----
fmt_p <- function(p) if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)
hg <- sample_metadata$Sample
# ---- t-test helper (raw and log2 scale) ----
ttest_both <- function(Z, g) {
  e <- Z[g, hg]; l <- log2(e + 1)
  data.frame(gene = g, Mean_Control = mean(e[1:3]), Mean_EGM = mean(e[4:6]),
             P_raw = t.test(e[4:6], e[1:3])$p.value, P_log2 = t.test(l[4:6], l[1:3])$p.value)
}
Z_PT_log <- log2(Z_PT_full + 1)
Z_Endo_log <- log2(Z_Endo_full + 1)
# ---- gene choice and t-tests ----
pt_v   <- c("Cdkn1b", "Igfbp7", "Ccnb1", "Cdk1")
endo_v <- c("Cdkn1a", "Cdkn2b", "Ccne1", "Mcm2")
stopifnot(all(pt_v %in% rownames(Z_PT_full)), all(endo_v %in% rownames(Z_Endo_full)))
# t-test는 Z 행렬에서 직접 계산 (Ccnb1/Cdk1/Ccne1/Mcm2는 aging set 밖이라 pt_all/endo_all에 없음)
tt_pt   <- bind_rows(lapply(pt_v,   function(g) ttest_both(Z_PT_full, g)))   %>% left_join(deg_list$PT   %>% select(gene, log2FoldChange, padj), by = "gene")
tt_endo <- bind_rows(lapply(endo_v, function(g) ttest_both(Z_Endo_full, g))) %>% left_join(deg_list$Endo %>% select(gene, log2FoldChange, padj), by = "gene")
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
# ---- row builder ----
make_row <- function(Z, genes, pvals, ylab, strip = FALSE, fun = plot_gene_violin_clean) {
  lapply(seq_along(genes), function(i) {
    p <- fun(Z, genes[i], sample_metadata, pval = pvals[i], ylab = ylab)
    if (strip) p <- strip_text(p)
    if (i != 1) p <- p + theme(axis.title.y = element_blank())   # strip 이후에 적용해야 텍스트 버전과 동일 레이아웃
    p })
}
# 원인 분리용 변형: A = strip만(행 제목은 검정), B = 행 제목만 투명
build_variant <- function(strip_panels, title_transparent, fun = plot_gene_violin_clean) {
  tcol <- if (title_transparent) "transparent" else "black"
  wrap_plots(c(list(row_title("Aging markers in Proximal tubular cells", tcol)),
               make_row(Z_PT_log, pt_v, tt_pt$P_log2, "Expression (log2 counts + 1)", strip_panels, fun),
               list(row_title("Aging markers in Endothelial cells", tcol)),
               make_row(Z_Endo_log, endo_v, tt_endo$P_log2, "Expression (log2 counts + 1)", strip_panels, fun)),
             design = design, heights = c(0.1, 1, 0.1, 1))
}
# ---- y-axis rule, theme, boxplot function; p_box / p_box_nt ----
# y축 눈금: 최대 3개, 0.5 단위로만 떨어지도록
breaks_half <- function(lims) {
  steps <- c(0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 15, 20)
  for (s in steps) { b <- seq(ceiling(lims[1] / s) * s, floor(lims[2] / s) * s, by = s); if (length(b) <= 3) return(b) }
  b
}
base_theme <- function() theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 16), plot.caption = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 11),
        axis.title.y = element_text(size = 12), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1), plot.margin = panel_margin)
gene_df <- function(Z_matrix, gene) {
  gene_expr <- Z_matrix[gene, ]
  data.frame(Sample = names(gene_expr), Expression = gene_expr,
             Group = sample_metadata$Group[match(names(gene_expr), sample_metadata$Sample)]) %>%
    filter(Group %in% c("HGPS_control", "HGPS_EGM")) %>% mutate(Group = factor(Group, levels = c("HGPS_control", "HGPS_EGM")))
}
plot_gene_box_clean <- function(Z_matrix, gene, sample_metadata, pval = NULL, ylab = "Expression (log2 TPM + 1)") {
  if (!gene %in% rownames(Z_matrix)) return(NULL)
  ggplot(gene_df(Z_matrix, gene), aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(width = 0.6, colour = "black", linewidth = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = point_size, alpha = 0.8) +
    scale_fill_manual(values = cond_colors) + scale_x_discrete(labels = c("Control", "EGM")) +
    scale_y_continuous(breaks = breaks_half) +
    labs(title = gene, x = "", y = ylab, caption = if (!is.null(pval)) fmt_p(pval) else NULL) + base_theme()
}
plot_gene_violin_clean <- function(Z_matrix, gene, sample_metadata, pval = NULL, ylab = "Expression (log2 TPM + 1)") {
  if (!gene %in% rownames(Z_matrix)) return(NULL)
  ggplot(gene_df(Z_matrix, gene), aes(x = Group, y = Expression, fill = Group)) +
    geom_violin(alpha = 1, trim = FALSE, colour = "black", linewidth = 0.5) + geom_jitter(width = 0.1, size = point_size, alpha = 0.8) +
    scale_fill_manual(values = cond_colors) + scale_x_discrete(labels = c("Control", "EGM")) +
    scale_y_continuous(breaks = breaks_half) +
    labs(title = gene, x = "", y = ylab, caption = if (!is.null(pval)) fmt_p(pval) else NULL) + base_theme()
}
p_box <- build_variant(FALSE, FALSE, plot_gene_box_clean); p_box_nt <- build_variant(TRUE, TRUE, plot_gene_box_clean)
# ---- output folder ----
final_dir <- file.path(root, "integrated_analysis_for_lizard_factors/integrated_analysis/kidney/final_aging")
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)
# ---- kidney violin/box 최종 저장 ----
ggsave(file.path(final_dir, "kidney_aging_boxplot_PT_Endo_FINAL.pdf"),        p_box,    width = 12, height = 8, dpi = 300)
ggsave(file.path(final_dir, "kidney_aging_boxplot_PT_Endo_FINAL.png"),        p_box,    width = 12, height = 8, dpi = 300)
ggsave(file.path(final_dir, "kidney_aging_boxplot_PT_Endo_FINAL_notext.pdf"), p_box_nt, width = 12, height = 8, dpi = 300)
stats_tbl <- bind_rows(tt_pt %>% mutate(CellType = "PT"), tt_endo %>% mutate(CellType = "Endo")) %>%
  select(CellType, gene, Mean_Control, Mean_EGM, log2FoldChange, padj_DESeq2 = padj, P_ttest_raw = P_raw, P_ttest_log2 = P_log2)
write.csv(stats_tbl, file.path(final_dir, "kidney_aging_boxplot_genes_stats.csv"), row.names = FALSE)
