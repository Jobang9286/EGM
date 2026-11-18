# ---- libraries, root ----
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ComplexHeatmap); library(circlize); library(ggplot2); library(patchwork)})
root    <- "EGM_project_for_deconvolution"
# ---- zone DESeq2 results and Z matrices ----
liver_dir <- file.path(root, "integrated_analysis_for_lizard_factors/liver/liver_deg_results")
liver_deg <- readRDS(file.path(liver_dir, "deg_results_list.rds"))
Z_liver <- readRDS(file.path(liver_dir, "Z_celltype_expression.rds"))
# ---- periportal matrix and sample groups ----
Z_peri     <- Z_liver$Periportal
Z_peri_log <- log2(Z_peri + 1)
liver_meta <- data.frame(Sample = colnames(Z_peri), Group = rep(c("HGPS_control", "HGPS_EGM"), each = 3))
# ---- periportal genes, pericentral matrix, t-test helper ----
portal_v  <- c("Cdkn1a", "Cdkn1c", "E2f1", "Yap1")
Z_cent <- Z_liver$Pericentral; Z_cent_log <- log2(Z_cent + 1)
tt_zone <- function(Z, g) { e <- Z[g, ]; l <- log2(e + 1)
  data.frame(gene = g, Mean_Control = mean(e[1:3]), Mean_EGM = mean(e[4:6]),
             P_raw = t.test(e[4:6], e[1:3])$p.value, P_log2 = t.test(l[4:6], l[1:3])$p.value) }
tt_portal  <- bind_rows(lapply(portal_v,  function(g) tt_zone(Z_peri, g))) %>%
  left_join(liver_deg$Periportal  %>% select(gene, log2FoldChange, padj), by = "gene")
# ---- pericentral genes ----
central_v <- c("Igfbp4", "Mmp14", "Mki67", "Nr5a2")
tt_central <- bind_rows(lapply(central_v, function(g) tt_zone(Z_cent, g))) %>%
  left_join(liver_deg$Pericentral %>% select(gene, log2FoldChange, padj), by = "gene")
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
# ---- p label, row builder, 2-row layout ----
fmt_p <- function(p) if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)
make_row_z <- function(Z, genes, pvals, strip = FALSE) lapply(seq_along(genes), function(i) {
  p <- plot_gene_box_clean(Z, genes[i], liver_meta, pval = pvals[i], ylab = "Expression (log2 counts + 1)")
  if (strip) p <- strip_text(p); if (i != 1) p <- p + theme(axis.title.y = element_blank()); p })
build_liver_box2 <- function(strip = FALSE) {
  tcol <- if (strip) "transparent" else "black"
  wrap_plots(c(list(row_title("Periportal hepatocytes", tcol)),  make_row_z(Z_peri_log, portal_v,  tt_portal$padj,  strip),
               list(row_title("Pericentral hepatocytes", tcol)), make_row_z(Z_cent_log, central_v, tt_central$padj, strip)),
             design = design, heights = c(0.1, 1, 0.1, 1))
}
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
# ---- build the plots ----
p_liver2 <- build_liver_box2(FALSE); p_liver2_nt <- build_liver_box2(TRUE)
# ---- output folder ----
liver_final <- file.path(root, "integrated_analysis_for_lizard_factors/integrated_analysis/liver/final_aging")
dir.create(liver_final, showWarnings = FALSE, recursive = TRUE)
# ---- save ----
ggsave(file.path(liver_final, "liver_aging_boxplot_FINAL.pdf"),        p_liver2,    width = 12, height = 8, dpi = 300)
ggsave(file.path(liver_final, "liver_aging_boxplot_FINAL.png"),        p_liver2,    width = 12, height = 8, dpi = 300)
ggsave(file.path(liver_final, "liver_aging_boxplot_FINAL_notext.pdf"), p_liver2_nt, width = 12, height = 8, dpi = 300)
liver_box_stats <- bind_rows(tt_portal %>% mutate(Zone = "Periportal"), tt_central %>% mutate(Zone = "Pericentral")) %>%
  select(Zone, gene, Mean_Control, Mean_EGM, log2FoldChange, padj_DESeq2 = padj, P_ttest_raw = P_raw, P_ttest_log2 = P_log2)
write.csv(liver_box_stats, file.path(liver_final, "liver_aging_boxplot_genes_stats.csv"), row.names = FALSE)
