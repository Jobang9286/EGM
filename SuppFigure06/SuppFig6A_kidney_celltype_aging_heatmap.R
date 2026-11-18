# ============================================================
# Kidney aging-feature heatmap - FINAL (26 genes, 9 categories, 2 cell-type groups)
# HGPS_EGM vs HGPS_control, cell type-specific DESeq2 on InstaPrism Z counts
# 입력: integrated_analysis_for_lizard_factors/kidney/deg_results/cell_type_specific/*_DEG_HGPS_EGM_vs_control.csv
# 출력: 이 폴더 (final_aging/) - 텍스트 있는 버전 / 없는 버전 PDF
# Date: 2026-09-06
# ============================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ComplexHeatmap); library(circlize)})
root      <- "EGM_project_for_deconvolution"
deg_dir   <- file.path(root, "integrated_analysis_for_lizard_factors/kidney/deg_results/cell_type_specific")
final_dir <- file.path(root, "integrated_analysis_for_lizard_factors/integrated_analysis/kidney/final_aging")
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 빠른 재현: 저장된 워크스페이스 ----
# load(file.path(final_dir, "kidney_aging_final_workspace.RData"))
# draw(ht_final2, column_title = "Kidney (HGPS_EGM vs HGPS control)", heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)

# ---- 1. 최종 유전자 / 카테고리 ----
final_cur <- tibble::tribble(
  ~gene, ~Fig,
  # aging features: 감소
  "Cdkn1a","Senescence", "Cdkn1b","Senescence", "Cdkn2b","Senescence",
  "Igfbp7","SASP", "Vegfa","SASP", "Igfbp5","SASP",
  "C3","Inflammation",
  "Thbs1","Fibrosis", "Col5a2","Fibrosis",
  "Spp1","Tubular_Injury", "Klf6","Tubular_Injury", "Cp","Tubular_Injury",
  "Atf4","Cell_Stress", "Nox4","Cell_Stress",
  # repair / function: 증가
  "Ccnb1","Proliferation_Repair", "Cdk1","Proliferation_Repair", "Ccne1","Proliferation_Repair", "Mcm2","Proliferation_Repair", "Igf1","Proliferation_Repair",
  "Hnf1b","PT_Function", "Hnf4a","PT_Function", "Slc22a6","PT_Function", "Slc22a12","PT_Function", "Ggt1","PT_Function",
  "Aplnr","Endothelial", "Ehd3","Endothelial")
final_order  <- c("Senescence","SASP","Inflammation","Fibrosis","Tubular_Injury","Cell_Stress","Proliferation_Repair","PT_Function","Endothelial")
final_colors <- c(Senescence = "#3C5488", SASP = "#E64B35", Inflammation = "#F39B7F", Fibrosis = "#4DBBD5", Tubular_Injury = "#91A546",
                  Cell_Stress = "#8491B4", Proliferation_Repair = "#DC0000", PT_Function = "#00A087", Endothelial = "#FFB300")
cell_types <- c("PT","Endo","LOH","DCT","CD-Trans","Podo","CD-PC","Macro")
col_order  <- c("PT","LOH","DCT","CD-Trans","CD-PC","Endo","Podo","Macro")
col_groups <- setNames(c(rep("Tubular epithelium", 5), rep("Vascular / Glomerular / Immune", 3)), col_order)
grp_levels <- c("Tubular epithelium", "Vascular / Glomerular / Immune")
grp_colors <- c("Tubular epithelium" = "#4DBBD5", "Vascular / Glomerular / Immune" = "#F39B7F")

# ---- 2. DESeq2 결과 -> log2FC / padj 행렬 (padj<0.05 셀만, 나머지 0) ----
deg_list <- setNames(lapply(cell_types, function(ct) read.csv(file.path(deg_dir, paste0(ct, "_DEG_HGPS_EGM_vs_control.csv")))), cell_types)
all_final <- bind_rows(lapply(cell_types, function(ct) deg_list[[ct]] %>% filter(gene %in% final_cur$gene) %>% mutate(CellType = ct)))
final_sig <- all_final %>% filter(padj < 0.05)
fc <- final_sig %>% select(gene, CellType, log2FoldChange) %>% pivot_wider(names_from = CellType, values_from = log2FoldChange) %>% tibble::column_to_rownames("gene") %>% as.matrix()
pv <- final_sig %>% select(gene, CellType, padj) %>% pivot_wider(names_from = CellType, values_from = padj) %>% tibble::column_to_rownames("gene") %>% as.matrix()
fc[is.na(fc)] <- 0
fc2 <- fc[final_cur$gene, col_order]; pv2 <- pv[final_cur$gene, col_order]

# ---- 3. Heatmap ----
star_fun <- function(j, i, x, y, w, h, fill) { p <- pv2[i, j]
  if (!is.na(p)) grid.text(if (p < 0.001) "***" else if (p < 0.01) "**" else "*", x, y, gp = gpar(fontsize = 8, col = "black", fontface = "bold")) }
col_fun <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
col_split <- factor(col_groups[col_order], levels = grp_levels)

ht_final2 <- Heatmap(fc2, name = "log2FC", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = col_split, column_gap = unit(2.5, "mm"),
  column_title = c("Tubular\nepithelium", "Vascular / Glomerular\n/ Immune"), column_title_gp = gpar(fontsize = 9.5, fontface = "bold"),
  top_annotation = HeatmapAnnotation(`Cell group` = col_split, col = list(`Cell group` = grp_colors), show_annotation_name = FALSE,
                                     annotation_legend_param = list(`Cell group` = list(title_gp = gpar(fontsize = 10, fontface = "bold"))),
                                     simple_anno_size = unit(3, "mm")),
  show_row_names = TRUE, show_column_names = TRUE, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 9),
  left_annotation = rowAnnotation(Category = factor(final_cur$Fig, levels = final_order), col = list(Category = final_colors),
                                  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"), width = unit(0.5, "cm")),
  rect_gp = gpar(col = "grey90", lwd = 0.3), cell_fun = star_fun,
  heatmap_legend_param = list(title = "log2FC", title_gp = gpar(fontsize = 10, fontface = "bold"), at = c(-3, -2, -1, 0, 1, 2, 3), legend_height = unit(4, "cm")),
  width = unit(6.4, "cm"), height = unit(0.5 * nrow(final_cur), "cm"), border = TRUE, border_gp = gpar(col = "black", lwd = 1))

# 텍스트 없는 버전: 유전자명 / cell type명 / 제목 / 범례 제거 (별표, 카테고리·그룹 색 띠만 유지)
ht_notext <- Heatmap(fc2, name = "log2FC", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = col_split, column_gap = unit(2.5, "mm"), column_title = NULL,
  top_annotation = HeatmapAnnotation(`Cell group` = col_split, col = list(`Cell group` = grp_colors),
                                     show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(3, "mm")),
  show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE,
  left_annotation = rowAnnotation(Category = factor(final_cur$Fig, levels = final_order), col = list(Category = final_colors),
                                  show_annotation_name = FALSE, show_legend = FALSE, width = unit(0.5, "cm")),
  rect_gp = gpar(col = "grey90", lwd = 0.3), cell_fun = star_fun,
  width = unit(6.4, "cm"), height = unit(0.5 * nrow(final_cur), "cm"), border = TRUE, border_gp = gpar(col = "black", lwd = 1))

# 오른쪽 범례만 별도 객체 (텍스트 있는 버전 / 없는 버전)
make_legends <- function(with_text = TRUE) {
  lab_cat <- if (with_text) final_order else rep("", length(final_order))
  lab_grp <- if (with_text) grp_levels else rep("", length(grp_levels))
  ttl <- function(x) if (with_text) x else ""
  packLegend(
    Legend(labels = lab_cat, legend_gp = gpar(fill = final_colors[final_order]), title = ttl("Category"),
           title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 10), grid_height = unit(4, "mm"), grid_width = unit(4, "mm")),
    Legend(labels = lab_grp, legend_gp = gpar(fill = grp_colors[grp_levels]), title = ttl("Cell group"),
           title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 10), grid_height = unit(4, "mm"), grid_width = unit(4, "mm")),
    Legend(col_fun = col_fun, title = ttl("log2FC"), at = c(-3, -2, -1, 0, 1, 2, 3),
           labels = if (with_text) c(-3, -2, -1, 0, 1, 2, 3) else rep("", 7),
           title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 10), legend_height = unit(4, "cm")),
    direction = "vertical", gap = unit(5, "mm"))
}
lg_text <- make_legends(TRUE); lg_notext <- make_legends(FALSE)

# ---- 4. 저장 ----
pdf(file.path(final_dir, "kidney_aging_heatmap_26genes_FINAL.pdf"), width = 9, height = 10)
draw(ht_final2, column_title = "Kidney (HGPS_EGM vs HGPS control)", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
     heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE); dev.off()
png(file.path(final_dir, "kidney_aging_heatmap_26genes_FINAL.png"), width = 9, height = 10, units = "in", res = 300)
draw(ht_final2, column_title = "Kidney (HGPS_EGM vs HGPS control)", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
     heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE); dev.off()
pdf(file.path(final_dir, "kidney_aging_heatmap_26genes_FINAL_notext.pdf"), width = 5, height = 9); draw(ht_notext); dev.off()
pdf(file.path(final_dir, "kidney_aging_heatmap_26genes_FINAL_legend.pdf"), width = 3.2, height = 6.5)
grid.newpage(); draw(lg_text, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "centre"); dev.off()
pdf(file.path(final_dir, "kidney_aging_heatmap_26genes_FINAL_legend_notext.pdf"), width = 2, height = 6.5)
grid.newpage(); draw(lg_notext, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "centre"); dev.off()
cat("saved to", final_dir, "\n")
