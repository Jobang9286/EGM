# ============================================================
# Kidney aging-feature marker screening (HGPS_EGM vs HGPS_control)
# - regeneration marker screening(kidney_deg.R / violin_plot_script.R)과 동일한 형식,
#   gene set만 aging 버전(Senescence, SASP, Inflammation, Fibrosis, Injury, ...)으로 교체
# - 입력: cell type별 DESeq2 결과(deg_results/cell_type_specific/*.csv),
#         Z_PT_full / Z_Endo_full / sample_metadata (deg_results_2/kidney_deconvolution_main_objects.RData)
# - 출력: 이 폴더(aging_screening/)의 csv / pdf / png
# Date: 2026-09-06
# ============================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ComplexHeatmap); library(circlize); library(ggplot2); library(patchwork)})
root    <- "EGM_project_for_deconvolution"
deg_dir <- file.path(root, "integrated_analysis_for_lizard_factors/kidney/deg_results/cell_type_specific")
out_dir <- file.path(root, "integrated_analysis_for_lizard_factors/kidney/aging_screening")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. 데이터 ----
cell_types <- c("PT","Endo","LOH","DCT","CD-Trans","Podo","CD-PC","Macro")
deg_list <- setNames(lapply(cell_types, function(ct)
  read.csv(file.path(deg_dir, paste0(ct, "_DEG_HGPS_EGM_vs_control.csv")))), cell_types)
load(file.path(root, "integrated_analysis_for_lizard_factors/kidney/deg_results_2/kidney_deconvolution_main_objects.RData"))
hg <- sample_metadata$Sample   # HGPS_control_1..3, HGPS_EGM_1..3

# ---- 2. Aging gene sets ----
aging_sets <- list(
  Senescence = c("Cdkn1a","Cdkn2a","Cdkn2b","Cdkn1b","Cdkn1c","Trp53","Glb1","Serpine1","Lmnb1","Hmgb1","Mdm2","Rb1"),
  SASP = c("Il6","Il1a","Il1b","Il18","Il1rn","Cxcl1","Cxcl2","Cxcl10","Cxcl12","Cxcl16","Ccl2","Ccl3","Ccl4","Ccl5","Ccl7","Ccl8","Ccl20",
           "Mmp3","Mmp9","Mmp12","Mmp13","Mmp14","Timp1","Igfbp2","Igfbp3","Igfbp4","Igfbp5","Igfbp7","Tnf","Vegfa","Areg","Ereg",
           "Csf1","Csf2","Fgf2","Hgf","Plaur","Lif","Osm","Il15","Icam1","Vcam1","Gdf15"),
  Inflammation = c("Nfkb1","Nfkb2","Rela","Relb","Ikbkb","Tlr2","Tlr4","Myd88","Stat1","Stat3","Irf1","Irf7","C3","C1qa","C1qb","C1qc","C4b",
                   "B2m","H2-K1","H2-D1","Cd74","Lgals3","Cd68","Adgre1","Ptprc","S100a8","S100a9","Nlrp3","Casp1","Cd44","Ifit1","Ifit3",
                   "Isg15","Oasl1","Cxcl9","Ccr2","Itgam","Trem2"),
  Fibrosis = c("Col1a1","Col1a2","Col3a1","Col4a1","Col4a2","Col5a1","Col5a2","Col6a1","Col6a3","Fn1","Acta2","Tagln","Tgfb1","Tgfb2","Tgfb3",
               "Ccn2","Ctgf","Postn","Vim","Lox","Loxl2","Sparc","Mmp2","Thbs1","Thbs2","Eln","Fbn1","Dcn","Lum","Pdgfrb","Pdgfb","Des","Fap",
               "Timp2","Serpinh1","Bgn","Tnc"),
  Injury_Damage = c("Havcr1","Lcn2","Spp1","Klf6","Atf3","Sox9","Krt20","Clu","Tnfrsf12a","Cst3","Gpnmb","Ccn1","Cyr61","Egr1","Fos","Jun","Junb",
                    "Lgals1","Anxa2","Anxa3","Cp","Tnfaip3","Nfkbia","Ier3"),
  Longevity_Nutrient = c("Kl","Klb","Sirt1","Sirt6","Sirt7","Foxo1","Foxo3","Foxo4","Mtor","Rptor","Rictor","Prkaa1","Prkaa2","Igf1","Igf1r","Insr",
                         "Irs1","Irs2","Fgf23","Fgf21","Nampt","Nmnat1","Cd38","Akt1","Tsc1","Tsc2","Ppard","Adipor1"),
  Tubular_Function = c("Umod","Egf","Lrp2","Cubn","Slc34a1","Slc5a2","Slc22a6","Slc22a8","Aqp1","Aqp2","Slc12a1","Slc12a3","Nphs1","Nphs2","Podxl",
                       "Wt1","Cdh16","Hnf1b","Hnf4a","Pax2","Pax8","Ggt1","Kap","Calb1","Pvalb","Fxyd2","Atp1a1"),
  DNA_Damage = c("Atm","Atr","Chek1","Chek2","H2ax","H2afx","Gadd45a","Gadd45b","Gadd45g","Parp1","Brca1","Brca2","Rad51","Xrcc1","Xrcc5","Xrcc6",
                 "Trp53bp1","Mre11a","Nbn","Rad50","Ddb2","Ercc1","Ercc2","Ogg1","Apex1","Lig4","Prkdc","Mgmt","Rpa1","Rad23b","Xpa","Xpc"),
  Oxidative_Stress = c("Sod1","Sod2","Sod3","Cat","Gpx1","Gpx3","Gpx4","Nfe2l2","Keap1","Hmox1","Nqo1","Txn1","Txnrd1","Prdx1","Prdx3","Prdx6",
                       "Gclc","Gclm","Gsr","Gstm1","Gsta4","Gstp1","Nox4","Cybb","Ncf1","Ncf2","Srxn1","Mt1","Mt2"),
  Mitochondrial = c("Ppargc1a","Ppargc1b","Tfam","Nrf1","Sirt3","Sirt4","Sirt5","Mfn1","Mfn2","Opa1","Dnm1l","Fis1","Cox4i1","Cox5a","Ndufa4","Ndufs1",
                    "Ndufb8","Sdha","Sdhb","Uqcrc1","Uqcrc2","Atp5a1","Atp5f1a","Cycs","Pink1","Prkn","Ucp2","Cpt1a","Acadm","Ppara","Esrra","Polg","Tomm20","Vdac1"),
  Autophagy_Proteostasis = c("Atg5","Atg7","Atg12","Atg16l1","Becn1","Map1lc3b","Map1lc3a","Sqstm1","Ulk1","Ulk2","Lamp1","Lamp2","Ctsb","Ctsd","Ctsl","Tfeb",
                             "Hspa1a","Hspa1b","Hspa5","Hspa8","Hsp90aa1","Hspb1","Hspd1","Dnajb1","Hsf1","Bnip3","Psmb5","Psmd1","Psma1","Ubb","Ubc",
                             "Atf4","Ddit3","Xbp1","Ern1","Eif2ak3","Nbr1","Wipi2","Gabarap"),
  Epigenetic = c("Dnmt1","Dnmt3a","Dnmt3b","Tet1","Tet2","Tet3","Hdac1","Hdac2","Hdac3","Hdac4","Kdm6a","Kdm6b","Kdm5b","Kdm4a","Ezh2","Kat2a","Kat6a",
                 "Setd7","Suv39h1","Cbx5","H3f3b","Brd4","Kmt2a","Ep300","Crebbp","Mecp2"),
  Nuclear_Lamina = c("Lmna","Lmnb2","Emd","Lbr","Sun1","Sun2","Zmpste24","Banf1","Tmpo","Syne1","Syne2","Nup153")
)
aging_markers <- bind_rows(lapply(names(aging_sets), function(n) data.frame(Category = n, Gene = aging_sets[[n]]))) %>%
  distinct(Gene, .keep_all = TRUE)
# 방향 주석: Pro-aging(노화 시 증가/유해) / Protective(항노화, 기능 유지) / Ambiguous
protective <- c("Lmnb1", aging_sets$Tubular_Function, aging_sets$Mitochondrial,
  "Sod1","Sod2","Sod3","Cat","Gpx1","Gpx3","Gpx4","Nfe2l2","Hmox1","Nqo1","Txn1","Txnrd1","Prdx1","Prdx3","Prdx6","Gclc","Gclm","Gsr","Gstm1","Gsta4","Gstp1","Srxn1","Mt1","Mt2",
  "Kl","Klb","Sirt1","Sirt6","Sirt7","Foxo1","Foxo3","Foxo4","Prkaa1","Prkaa2","Nampt","Nmnat1","Tsc1","Tsc2","Ppard","Adipor1","Fgf21",
  "Atg5","Atg7","Atg12","Atg16l1","Becn1","Map1lc3b","Map1lc3a","Ulk1","Ulk2","Lamp1","Lamp2","Ctsb","Ctsd","Ctsl","Tfeb","Hspa1a","Hspa1b","Hspa5","Hspa8",
  "Hsp90aa1","Hspb1","Hspd1","Dnajb1","Hsf1","Psmb5","Psmd1","Psma1","Ubb","Ubc","Nbr1","Wipi2","Gabarap",
  "Brca1","Brca2","Rad51","Xrcc1","Xrcc5","Xrcc6","Trp53bp1","Mre11a","Nbn","Rad50","Ddb2","Ercc1","Ercc2","Ogg1","Apex1","Lig4","Prkdc","Mgmt","Rpa1","Rad23b","Xpa","Xpc",
  "Lmnb2","Emd","Lbr","Sun1","Sun2","Zmpste24","Banf1","Tmpo","Syne1","Syne2","Nup153")
ambiguous <- c(aging_sets$Epigenetic, "Atm","Atr","Chek1","Chek2","H2ax","H2afx","Gadd45a","Gadd45b","Gadd45g","Parp1","Lmna","Hmgb1","Trp53","Rb1","Mdm2","Cd38")
aging_markers$Aging_Role <- ifelse(aging_markers$Gene %in% ambiguous, "Ambiguous",
                             ifelse(aging_markers$Gene %in% protective, "Protective", "Pro-aging"))
write.csv(aging_markers, file.path(out_dir, "aging_gene_sets.csv"), row.names = FALSE)

# ---- 3. Screening (padj < 0.05) ----
aging_res <- bind_rows(lapply(cell_types, function(ct)
  deg_list[[ct]] %>% filter(gene %in% aging_markers$Gene) %>% mutate(CellType = ct))) %>%
  left_join(aging_markers, by = c("gene" = "Gene")) %>%
  mutate(Effect = case_when(
    is.na(padj) | padj >= 0.05 ~ "NS",
    Aging_Role == "Pro-aging"  & log2FoldChange < 0 ~ "Improved",   # 노화 마커 감소
    Aging_Role == "Protective" & log2FoldChange > 0 ~ "Improved",   # 보호 인자 증가
    Aging_Role == "Ambiguous" ~ "Ambiguous",
    TRUE ~ "Worsened"))
write.csv(aging_res, file.path(out_dir, "aging_markers_all_celltypes.csv"), row.names = FALSE)
sig <- aging_res %>% filter(padj < 0.05)
write.csv(sig %>% arrange(Category, gene, CellType), file.path(out_dir, "aging_markers_significant_padj0.05.csv"), row.names = FALSE)
gene_tbl <- sig %>% group_by(Category, gene, Aging_Role) %>%
  summarise(nCT = n(), meanFC = round(mean(log2FoldChange), 2), minPadj = signif(min(padj), 2),
            CTs = paste0(CellType, "(", ifelse(log2FoldChange > 0, "+", "-"), ")", collapse = " "), .groups = "drop") %>%
  arrange(Category, desc(nCT), minPadj)
write.csv(gene_tbl, file.path(out_dir, "aging_markers_significant_by_gene.csv"), row.names = FALSE)

# ---- 4. Heatmap (kidney_deg.R 스타일) ----
cat_order <- c("Senescence","SASP","Inflammation","Fibrosis","Injury_Damage","Longevity_Nutrient","Oxidative_Stress",
               "Mitochondrial","Autophagy_Proteostasis","DNA_Damage","Epigenetic","Nuclear_Lamina","Tubular_Function")
cat_colors <- c(Senescence = "#3C5488FF", SASP = "#E64B35FF", Inflammation = "#F39B7FFF", Fibrosis = "#4DBBD5FF",
                Injury_Damage = "#91A54667", Longevity_Nutrient = "#00A087FF", Oxidative_Stress = "#FFB300FF",
                Mitochondrial = "#7E6148FF", Autophagy_Proteostasis = "#8491B4FF", DNA_Damage = "#B09C85FF",
                Epigenetic = "#DC0000FF", Nuclear_Lamina = "#631879FF", Tubular_Function = "#66C2A5FF")
build_mats <- function(genes) {
  d <- sig %>% filter(gene %in% genes)
  fc <- d %>% select(gene, CellType, log2FoldChange) %>% pivot_wider(names_from = CellType, values_from = log2FoldChange) %>% tibble::column_to_rownames("gene") %>% as.matrix()
  pv <- d %>% select(gene, CellType, padj) %>% pivot_wider(names_from = CellType, values_from = padj) %>% tibble::column_to_rownames("gene") %>% as.matrix()
  fc[is.na(fc)] <- 0
  cts <- cell_types[cell_types %in% colnames(fc)]
  list(fc = fc[genes, cts, drop = FALSE], pv = pv[genes, cts, drop = FALSE])
}
make_ht <- function(genes, title, height_cm, fontsize = 10) {
  m <- build_mats(genes)
  gcat <- aging_markers$Category[match(genes, aging_markers$Gene)]
  add_sig <- function(j, i, x, y, width, height, fill) {
    p <- m$pv[i, j]
    if (!is.na(p)) { s <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
      grid.text(s, x, y, gp = gpar(fontsize = 8, col = "black", fontface = "bold")) }
  }
  row_ha <- rowAnnotation(Category = factor(gcat, levels = cat_order), col = list(Category = cat_colors),
                          annotation_name_gp = gpar(fontsize = 10, fontface = "bold"), width = unit(0.5, "cm"))
  Heatmap(m$fc, name = "log2FC", col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
          cluster_rows = FALSE, cluster_columns = TRUE, show_row_names = TRUE, show_column_names = TRUE,
          row_names_gp = gpar(fontsize = fontsize), column_names_gp = gpar(fontsize = 9),
          left_annotation = row_ha, rect_gp = gpar(col = "grey90", lwd = 0.3), cell_fun = add_sig,
          column_title = title, column_title_gp = gpar(fontsize = 12, fontface = "bold"),
          heatmap_legend_param = list(title = "log2FC", title_gp = gpar(fontsize = 10, fontface = "bold"),
                                      at = c(-3, -2, -1, 0, 1, 2, 3), legend_height = unit(4, "cm")),
          width = unit(6, "cm"), height = unit(height_cm, "cm"), border = TRUE, border_gp = gpar(col = "black", lwd = 1))
}
save_ht <- function(ht, base, w, h) {
  pdf(file.path(out_dir, paste0(base, ".pdf")), width = w, height = h); draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right"); dev.off()
  png(file.path(out_dir, paste0(base, ".png")), width = w, height = h, units = "in", res = 200); draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right"); dev.off()
}
title_txt <- "Kidney\n(HGPS_EGM vs HGPS control)"

# 4a. 전체 유의 유전자 (supplementary)
genes_full <- gene_tbl %>% mutate(Category = factor(Category, levels = cat_order)) %>% arrange(Category, desc(nCT), minPadj) %>% pull(gene)
pdf(file.path(out_dir, "kidney_aging_heatmap_ALL_significant.pdf"), width = 10, height = 0.14 * length(genes_full) + 4)
draw(make_ht(genes_full, title_txt, 0.32 * length(genes_full), fontsize = 7), heatmap_legend_side = "right", annotation_legend_side = "right"); dev.off()

# 4b. 카테고리별 top-3 (nCT, padj) 스크리닝 버전
sel <- gene_tbl %>% mutate(Category = factor(Category, levels = cat_order)) %>% arrange(Category, desc(nCT), minPadj) %>%
  group_by(Category) %>% slice_head(n = 3) %>% ungroup()
write.csv(sel, file.path(out_dir, "kidney_aging_heatmap_top3_per_category_genes.csv"), row.names = FALSE)
save_ht(make_ht(sel$gene, title_txt, 0.55 * nrow(sel)), "kidney_aging_heatmap_top3_per_category", 10, 0.28 * nrow(sel) + 4)

# 4c. 핵심 aging 카테고리 27-gene 버전 (재생 히트맵 크기에 맞춤) - 수동 큐레이션 대상
core_genes <- c("Cdkn1a","Cdkn1b","Cdkn2b",            # Senescence
                "Igfbp7","Vegfa","Igfbp5","Cxcl16",    # SASP
                "Lgals3","Irf1","Relb","C3",           # Inflammation
                "Serpinh1","Thbs1","Col5a2","Acta2",   # Fibrosis
                "Spp1","Klf6","Junb","Egr1",           # Injury_Damage
                "Kl","Foxo3","Nampt","Igf1",           # Longevity_Nutrient
                "Hnf1b","Hnf4a","Ggt1","Egf")          # Tubular_Function
core_genes <- core_genes[core_genes %in% gene_tbl$gene]
save_ht(make_ht(core_genes, title_txt, 0.55 * length(core_genes)), "kidney_aging_heatmap_core_27genes", 10, 12)
write.csv(gene_tbl %>% filter(gene %in% core_genes) %>% arrange(match(gene, core_genes)),
          file.path(out_dir, "kidney_aging_heatmap_core_27genes_table.csv"), row.names = FALSE)

# ---- 5. Violin plots (violin_plot_script.R 스타일) + t-test ----
fmt_p <- function(p) if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)
plot_gene_violin_clean <- function(Z_matrix, gene, sample_metadata, pval = NULL, ylab = "Expression (log2 TPM + 1)") {
  if (!gene %in% rownames(Z_matrix)) return(NULL)
  gene_expr <- Z_matrix[gene, ]
  plot_df <- data.frame(Sample = names(gene_expr), Expression = gene_expr,
                        Group = sample_metadata$Group[match(names(gene_expr), sample_metadata$Sample)]) %>%
    filter(Group %in% c("HGPS_control", "HGPS_EGM")) %>% mutate(Group = factor(Group, levels = c("HGPS_control", "HGPS_EGM")))
  ggplot(plot_df, aes(x = Group, y = Expression, fill = Group)) +
    geom_violin(alpha = 0.7, trim = FALSE) + geom_jitter(width = 0.1, size = 2.5, alpha = 0.7) +
    scale_fill_manual(values = c("HGPS_control" = "#E64B35", "HGPS_EGM" = "#4DBBD5")) +
    scale_x_discrete(labels = c("Control", "EGM")) +
    labs(title = gene, x = "", y = ylab, caption = if (!is.null(pval)) fmt_p(pval) else NULL) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 16), plot.caption = element_text(hjust = 0.5, size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 11),
          axis.title.y = element_text(size = 12), legend.position = "none",
          panel.grid.major = element_line(color = "grey90"), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", linewidth = 1))
}
ttest_both <- function(Z, g) {   # raw Z count 및 log2(Z+1) 두 스케일의 Welch t-test
  e <- Z[g, hg]; l <- log2(e + 1)
  data.frame(gene = g, Mean_Control = mean(e[1:3]), Mean_EGM = mean(e[4:6]),
             P_raw = t.test(e[4:6], e[1:3])$p.value, P_log2 = t.test(l[4:6], l[1:3])$p.value)
}
cand_all <- function(ct, Z) {
  d <- sig %>% filter(CellType == ct) %>% select(gene, Category, Aging_Role, Effect, log2FoldChange, padj)
  d %>% left_join(bind_rows(lapply(d$gene, function(g) ttest_both(Z, g))), by = "gene") %>% arrange(P_raw)
}
pt_all <- cand_all("PT", Z_PT_full); endo_all <- cand_all("Endo", Z_Endo_full)
write.csv(pt_all,   file.path(out_dir, "ProximalTubule_all_sig_aging_genes_ttest.csv"), row.names = FALSE)
write.csv(endo_all, file.path(out_dir, "Endothelial_all_sig_aging_genes_ttest.csv"), row.names = FALSE)

mk_panel2 <- function(Z, genes, pvals, ncol, file, w, h, title, ylab) {
  ps <- lapply(seq_along(genes), function(i) plot_gene_violin_clean(Z, genes[i], sample_metadata, pval = pvals[i], ylab = ylab))
  p <- wrap_plots(ps[!sapply(ps, is.null)], ncol = ncol) +
    plot_annotation(title = title, theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
  ggsave(file.path(out_dir, file), p, width = w, height = h, dpi = 300)
  ggsave(file.path(out_dir, sub("\\.pdf$", ".png", file)), p, width = w, height = h, dpi = 150)
}
# 선택 규칙: Improved 방향, 핵심 카테고리 우선, t-test p 낮은 순 (PT 6개, Endo 3개)
core_cats <- c("Senescence","SASP","Inflammation","Fibrosis","Injury_Damage","Longevity_Nutrient","Tubular_Function")
pick2 <- function(tab, n, pcol) { x <- tab %>% filter(Effect == "Improved", .data[[pcol]] < 0.05) %>%
  mutate(core = Category %in% core_cats) %>% arrange(desc(core), .data[[pcol]]); head(x$gene, n) }
gp <- function(tab, g, col) tab[[col]][match(g, tab$gene)]
pt_raw <- pick2(pt_all, 6, "P_raw"); endo_raw <- pick2(endo_all, 3, "P_raw")
pt_log <- pick2(pt_all, 6, "P_log2"); endo_log <- pick2(endo_all, 3, "P_log2")
# 5a. raw Z 스케일 (기존 regeneration violin과 동일 형식)
mk_panel2(Z_PT_full,   pt_raw,   gp(pt_all, pt_raw, "P_raw"),     3, "ProximalTubule_aging_violin_plots.pdf", 12, 8.5, "Aging markers in Proximal tubular cells", "Expression (log2 TPM + 1)")
mk_panel2(Z_Endo_full, endo_raw, gp(endo_all, endo_raw, "P_raw"), 3, "Endothelial_aging_violin_plots.pdf",    12, 4.5, "Aging markers in Endothelial cells",      "Expression (log2 TPM + 1)")
# 5b. log2(Z+1) 스케일 (control_1 편차 영향이 줄어 하향 유전자도 유의)
Z_PT_log <- log2(Z_PT_full + 1); Z_Endo_log <- log2(Z_Endo_full + 1)
mk_panel2(Z_PT_log,   pt_log,   gp(pt_all, pt_log, "P_log2"),     3, "ProximalTubule_aging_violin_plots_log2.pdf", 12, 8.5, "Aging markers in Proximal tubular cells", "Expression (log2 counts + 1)")
mk_panel2(Z_Endo_log, endo_log, gp(endo_all, endo_log, "P_log2"), 3, "Endothelial_aging_violin_plots_log2.pdf",    12, 4.5, "Aging markers in Endothelial cells",      "Expression (log2 counts + 1)")
cat("done. outputs in", out_dir, "\n")
