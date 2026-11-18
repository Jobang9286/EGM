source("Fig2P-Q_00_aging_gene_sets_setup.R")   # root, aging_sets, protective, ambiguous

## input loading.
brain_dir <- file.path(root, "integrated_analysis_for_lizard_factors/brain")
bd <- readRDS(file.path(brain_dir, "DEG_results_InstaPrism/Brain_DEG_Complete_Results.rds"))

## brain aging screen gene sets.
brain_ct <- c("Excitatory Neurons","Inhibitory Neurons","Astrocytes","Microglia","ODC")
# 뇌 특이 기능/재생 마커 (kidney의 Tubular_Function, liver의 Hepatocyte_Function 자리)
brain_sets <- c(aging_sets[setdiff(names(aging_sets), "Tubular_Function")], list(
  Neuronal_Function = c("Syn1","Syp","Snap25","Dlg4","Grin1","Grin2a","Grin2b","Gria1","Gria2","Camk2a","Camk2b","Bdnf","Ntrk2","Nefl","Nefm","Nefh",
                        "Map2","Tubb3","Rbfox3","Sv2a","Vamp2","Stx1a","Syt1","Nrgn","Calm1","Arc","Homer1","Shank3","Nlgn1","Nrxn1","Pvalb","Sst","Vip","Gad1","Gad2","Slc17a7","Slc32a1"),
  Glia_Myelin = c("Mbp","Plp1","Mog","Mag","Cnp","Mal","Sox10","Olig1","Olig2","Cldn11","Mobp","Ugt8a","Gfap","Aqp4","Slc1a2","Slc1a3","Gja1","Aldh1l1","S100b","Kcnj10","Glul","Vim"),
  Neuro_Regeneration = c("Gap43","Sox2","Sox4","Sox11","Dcx","Neurod1","Ascl1","Nes","Pax6","Tbr1","Eomes","Mki67","Pcna","Top2a","Mcm2","Ccnd1","Ccnb1","Cdk1","E2f1","Notch1","Hes5","Dll1",
                         "Igf1","Fgf2","Ngf","Ntf3","Vegfa","Cntf","Lif","Stat3","Creb1","Atf3","Klf6","Klf7","Sprr1a","Gal","Vgf")))
brain_markers <- bind_rows(lapply(names(brain_sets), function(n) data.frame(Category = n, Gene = brain_sets[[n]]))) %>% distinct(Gene, .keep_all = TRUE)
protective_brain <- c(protective, brain_sets$Neuronal_Function, brain_sets$Glia_Myelin)
brain_markers$Aging_Role <- ifelse(brain_markers$Gene %in% ambiguous, "Ambiguous",
                            ifelse(brain_markers$Gene %in% protective_brain, "Protective", "Pro-aging"))
brain_markers$Aging_Role[brain_markers$Category == "Neuro_Regeneration"] <- "Protective"
brain_out <- file.path(root, "integrated_analysis_for_lizard_factors/brain/aging_screening")
dir.create(brain_out, showWarnings = FALSE, recursive = TRUE)
write.csv(brain_markers, file.path(brain_out, "brain_aging_gene_sets.csv"), row.names = FALSE)
cat("gene set:", nrow(brain_markers), "genes /", length(brain_sets), "categories\n")

## per-cell-type hits for the aging gene sets.
brain_deg <- setNames(lapply(brain_ct, function(ct) bd$deg_results[[ct]]$full_results %>% rename(gene = Gene)), brain_ct)
brain_res <- bind_rows(lapply(brain_ct, function(ct) brain_deg[[ct]] %>% filter(gene %in% brain_markers$Gene) %>% mutate(CellType = ct))) %>%
  left_join(brain_markers, by = c("gene" = "Gene")) %>%
  mutate(Effect = case_when(is.na(padj) | padj >= 0.05 ~ "NS",
                            Aging_Role == "Pro-aging"  & log2FoldChange < 0 ~ "Improved",
                            Aging_Role == "Protective" & log2FoldChange > 0 ~ "Improved",
                            Aging_Role == "Ambiguous" ~ "Ambiguous", TRUE ~ "Worsened"))
write.csv(brain_res, file.path(brain_out, "brain_aging_markers_all_celltypes.csv"), row.names = FALSE)
bsig <- brain_res %>% filter(padj < 0.05)
cat("detected:", n_distinct(brain_res$gene), "/", nrow(brain_markers), "| significant genes:", n_distinct(bsig$gene), "| hits:", nrow(bsig), "\n\n")
print(as.data.frame(bsig %>% group_by(Category) %>% summarise(n_genes = n_distinct(gene), hits = n(), Up = sum(log2FoldChange > 0), Down = sum(log2FoldChange < 0),
      Improved = sum(Effect == "Improved"), Worsened = sum(Effect == "Worsened"), .groups = "drop") %>% arrange(desc(hits))), row.names = FALSE)
cat("\n--- by cell type ---\n")
print(as.data.frame(bsig %>% group_by(CellType) %>% summarise(n = n_distinct(gene), Up = sum(log2FoldChange > 0), Down = sum(log2FoldChange < 0),
      Improved = sum(Effect == "Improved"), Worsened = sum(Effect == "Worsened"), .groups = "drop")), row.names = FALSE)
bgene_tbl <- bsig %>% group_by(Category, gene, Aging_Role) %>%
  summarise(nCT = n(), meanFC = round(mean(log2FoldChange), 2), minPadj = signif(min(padj), 2), Effect = names(sort(table(Effect), decreasing = TRUE))[1],
            CTs = paste0(substr(CellType, 1, 4), "(", ifelse(log2FoldChange > 0, "+", "-"), ")", collapse = " "), .groups = "drop") %>%
  arrange(Category, desc(nCT), minPadj)
write.csv(bgene_tbl, file.path(brain_out, "brain_aging_markers_significant_by_gene.csv"), row.names = FALSE)
cat("\n--- significant genes by category ---\n"); print(as.data.frame(bgene_tbl), row.names = FALSE)

## first curated list; defines the column order (bcols) and the Neurons / Glia grouping.
brain_cur <- tibble::tribble(
  ~gene, ~Fig,
  # ---- 감소 ----
  "Cdkn1a","Senescence", "Cdkn1c","Senescence", "Serpine1","Senescence",
  "Igfbp4","SASP", "Igfbp3","SASP",
  "Cd68","Neuroinflammation", "Cd74","Neuroinflammation", "Apoe","Neuroinflammation",
  "Ccn2","ECM_Fibrosis", "Fn1","ECM_Fibrosis", "Col1a2","ECM_Fibrosis", "Col1a1","ECM_Fibrosis",
  "Junb","Injury_Damage", "Spp1","Injury_Damage", "Nfkbia","Injury_Damage",
  # ---- 증가 ----
  "Kl","Longevity",
  "Plp1","Myelination", "Opalin","Myelination", "Bcas1","Myelination", "Gjb1","Myelination", "Gpr37","Myelination",
  "Cx3cr1","Homeostatic_Microglia", "P2ry12","Homeostatic_Microglia", "Hexb","Homeostatic_Microglia", "Rgs10","Homeostatic_Microglia",
  "Gria1","Synaptic_Function", "Nrxn1","Synaptic_Function", "Kcnip4","Synaptic_Function")
brain_order <- c("Senescence","SASP","Neuroinflammation","ECM_Fibrosis","Injury_Damage","Longevity","Myelination","Homeostatic_Microglia","Synaptic_Function")
brain_colors <- c(Senescence = "#3C5488", SASP = "#E64B35", Neuroinflammation = "#F39B7F", ECM_Fibrosis = "#4DBBD5", Injury_Damage = "#91A546",
                  Longevity = "#B09C85", Myelination = "#00A087", Homeostatic_Microglia = "#FFB300", Synaptic_Function = "#631879")
ball <- bind_rows(lapply(brain_ct, function(ct) brain_deg[[ct]] %>% filter(gene %in% brain_cur$gene) %>% mutate(CellType = ct)))
bs <- ball %>% filter(padj < 0.05)
cat("유의 셀이 없는 유전자:", setdiff(brain_cur$gene, bs$gene), "\n")
bfc <- bs %>% select(gene, CellType, log2FoldChange) %>% pivot_wider(names_from = CellType, values_from = log2FoldChange) %>% tibble::column_to_rownames("gene") %>% as.matrix()
bpv <- bs %>% select(gene, CellType, padj) %>% pivot_wider(names_from = CellType, values_from = padj) %>% tibble::column_to_rownames("gene") %>% as.matrix()
bfc[is.na(bfc)] <- 0
bcols <- brain_ct[brain_ct %in% colnames(bfc)]
bfc <- bfc[brain_cur$gene, bcols]; bpv <- bpv[brain_cur$gene, bcols]
brain_grp <- setNames(c("Neurons","Neurons","Glia","Glia","Glia"), brain_ct)[bcols]
grp_lev_b <- c("Neurons", "Glia"); grp_col_b <- c(Neurons = "#4DBBD5", Glia = "#F39B7F")

## column split used by every later heatmap.
col_split_b <- factor(brain_grp, levels = grp_lev_b)

## merge Synaptic + Neurite, add Neuroprotection (Kl moved there).
brain_cur3 <- tibble::tribble(
  ~gene, ~Fig,
  "Cdkn1a","Senescence", "Cdkn1c","Senescence", "Serpine1","Senescence",
  "Igfbp4","SASP", "Igfbp3","SASP",
  "Cd68","Neuroinflammation", "Cd74","Neuroinflammation", "Apoe","Neuroinflammation",
  "Junb","Injury_Damage", "Spp1","Injury_Damage", "Nfkbia","Injury_Damage",
  "Plp1","Myelination", "Opalin","Myelination", "Bcas1","Myelination", "Gjb1","Myelination", "Gpr37","Myelination",
  "Cx3cr1","Homeostatic_Microglia", "P2ry12","Homeostatic_Microglia", "Hexb","Homeostatic_Microglia", "Rgs10","Homeostatic_Microglia",
  "Kl","Neuroprotection", "Penk","Neuroprotection", "Cryab","Neuroprotection", "Nptx2","Neuroprotection",   # Kl을 Myelination에서 이동
  "Gria1","Synaptic_Neurite", "Nrxn1","Synaptic_Neurite", "Kcnip4","Synaptic_Neurite",
  "Shtn1","Synaptic_Neurite", "Marcksl1","Synaptic_Neurite", "Nrp1","Synaptic_Neurite", "Plxna4","Synaptic_Neurite")
brain_order3 <- c("Senescence","SASP","Neuroinflammation","Injury_Damage","Myelination","Homeostatic_Microglia","Neuroprotection","Synaptic_Neurite")
brain_colors3 <- c(Senescence = "#3C5488", SASP = "#E64B35", Neuroinflammation = "#F39B7F", Injury_Damage = "#91A546",
                   Myelination = "#00A087", Homeostatic_Microglia = "#FFB300", Neuroprotection = "#B09C85", Synaptic_Neurite = "#631879")

## remove Myelination -> final 26 genes / 7 categories.
brain_cur4 <- brain_cur3 %>% filter(Fig != "Myelination")
brain_order4 <- setdiff(brain_order3, "Myelination")
brain_colors4 <- brain_colors3[brain_order4]
ball4 <- bind_rows(lapply(brain_ct, function(ct) brain_deg[[ct]] %>% filter(gene %in% brain_cur4$gene) %>% mutate(CellType = ct)))
bs4 <- ball4 %>% filter(padj < 0.05)
bfc4 <- bs4 %>% select(gene, CellType, log2FoldChange) %>% pivot_wider(names_from = CellType, values_from = log2FoldChange) %>% tibble::column_to_rownames("gene") %>% as.matrix()
bpv4 <- bs4 %>% select(gene, CellType, padj) %>% pivot_wider(names_from = CellType, values_from = padj) %>% tibble::column_to_rownames("gene") %>% as.matrix()
bfc4[is.na(bfc4)] <- 0; bfc4 <- bfc4[brain_cur4$gene, bcols]; bpv4 <- bpv4[brain_cur4$gene, bcols]
star_b4 <- function(j, i, x, y, w, h, fill) { p <- bpv4[i, j]
  if (!is.na(p)) grid.text(if (p < 0.001) "***" else if (p < 0.01) "**" else "*", x, y, gp = gpar(fontsize = 8, col = "black", fontface = "bold")) }

## final save.
brain_final <- file.path(root, "integrated_analysis_for_lizard_factors/integrated_analysis/brain/final_aging")
dir.create(brain_final, showWarnings = FALSE, recursive = TRUE)
col_fun_b <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
# 히트맵 (텍스트 있음, 범례 제외)
ht_b <- Heatmap(bfc4, name = "log2FC", col = col_fun_b, cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = col_split_b, column_gap = unit(2.5, "mm"),
  column_title = c("Neurons", "Glia"), column_title_gp = gpar(fontsize = 9.5, fontface = "bold"),
  top_annotation = HeatmapAnnotation(`Cell group` = col_split_b, col = list(`Cell group` = grp_col_b),
                                     show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(3, "mm")),
  show_row_names = TRUE, show_column_names = TRUE, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 9),
  left_annotation = rowAnnotation(Category = factor(brain_cur4$Fig, levels = brain_order4), col = list(Category = brain_colors4),
                                  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"), show_legend = FALSE, width = unit(0.5, "cm")),
  rect_gp = gpar(col = "grey90", lwd = 0.3), cell_fun = star_b4, show_heatmap_legend = FALSE,
  width = unit(5, "cm"), height = unit(0.5 * nrow(brain_cur4), "cm"), border = TRUE, border_gp = gpar(col = "black", lwd = 1))
# 히트맵 (텍스트 없음)
ht_b_nt <- Heatmap(bfc4, name = "log2FC", col = col_fun_b, cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = col_split_b, column_gap = unit(2.5, "mm"), column_title = NULL,
  top_annotation = HeatmapAnnotation(`Cell group` = col_split_b, col = list(`Cell group` = grp_col_b),
                                     show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(3, "mm")),
  show_row_names = FALSE, show_column_names = FALSE,
  left_annotation = rowAnnotation(Category = factor(brain_cur4$Fig, levels = brain_order4), col = list(Category = brain_colors4),
                                  show_annotation_name = FALSE, show_legend = FALSE, width = unit(0.5, "cm")),
  rect_gp = gpar(col = "grey90", lwd = 0.3), cell_fun = star_b4, show_heatmap_legend = FALSE,
  width = unit(5, "cm"), height = unit(0.5 * nrow(brain_cur4), "cm"), border = TRUE, border_gp = gpar(col = "black", lwd = 1))
# 범례
make_legends_b <- function(with_text = TRUE) {
  ttl <- function(x) if (with_text) x else ""
  packLegend(
    Legend(labels = if (with_text) brain_order4 else rep("", length(brain_order4)), legend_gp = gpar(fill = brain_colors4[brain_order4]),
           title = ttl("Category"), title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 10), grid_height = unit(4, "mm"), grid_width = unit(4, "mm")),
    Legend(labels = if (with_text) grp_lev_b else rep("", length(grp_lev_b)), legend_gp = gpar(fill = grp_col_b[grp_lev_b]),
           title = ttl("Cell group"), title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 10), grid_height = unit(4, "mm"), grid_width = unit(4, "mm")),
    Legend(col_fun = col_fun_b, title = ttl("log2FC"), at = c(-3, -2, -1, 0, 1, 2, 3),
           labels = if (with_text) c(-3, -2, -1, 0, 1, 2, 3) else rep("", 7),
           title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 10), legend_height = unit(4, "cm")),
    direction = "vertical", gap = unit(5, "mm"))
}
lg_b <- make_legends_b(TRUE); lg_b_nt <- make_legends_b(FALSE)
pdf(file.path(brain_final, "brain_aging_heatmap_FINAL.pdf"), width = 6, height = 9)
draw(ht_b, column_title = "Brain (HGPS_EGM vs HGPS control)", column_title_gp = gpar(fontsize = 12, fontface = "bold")); dev.off()
png(file.path(brain_final, "brain_aging_heatmap_FINAL.png"), width = 6, height = 9, units = "in", res = 300)
draw(ht_b, column_title = "Brain (HGPS_EGM vs HGPS control)", column_title_gp = gpar(fontsize = 12, fontface = "bold")); dev.off()
pdf(file.path(brain_final, "brain_aging_heatmap_FINAL_notext.pdf"), width = 4.5, height = 8.5); draw(ht_b_nt); dev.off()
pdf(file.path(brain_final, "brain_aging_heatmap_FINAL_legend.pdf"), width = 3.2, height = 6)
grid.newpage(); draw(lg_b, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "centre"); dev.off()
pdf(file.path(brain_final, "brain_aging_heatmap_FINAL_legend_notext.pdf"), width = 2, height = 6)
grid.newpage(); draw(lg_b_nt, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "centre"); dev.off()
saveRDS(bfc4, file.path(brain_final, "brain_aging_log2fc_matrix_final.rds"))
saveRDS(bpv4, file.path(brain_final, "brain_aging_padj_matrix_final.rds"))
saveRDS(brain_cur4, file.path(brain_final, "brain_aging_gene_categories_final.rds"))
write.csv(as.data.frame(bfc4) %>% tibble::rownames_to_column("gene") %>% mutate(Category = brain_cur4$Fig, .after = gene),
          file.path(brain_final, "brain_aging_heatmap_log2FC_matrix.csv"), row.names = FALSE)
write.csv(ball4 %>% left_join(brain_cur4, by = "gene") %>% mutate(Fig = factor(Fig, levels = brain_order4), CellType = factor(CellType, levels = bcols)) %>%
            arrange(Fig, match(gene, brain_cur4$gene), CellType) %>% select(Category = Fig, gene, CellType, baseMean, log2FoldChange, pvalue, padj),
          file.path(brain_final, "brain_aging_heatmap_DESeq2_results.csv"), row.names = FALSE)
save(bfc4, bpv4, brain_cur4, brain_order4, brain_colors4, bcols, brain_grp, col_split_b, grp_lev_b, grp_col_b, col_fun_b, star_b4,
     ht_b, ht_b_nt, lg_b, lg_b_nt, brain_res, bsig, bgene_tbl, brain_markers, brain_deg,
     file = file.path(brain_final, "brain_aging_final_workspace.RData"))
write.csv(brain_res, file.path(brain_out, "brain_aging_markers_all_celltypes.csv"), row.names = FALSE)
print(list.files(brain_final))

