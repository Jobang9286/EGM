source("Fig2P-Q_00_aging_gene_sets_setup.R")   # root, aging_sets, protective, ambiguous

# ---- liver: 입력 데이터 확인 ----
liver_dir <- file.path(root, "integrated_analysis_for_lizard_factors/liver/liver_deg_results")
print(list.files(liver_dir, pattern = "csv$|rds$"))
liver_deg <- readRDS(file.path(liver_dir, "deg_results_list.rds"))
cat("\nzones:", names(liver_deg), "\n")
for (z in names(liver_deg)) cat(sprintf("  %-12s genes=%d\n", z, nrow(liver_deg[[z]])))
str(liver_deg[[1]], max.level = 1)
Z_liver <- readRDS(file.path(liver_dir, "Z_celltype_expression.rds"))
cat("\nZ_celltype_expression:", class(Z_liver), names(Z_liver), "\n")
for (z in names(Z_liver)) cat(sprintf("  %-12s dim=%s  cols=%s\n", z, paste(dim(Z_liver[[z]]), collapse="x"), paste(colnames(Z_liver[[z]]), collapse=",")))

## liver aging screen: aging_sets + liver function / regeneration markers, per-zone DESeq2 hits.
options(width = 230, max.print = 8000)
liver_out <- file.path(root, "integrated_analysis_for_lizard_factors/liver/aging_screening")
dir.create(liver_out, showWarnings = FALSE, recursive = TRUE)
zones <- c("Periportal", "Midzonal", "Pericentral")

# 간 특이 기능/재생 마커 추가 (kidney의 Tubular_Function 자리)
liver_sets <- c(aging_sets[setdiff(names(aging_sets), "Tubular_Function")], list(
  Hepatocyte_Function = c("Alb","Ttr","Apoa1","Apob","Apoc3","Serpina1a","Serpina1b","Trf","Tf","F2","F7","Fga","Fgb","Fgg","Ahsg","Hp","Hpx",
                          "Cyp1a2","Cyp2e1","Cyp3a11","Cyp7a1","Cyp8b1","Cyp2c29","Ugt1a1","Gsta3","Nr1h4","Slco1b2","Abcb11","Slc10a1",
                          "G6pc","Pck1","Glul","Cps1","Otc","Arg1","Ass1","Asl","Tat","Hnf4a","Hnf1a","Foxa2","Nr1i2","Nr1i3","Ppara"),
  Liver_Regeneration = c("Hgf","Met","Egfr","Areg","Hbegf","Il6","Stat3","Ccnd1","Myc","Foxm1","Yap1","Ctgf","Ccn2","Sox9","Epcam","Krt19","Krt7",
                         "Prom1","Lgr5","Afp","Tbx3","Wnt2","Rspo3","Axin2","Glul","Notch2","Jag1","Onecut1","Cebpa","Cebpb","Nr5a2","Ppargc1a")))
liver_markers <- bind_rows(lapply(names(liver_sets), function(n) data.frame(Category = n, Gene = liver_sets[[n]]))) %>% distinct(Gene, .keep_all = TRUE)
protective_liver <- c(protective, liver_sets$Hepatocyte_Function)
liver_markers$Aging_Role <- ifelse(liver_markers$Gene %in% ambiguous, "Ambiguous",
                            ifelse(liver_markers$Gene %in% protective_liver, "Protective", "Pro-aging"))
liver_markers$Aging_Role[liver_markers$Category == "Liver_Regeneration"] <- "Protective"
write.csv(liver_markers, file.path(liver_out, "liver_aging_gene_sets.csv"), row.names = FALSE)
cat("gene set:", nrow(liver_markers), "genes /", length(liver_sets), "categories\n")

liver_res <- bind_rows(lapply(zones, function(z) liver_deg[[z]] %>% filter(gene %in% liver_markers$Gene) %>% mutate(Zone = z))) %>%
  left_join(liver_markers, by = c("gene" = "Gene")) %>%
  mutate(Effect = case_when(is.na(padj) | padj >= 0.05 ~ "NS",
                            Aging_Role == "Pro-aging"  & log2FoldChange < 0 ~ "Improved",
                            Aging_Role == "Protective" & log2FoldChange > 0 ~ "Improved",
                            Aging_Role == "Ambiguous" ~ "Ambiguous", TRUE ~ "Worsened"))
write.csv(liver_res, file.path(liver_out, "liver_aging_markers_all_zones.csv"), row.names = FALSE)
lsig <- liver_res %>% filter(padj < 0.05)
cat("detected:", n_distinct(liver_res$gene), "/", nrow(liver_markers), "| significant genes:", n_distinct(lsig$gene), "| hits:", nrow(lsig), "\n\n")
print(as.data.frame(lsig %>% group_by(Category) %>% summarise(n_genes = n_distinct(gene), hits = n(), Up = sum(log2FoldChange > 0), Down = sum(log2FoldChange < 0),
      Improved = sum(Effect == "Improved"), Worsened = sum(Effect == "Worsened"), .groups = "drop") %>% arrange(desc(hits))), row.names = FALSE)
cat("\n--- by zone ---\n")
print(as.data.frame(lsig %>% group_by(Zone) %>% summarise(n = n_distinct(gene), Up = sum(log2FoldChange > 0), Down = sum(log2FoldChange < 0),
      Improved = sum(Effect == "Improved"), Worsened = sum(Effect == "Worsened"), .groups = "drop")), row.names = FALSE)

lgene_tbl <- lsig %>% group_by(Category, gene, Aging_Role) %>%
  summarise(nZone = n(), meanFC = round(mean(log2FoldChange), 2), minPadj = signif(min(padj), 2), Effect = names(sort(table(Effect), decreasing = TRUE))[1],
            Zones = paste0(substr(Zone, 1, 4), "(", ifelse(log2FoldChange > 0, "+", "-"), ")", collapse = " "), .groups = "drop") %>%
  arrange(Category, desc(nZone), minPadj)
write.csv(lgene_tbl, file.path(liver_out, "liver_aging_markers_significant_by_gene.csv"), row.names = FALSE)
cat("\n--- significant genes by category ---\n"); print(as.data.frame(lgene_tbl), row.names = FALSE)

## first curated gene list and the log2FC / padj matrices.
liver_cur <- tibble::tribble(
  ~gene, ~Fig,
  # ---- aging features: 감소 ----
  "Cdkn1a","Senescence", "Cdkn1c","Senescence",
  "Igfbp4","SASP", "Mmp14","SASP", "Igfbp7","SASP",
  "Stat3","Inflammation",
  "Lcn2","Liver_Injury", "Cp","Liver_Injury", "Cst3","Liver_Injury", "Tnfrsf12a","Liver_Injury",
  "Bnip3","Cell_Stress", "Atf4","Cell_Stress",
  "Akt1","Nutrient_Sensing", "Igf1","Nutrient_Sensing",
  # ---- repair / function: 증가 ----
  "Ulk1","Autophagy", "Atg5","Autophagy", "Nbr1","Autophagy",
  "Rad23b","DNA_Repair", "Xpc","DNA_Repair",
  "Ppara","Mitochondrial", "Ppargc1b","Mitochondrial", "Cpt1a","Mitochondrial",
  "Gclc","Antioxidant", "Cat","Antioxidant",
  "Klb","Longevity", "Fgf21","Longevity",
  "Cyp8b1","Hepatocyte_Function", "Apob","Hepatocyte_Function", "Nr1h4","Hepatocyte_Function",
  "Cebpa","Liver_Regeneration", "Nr5a2","Liver_Regeneration", "Yap1","Liver_Regeneration")
liver_order  <- c("Senescence","SASP","Inflammation","Liver_Injury","Cell_Stress","Nutrient_Sensing",
                  "Autophagy","DNA_Repair","Mitochondrial","Antioxidant","Longevity","Hepatocyte_Function","Liver_Regeneration")
liver_colors <- c(Senescence = "#3C5488", SASP = "#E64B35", Inflammation = "#F39B7F", Liver_Injury = "#91A546", Cell_Stress = "#8491B4",
                  Nutrient_Sensing = "#B09C85", Autophagy = "#DC0000", DNA_Repair = "#7E6148", Mitochondrial = "#00A087",
                  Antioxidant = "#FFB300", Longevity = "#4DBBD5", Hepatocyte_Function = "#66C2A5", Liver_Regeneration = "#631879")
stopifnot(all(liver_cur$gene %in% lgene_tbl$gene))
lfc <- lsig %>% filter(gene %in% liver_cur$gene) %>% select(gene, Zone, log2FoldChange) %>% pivot_wider(names_from = Zone, values_from = log2FoldChange) %>% tibble::column_to_rownames("gene") %>% as.matrix()
lpv <- lsig %>% filter(gene %in% liver_cur$gene) %>% select(gene, Zone, padj) %>% pivot_wider(names_from = Zone, values_from = padj) %>% tibble::column_to_rownames("gene") %>% as.matrix()
lfc[is.na(lfc)] <- 0
lfc <- lfc[liver_cur$gene, zones]; lpv <- lpv[liver_cur$gene, zones]

## remove Inflammation and Nutrient_Sensing.
liver_cur2 <- liver_cur %>% filter(!Fig %in% c("Inflammation", "Nutrient_Sensing"))
liver_order2 <- setdiff(liver_order, c("Inflammation", "Nutrient_Sensing"))
cat(nrow(liver_cur2), "genes,", length(liver_order2), "categories (제거:", setdiff(liver_cur$gene, liver_cur2$gene), ")\n\n")

## zone colours.
zone_colors <- c(Periportal = "#66C2A5", Midzonal = "#FC8D62", Pericentral = "#8DA0CB")

## output folder and colour function.
liver_final <- file.path(root, "integrated_analysis_for_lizard_factors/integrated_analysis/liver/final_aging")
dir.create(liver_final, showWarnings = FALSE, recursive = TRUE)
col_fun_l <- colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))

## proliferation / hepatic progenitor marker screen.
prolif_stem <- list(
  Proliferation = c("Mki67","Pcna","Top2a","Ccnb1","Ccnb2","Ccnd1","Ccne1","Ccne2","Cdk1","Cdk2","Cdk4","Cdk6","Ccna2","Birc5","Foxm1","Aurka","Aurkb","Plk1",
                    "E2f1","E2f2","Mcm2","Mcm3","Mcm4","Mcm5","Mcm6","Mcm7","Rrm2","Tk1","Tyms","Cdc20","Bub1","Kif11","Nusap1","Ube2c","Cenpf","Cenpa"),
  Hepatic_Progenitor = c("Sox9","Epcam","Krt19","Krt7","Krt8","Krt18","Prom1","Lgr5","Lgr4","Afp","Cd24a","Onecut1","Onecut2","Hnf1b","Tbx3","Axin2","Rspo3",
                         "Notch2","Jag1","Spp1","Cd44","Thy1","Ncam1","Dlk1","Trop2","Tacstd2","Mfge8","Clu","Anxa2","Sall4","Yap1","Wwtr1","Ctgf","Ccn2","Cyr61","Ccn1"))
ps_markers <- bind_rows(lapply(names(prolif_stem), function(n) data.frame(Category = n, Gene = prolif_stem[[n]]))) %>% distinct(Gene, .keep_all = TRUE)
ps_res <- bind_rows(lapply(zones, function(z) liver_deg[[z]] %>% filter(gene %in% ps_markers$Gene) %>% mutate(Zone = z))) %>% left_join(ps_markers, by = c("gene" = "Gene"))

## replace Longevity (Klb, Fgf21) by Proliferation (Mcm5, E2f1, Mki67).
liver_cur3 <- liver_cur2 %>% filter(Fig != "Longevity") %>%
  bind_rows(tibble::tibble(gene = c("Mcm5","E2f1","Mki67"), Fig = "Proliferation")) %>%
  mutate(Fig = factor(Fig, levels = c("Senescence","SASP","Liver_Injury","Cell_Stress","Autophagy","DNA_Repair",
                                      "Mitochondrial","Antioxidant","Proliferation","Hepatocyte_Function","Liver_Regeneration"))) %>%
  arrange(Fig)
liver_order3 <- levels(liver_cur3$Fig)
liver_colors3 <- c(liver_colors[setdiff(liver_order3, "Proliferation")], Proliferation = "#4DBBD5")   # Longevity가 쓰던 색 재사용
lsig_all <- bind_rows(lsig, ps_res %>% filter(padj < 0.05, gene %in% c("Mcm5","E2f1","Mki67")) %>% mutate(Aging_Role = "Protective", Effect = "Improved"))
lfc3 <- lsig_all %>% filter(gene %in% liver_cur3$gene) %>% select(gene, Zone, log2FoldChange) %>% pivot_wider(names_from = Zone, values_from = log2FoldChange) %>% tibble::column_to_rownames("gene") %>% as.matrix()
lpv3 <- lsig_all %>% filter(gene %in% liver_cur3$gene) %>% select(gene, Zone, padj) %>% pivot_wider(names_from = Zone, values_from = padj) %>% tibble::column_to_rownames("gene") %>% as.matrix()
lfc3[is.na(lfc3)] <- 0; lfc3 <- lfc3[liver_cur3$gene, zones]; lpv3 <- lpv3[liver_cur3$gene, zones]

## drop Cell_Stress; merge Proliferation + Liver_Regeneration into Hepatic_Progenitor_Stem (+ Hnf1b).
liver_cur4 <- liver_cur3 %>% mutate(Fig = as.character(Fig)) %>%
  filter(Fig != "Cell_Stress") %>%
  mutate(Fig = ifelse(Fig %in% c("Proliferation", "Liver_Regeneration"), "Hepatic_Progenitor_Stem", Fig))
# Hnf1b: hepatic progenitor/biliary TF, 3 zone 모두 유의(+0.60) -> 새 카테고리 이름에 맞아 추가 (원치 않으면 빼기)
liver_cur4 <- bind_rows(liver_cur4, tibble::tibble(gene = "Hnf1b", Fig = "Hepatic_Progenitor_Stem"))
liver_order4 <- c("Senescence","SASP","Liver_Injury","Autophagy","DNA_Repair","Mitochondrial","Antioxidant","Hepatocyte_Function","Hepatic_Progenitor_Stem")
liver_cur4 <- liver_cur4 %>% mutate(Fig = factor(Fig, levels = liver_order4)) %>% arrange(Fig)
liver_colors4 <- c(liver_colors[c("Senescence","SASP","Liver_Injury","Autophagy","DNA_Repair","Mitochondrial","Antioxidant","Hepatocyte_Function")],
                   Hepatic_Progenitor_Stem = "#631879")
lsig_all2 <- bind_rows(lsig_all, ps_res %>% filter(padj < 0.05, gene == "Hnf1b") %>% mutate(Aging_Role = "Protective", Effect = "Improved")) %>% distinct(gene, Zone, .keep_all = TRUE)
lfc4 <- lsig_all2 %>% filter(gene %in% liver_cur4$gene) %>% select(gene, Zone, log2FoldChange) %>% pivot_wider(names_from = Zone, values_from = log2FoldChange) %>% tibble::column_to_rownames("gene") %>% as.matrix()
lpv4 <- lsig_all2 %>% filter(gene %in% liver_cur4$gene) %>% select(gene, Zone, padj) %>% pivot_wider(names_from = Zone, values_from = padj) %>% tibble::column_to_rownames("gene") %>% as.matrix()
lfc4[is.na(lfc4)] <- 0; lfc4 <- lfc4[liver_cur4$gene, zones]; lpv4 <- lpv4[liver_cur4$gene, zones]
star_fun_l4 <- function(j, i, x, y, w, h, fill) { p <- lpv4[i, j]
  if (!is.na(p)) grid.text(if (p < 0.001) "***" else if (p < 0.01) "**" else "*", x, y, gp = gpar(fontsize = 8, col = "black", fontface = "bold")) }
ht_liver4 <- Heatmap(lfc4, name = "log2FC", col = col_fun_l, cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 9),
  top_annotation = HeatmapAnnotation(Zonation = factor(zones, levels = zones), col = list(Zonation = zone_colors),
                                     show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(3, "mm")),
  column_title = "Hepatocyte zonation\n(portal → central)", column_title_gp = gpar(fontsize = 9.5, fontface = "bold"),
  left_annotation = rowAnnotation(Category = liver_cur4$Fig, col = list(Category = liver_colors4),
                                  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"), show_legend = FALSE, width = unit(0.5, "cm")),
  rect_gp = gpar(col = "grey90", lwd = 0.3), cell_fun = star_fun_l4, show_heatmap_legend = FALSE,
  width = unit(3.6, "cm"), height = unit(0.5 * nrow(liver_cur4), "cm"), border = TRUE, border_gp = gpar(col = "black", lwd = 1))

## final save.
ht_liver_nt4 <- Heatmap(lfc4, name = "log2FC", col = col_fun_l, cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = FALSE, column_title = NULL,
  top_annotation = HeatmapAnnotation(Zonation = factor(zones, levels = zones), col = list(Zonation = zone_colors),
                                     show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(3, "mm")),
  left_annotation = rowAnnotation(Category = liver_cur4$Fig, col = list(Category = liver_colors4),
                                  show_annotation_name = FALSE, show_legend = FALSE, width = unit(0.5, "cm")),
  rect_gp = gpar(col = "grey90", lwd = 0.3), cell_fun = star_fun_l4, show_heatmap_legend = FALSE,
  width = unit(3.6, "cm"), height = unit(0.5 * nrow(liver_cur4), "cm"), border = TRUE, border_gp = gpar(col = "black", lwd = 1))
make_legends_l4 <- function(with_text = TRUE) {
  ttl <- function(x) if (with_text) x else ""
  packLegend(
    Legend(labels = if (with_text) liver_order4 else rep("", length(liver_order4)), legend_gp = gpar(fill = liver_colors4[liver_order4]),
           title = ttl("Category"), title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 10), grid_height = unit(4, "mm"), grid_width = unit(4, "mm")),
    Legend(labels = if (with_text) zones else rep("", length(zones)), legend_gp = gpar(fill = zone_colors[zones]),
           title = ttl("Zonation"), title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 10), grid_height = unit(4, "mm"), grid_width = unit(4, "mm")),
    Legend(col_fun = col_fun_l, title = ttl("log2FC"), at = c(-3, -2, -1, 0, 1, 2, 3),
           labels = if (with_text) c(-3, -2, -1, 0, 1, 2, 3) else rep("", 7),
           title_gp = gpar(fontsize = 10, fontface = "bold"), labels_gp = gpar(fontsize = 10), legend_height = unit(4, "cm")),
    direction = "vertical", gap = unit(5, "mm"))
}
lg_l4 <- make_legends_l4(TRUE); lg_l4_nt <- make_legends_l4(FALSE)
# 이전 29-gene(Longevity 포함) 버전 파일 제거 후 새 버전 저장
# file.remove(list.files(liver_final, pattern = "^liver_aging_heatmap_29genes", full.names = TRUE))
pdf(file.path(liver_final, "liver_aging_heatmap_FINAL.pdf"), width = 5.5, height = 9.5)
draw(ht_liver4, column_title = "Liver (HGPS_EGM vs HGPS control)", column_title_gp = gpar(fontsize = 12, fontface = "bold")); dev.off()
png(file.path(liver_final, "liver_aging_heatmap_FINAL.png"), width = 5.5, height = 9.5, units = "in", res = 300)
draw(ht_liver4, column_title = "Liver (HGPS_EGM vs HGPS control)", column_title_gp = gpar(fontsize = 12, fontface = "bold")); dev.off()
pdf(file.path(liver_final, "liver_aging_heatmap_FINAL_notext.pdf"), width = 4, height = 9); draw(ht_liver_nt4); dev.off()
pdf(file.path(liver_final, "liver_aging_heatmap_FINAL_legend.pdf"), width = 3.2, height = 6.5)
grid.newpage(); draw(lg_l4, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "centre"); dev.off()
pdf(file.path(liver_final, "liver_aging_heatmap_FINAL_legend_notext.pdf"), width = 2, height = 6.5)
grid.newpage(); draw(lg_l4_nt, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "centre"); dev.off()
saveRDS(lfc4, file.path(liver_final, "liver_aging_log2fc_matrix_final.rds"))
saveRDS(lpv4, file.path(liver_final, "liver_aging_padj_matrix_final.rds"))
saveRDS(liver_cur4, file.path(liver_final, "liver_aging_gene_categories_final.rds"))
write.csv(as.data.frame(lfc4) %>% tibble::rownames_to_column("gene") %>% mutate(Category = as.character(liver_cur4$Fig), .after = gene),
          file.path(liver_final, "liver_aging_heatmap_log2FC_matrix.csv"), row.names = FALSE)
all_liver_final <- bind_rows(lapply(zones, function(z) liver_deg[[z]] %>% filter(gene %in% liver_cur4$gene) %>% mutate(Zone = z))) %>%
  left_join(liver_cur4, by = "gene") %>% mutate(Zone = factor(Zone, levels = zones)) %>%
  arrange(Fig, match(gene, liver_cur4$gene), Zone) %>% select(Category = Fig, gene, Zone, baseMean, log2FoldChange, pvalue, padj)
write.csv(all_liver_final, file.path(liver_final, "liver_aging_heatmap_DESeq2_results.csv"), row.names = FALSE)
save(lfc4, lpv4, liver_cur4, liver_order4, liver_colors4, zones, zone_colors, col_fun_l, star_fun_l4,
     ht_liver4, ht_liver_nt4, lg_l4, lg_l4_nt, liver_res, lsig_all2, lgene_tbl, liver_markers,
     file = file.path(liver_final, "liver_aging_final_workspace.RData"))
print(list.files(liver_final))

