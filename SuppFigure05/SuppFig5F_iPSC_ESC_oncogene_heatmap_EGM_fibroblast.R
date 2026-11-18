# OncoKB core oncogenes in egm_blood: EGM vs GFP (Control)
# Same approach as EGM_tissue/oncogene_analysis/Mouse/oncogene_heatmap/oncogene_heatmap.R
# - Counts   : counts_symbol.csv (STAR -> RSEM, mouse gene symbols); GFP-1..3 vs EGM-1..3 only
# - Mapping  : human OncoKB -> mouse via Ensembl orthologs (reuses the tissue analysis mapping table)
# - DE       : DESeq2, filter >=10 counts in >=3 samples
# - Enrichment: % oncogenes up (padj<0.05) vs expression-matched non-oncogene background (10 bins, 2000 draws)
# - Core set : OncoKB ONCOGENE listed in all 7 aggregated resources (n = 36); not-expressed genes dropped

suppressMessages({
  library(DESeq2); library(dplyr); library(ComplexHeatmap); library(circlize)
})
set.seed(1)
bbase <- "4FP_plasma/egm_blood"
bout  <- file.path(bbase, "oncogene_heatmap")
tis   <- "EGM_tissue"
ht_opt$message <- FALSE

okb   <- read.delim(file.path(tis, "cancerGeneList.tsv"), check.names = FALSE)
core7 <- okb$`Hugo Symbol`[okb$`Gene Type` == "ONCOGENE" & okb$`# of occurrence within resources (Column K-P)` == 7]
om2   <- read.csv(file.path(tis, "oncogene_analysis/Mouse/oncogene_heatmap/oncokb_oncogene_mouse_mapping.csv"))

# ---- DESeq2 ----
bc    <- read.csv(file.path(bbase, "counts_symbol.csv"), row.names = 1, check.names = FALSE)
bcols <- c(paste0("GFP-", 1:3), paste0("EGM-", 1:3))
meta  <- data.frame(group = factor(rep(c("Control", "EGM"), each = 3), levels = c("Control", "EGM")), row.names = bcols)
dds <- DESeqDataSetFromMatrix(round(as.matrix(bc[, bcols])), meta, ~group)
dds <- DESeq(dds[rowSums(counts(dds) >= 10) >= 3, ], quiet = TRUE)
br  <- as.data.frame(results(dds, contrast = c("group", "EGM", "Control"))); br$mm_symbol <- rownames(br)
vsd <- vst(dds, blind = FALSE)

# ---- oncogene results + expression-matched enrichment ----
r <- br %>% filter(!is.na(padj)) %>%
  mutate(onc = mm_symbol %in% om2$mm_symbol,
         dir = ifelse(padj < 0.05 & log2FoldChange > 0, "Up", ifelse(padj < 0.05 & log2FoldChange < 0, "Down", "NS")),
         bin = cut(log10(baseMean + 1), quantile(log10(baseMean + 1), 0:10 / 10), include.lowest = TRUE))
ro <- r %>% filter(onc) %>% inner_join(om2, by = "mm_symbol") %>%
  mutate(same = toupper(mm_symbol) == hs_symbol) %>% group_by(hs_symbol) %>%
  arrange(desc(same), desc(baseMean), .by_group = TRUE) %>% slice(1) %>% ungroup()
bg <- r %>% filter(!onc)
permtest <- function(ids) {
  tab <- table(r$bin[match(ids, r$mm_symbol)])
  p <- replicate(2000, { s <- match(unlist(lapply(names(tab), function(b) sample(bg$mm_symbol[bg$bin == b], tab[b]))), r$mm_symbol)
                         mean(r$dir[s] == "Up") })
  obs <- mean(r$dir[match(ids, r$mm_symbol)] == "Up"); c(obs = 100 * obs, bg = 100 * mean(p), p = mean(p >= obs))
}
cb <- ro %>% filter(hs_symbol %in% core7)
a <- permtest(ro$mm_symbol); b <- permtest(cb$mm_symbol)
summ <- data.frame(set = c("all_oncogenes", "core36"), n = c(nrow(ro), nrow(cb)),
  up = c(sum(ro$dir == "Up"), sum(cb$dir == "Up")), down = c(sum(ro$dir == "Down"), sum(cb$dir == "Down")),
  ns = c(sum(ro$dir == "NS"), sum(cb$dir == "NS")),
  up_pct = round(c(a["obs"], b["obs"]), 1), matched_bg_up_pct = round(c(a["bg"], b["bg"]), 1), perm_p = c(a["p"], b["p"]))
out_cols <- c("hs_symbol", "mm_symbol", "type", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "dir")
write.csv(ro[, out_cols] %>% arrange(padj), file.path(bout, "oncogene_DE_EGM_vs_GFP.csv"), row.names = FALSE)
write.csv(cb[, out_cols] %>% arrange(hs_symbol), file.path(bout, "core_oncogenes_DE_EGM_vs_GFP.csv"), row.names = FALSE)
write.csv(summ, file.path(bout, "oncogene_enrichment_summary.csv"), row.names = FALSE)

# ---- core oncogene heatmap (same style as tissue figures) ----
cb <- cb %>% arrange(mm_symbol)
v  <- assay(vsd)[cb$mm_symbol, bcols, drop = FALSE]
R  <- v - rowMeans(v[, 1:3])                          # VST minus GFP mean
d  <- table(factor(cb$dir, levels = c("Up", "NS", "Down")))
grp <- factor(rep(c("GFP (Control)", "EGM"), each = 3), levels = c("GFP (Control)", "EGM"))
vir <- colorRamp2(seq(-2, 2, length.out = 256), viridisLite::viridis(256))
hb <- Heatmap(R, name = "log2 ratio\nvs. GFP mean", col = vir,
  cluster_rows = FALSE, cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation(Group = grp, col = list(Group = c("GFP (Control)" = "#d1d2d7", EGM = "#6d98f8")),
    show_annotation_name = FALSE, simple_anno_size = unit(3, "mm"),
    annotation_legend_param = list(Group = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)))),
  bottom_annotation = HeatmapAnnotation(summary = anno_block(labels = sprintf("Up %d · Down %d\nNS %d", d["Up"], d["Down"], d["NS"]),
    labels_gp = gpar(fontsize = 7), gp = gpar(fill = NA, col = NA), height = unit(8, "mm")), show_annotation_name = FALSE),
  show_column_names = FALSE, row_names_side = "right", row_names_gp = gpar(fontsize = 9.5),
  width = unit(6 * 4.2, "mm"), height = unit(nrow(R) * 3.6, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.4),
  heatmap_legend_param = list(at = -2:2, title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)))
f <- file.path(bout, "Fig_core_oncogenes_heatmap_egm_blood")
for (dev in c("pdf", "png")) {
  if (dev == "pdf") pdf(paste0(f, ".pdf"), width = 3.6, height = 6.6) else png(paste0(f, ".png"), width = 3.6, height = 6.6, units = "in", res = 300)
  draw(hb, column_title = sprintf("Core oncogenes\n(OncoKB consensus, n = %d)\nEGM vs GFP", nrow(R)),
       column_title_gp = gpar(fontsize = 11), merge_legend = TRUE, padding = unit(c(2, 4, 2, 6), "mm"))
  dev.off()
}
print(summ)

# ---- iPSC/ESC-associated oncogenes (in vitro readout) ----
# 40 genes (`cancer_genes` in code_SOT1.R):
# iPSC/ESC markers (SOT C1 OSK reprogramming data, FindMarkers ips+esc vs rest; ipsc_esc.csv) that are in the OncoKB list
sot_genes <- c("Kit", "Sf3b1", "Ccne1", "Mdm4", "Mycn", "Pim1", "Srsf2", "Xpo1", "Eif4a2", "Fgf4", "Etv5", "Hdac1", "Sgk1",
  "Cenpa", "Dek", "Hsp90aa1", "Ncoa3", "Nup98", "Set", "Cad", "Gab1", "Hsp90ab1", "Srsf3", "Bub1b", "Ccnb1ip1", "Clp1",
  "Cnbp", "Folr1", "Inpp5d", "Kat6b", "Kdm2b", "Mthfd2", "Naca", "Parp2", "Pum1", "Rnaseh2a", "Sall4", "Tcea1", "Tfrc", "Top2a")
sg <- r %>% filter(mm_symbol %in% sot_genes) %>% arrange(mm_symbol)   # expressed & tested genes only
s2 <- permtest(sg$mm_symbol)
summ_sot <- data.frame(set = "ipsc_esc_oncogenes", n_list = length(sot_genes), n_expressed = nrow(sg),
  not_expressed = paste(setdiff(sot_genes, sg$mm_symbol), collapse = ";"),
  up = sum(sg$dir == "Up"), down = sum(sg$dir == "Down"), ns = sum(sg$dir == "NS"),
  up_pct = round(s2["obs"], 1), matched_bg_up_pct = round(s2["bg"], 1), perm_p = s2["p"], row.names = NULL)
write.csv(sg %>% select(mm_symbol, baseMean, log2FoldChange, lfcSE, pvalue, padj, dir),
  file.path(bout, "ipsc_esc_oncogenes_DE_EGM_vs_GFP.csv"), row.names = FALSE)
write.csv(summ_sot, file.path(bout, "ipsc_esc_oncogenes_summary.csv"), row.names = FALSE)

v2 <- assay(vsd)[sg$mm_symbol, bcols, drop = FALSE]
R2 <- v2 - rowMeans(v2[, 1:3])                        # VST minus GFP mean
d2 <- table(factor(sg$dir, levels = c("Up", "NS", "Down")))
hs2 <- Heatmap(R2, name = "log2 ratio\nvs. GFP mean", col = vir,
  cluster_rows = FALSE, cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation(Group = grp, col = list(Group = c("GFP (Control)" = "#d1d2d7", EGM = "#6d98f8")),
    show_annotation_name = FALSE, simple_anno_size = unit(3, "mm"),
    annotation_legend_param = list(Group = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)))),
  bottom_annotation = HeatmapAnnotation(summary = anno_block(labels = sprintf("Up %d · Down %d\nNS %d", d2["Up"], d2["Down"], d2["NS"]),
    labels_gp = gpar(fontsize = 7), gp = gpar(fill = NA, col = NA), height = unit(8, "mm")), show_annotation_name = FALSE),
  show_column_names = FALSE, row_names_side = "right", row_names_gp = gpar(fontsize = 9.5),
  width = unit(6 * 4.2, "mm"), height = unit(nrow(R2) * 3.6, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.4),
  heatmap_legend_param = list(at = -2:2, title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)))
f2 <- file.path(bout, "Fig_ipsc_esc_oncogenes_heatmap_egm_blood")
for (dev in c("pdf", "png")) {
  if (dev == "pdf") pdf(paste0(f2, ".pdf"), width = 3.6, height = 7.0) else png(paste0(f2, ".png"), width = 3.6, height = 7.0, units = "in", res = 300)
  draw(hs2, column_title = sprintf("iPSC/ESC-associated\noncogenes (n = %d)\nEGM vs GFP", nrow(R2)),
       column_title_gp = gpar(fontsize = 11), merge_legend = TRUE, padding = unit(c(2, 4, 2, 6), "mm"))
  dev.off()
}
print(summ_sot)

# ---- rotated (90° counter-clockwise) version: genes as columns (labels at bottom), EGM rows on top ----
draw_rot <- function(M, legend_title, file) {
  Mt <- t(M)[rev(bcols), ]
  gr <- factor(ifelse(grepl("^EGM", rownames(Mt)), "EGM", "GFP (Control)"), levels = c("GFP (Control)", "EGM"))
  h <- Heatmap(Mt, name = legend_title, col = vir, cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, column_names_side = "bottom", column_names_rot = 90, column_names_gp = gpar(fontsize = 9.5),
    left_annotation = rowAnnotation(Group = gr, col = list(Group = c("GFP (Control)" = "#d1d2d7", EGM = "#6d98f8")),
      show_annotation_name = FALSE, simple_anno_size = unit(3, "mm"),
      annotation_legend_param = list(Group = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)))),
    width = unit(ncol(Mt) * 4, "mm"), height = unit(nrow(Mt) * 4.2, "mm"), rect_gp = gpar(col = "white", lwd = 0.4),
    heatmap_legend_param = list(at = -2:2, title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)))
  for (dev in c("pdf", "png")) {
    if (dev == "pdf") pdf(paste0(file, ".pdf"), width = 7.6, height = 2.4) else png(paste0(file, ".png"), width = 7.6, height = 2.4, units = "in", res = 300)
    draw(h, merge_legend = TRUE, padding = unit(c(2, 4, 2, 4), "mm"))
    dev.off()
  }
}
draw_rot(R2, "log2 ratio\nvs. GFP mean", file.path(bout, "Fig_ipsc_esc_oncogenes_heatmap_egm_blood_rotated"))
