## ---- (1) setup -------------------------------------------------------------
suppressPackageStartupMessages({library(ComplexHeatmap); library(circlize); library(grid); library(ggplot2); library(ggpubr); library(patchwork)})
ROOT <- "human_egm_261007/egm_score"; FIG <- file.path(ROOT, "figures"); DAT <- file.path(ROOT, "data")
SRC  <- file.path(dirname(ROOT), "human_egm_results")
tpm <- read.delim(file.path(SRC, "human_EGM_TPM.tsv"), row.names = 1, check.names = FALSE)
map <- read.delim(file.path(SRC, "human_gene_id_to_symbol.tsv")); sym <- setNames(map$Gene.Name, map$Gene.ID)
de  <- read.csv(file.path(DAT, "EGM_vs_Control_signature.csv"))
tpm$symbol <- unname(sym[rownames(tpm)]); tpm <- tpm[!is.na(tpm$symbol) & tpm$symbol != "" & !grepl("^ENSG", tpm$symbol), ]
tpm$mean <- rowMeans(tpm[, 1:6]); tpm <- tpm[order(tpm$symbol, -tpm$mean), ]; tpm <- tpm[!duplicated(tpm$symbol), ]
TPM <- as.matrix(tpm[, 1:6]); rownames(TPM) <- tpm$symbol
grp <- factor(ifelse(grepl("^Control", colnames(TPM)), "Control", "EGM"), c("Control", "EGM"))

## ---- (1) colours / legend ----------------------------------------------------
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#2166AC", "white", "#B2182B"))
catcol  <- c("Senescence/SASP" = "#8C6BB1", "Fibrosis" = "#D08C3A", "Proliferation" = "#4D9F6B")
grpcol  <- c(Control = "#9A9A9A", EGM = "#3E8FBF")
leg <- list(title = "row z-score\nlog2(TPM+1)", direction = "horizontal", title_position = "topcenter", legend_width = unit(3, "cm"),
            title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7), at = seq(-1.5, 1.5, .5), labels = c("-1.5","-1","","0","","1","1.5"))

## ---- (2) gene blocks + gene table --------------------------------------------
genes6 <- list("Senescence/SASP"               = c("CDKN1A","CDKN2A","CDKN2B","TP53","GDF15","IL1A","CXCL2","SPP1","VEGFA"),
               "Stress / apoptosis\n(ISR-CHOP)" = c("DDIT3","TRIB3","ATF3","BCL2L11","PMAIP1","BBC3"),
               "Fibrosis"                      = c("TGFB2","LTBP1","ITGA11","COL15A1","P4HA2"),
               "Proliferation"                 = c("MKI67","TOP2A","CCNB1","CDK1","E2F1"))
cc6 <- setNames(c("#8C6BB1", "#B5533C", "#D08C3A", "#4D9F6B"), names(genes6))
g <- unlist(genes6, use.names = FALSE); cat_ <- factor(rep(names(genes6), lengths(genes6)), names(genes6))
Zc <- pmax(pmin(t(scale(t(log2(TPM[g, ] + 1)))), 1.5), -1.5)
lfc <- de$lfc[match(g, de$id)]; padj <- de$padj[match(g, de$id)]
star <- ifelse(padj < 1e-3, "***", ifelse(padj < 1e-2, "**", ifelse(padj < .05, "*", "ns")))
write.csv(data.frame(category = gsub("\n", " ", cat_), gene = g, log2FC_apeglm = round(lfc, 2), padj = signif(padj, 2), sig = star,
                     meanTPM_Control = round(rowMeans(TPM[g, 1:3]), 1), meanTPM_EGM = round(rowMeans(TPM[g, 4:6]), 1)),
          file.path(DAT, "heatmap_4blocks_genes.csv"), row.names = FALSE)

## ---- (3) final panel, EGM rows on top ----------------------------------------
g <- unlist(genes6, use.names = FALSE); cat_ <- factor(rep(names(genes6), lengths(genes6)), names(genes6))
Zc <- pmax(pmin(t(scale(t(log2(TPM[g, ] + 1)))), 1.5), -1.5); cs <- 5.2; GAP <- unit(1, "mm")
grp_top <- factor(as.character(grp), c("EGM", "Control"))          # EGM rows first
ord <- order(grp_top)                                               # EGM-1..3 then Control-1..3
c1 <- Heatmap(t(Zc)[ord, ], name = "z", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = cat_, row_split = grp_top[ord], column_title = NULL, row_title = NULL,
  top_annotation = HeatmapAnnotation(Category = anno_block(gp = gpar(fill = cc6, col = NA), labels = names(genes6), labels_gp = gpar(col = "white", fontsize = 8.5, fontface = "bold", lineheight = 0.85), height = unit(7, "mm")), show_annotation_name = FALSE),
  left_annotation = rowAnnotation(Group = anno_block(gp = gpar(fill = grpcol[c("EGM", "Control")], col = NA), labels = c("EGM", "Control"), labels_gp = gpar(col = "white", fontsize = 8, fontface = "bold"), labels_rot = 90, width = unit(4.5, "mm")), show_annotation_name = FALSE),
  column_names_gp = gpar(fontsize = 9, fontface = "italic"), column_names_rot = 45, row_names_gp = gpar(fontsize = 8), row_names_side = "right",
  rect_gp = gpar(col = "white", lwd = 0.8), column_gap = GAP, row_gap = GAP,
  width = unit(length(g) * cs, "mm"), height = unit(6 * cs, "mm"), heatmap_legend_param = leg)
row_labels <- c("Senescence/SASP", "Apoptosis\n(ISR-CHOP)", "Fibrosis", "Proliferation")
c2 <- Heatmap(Zc[, ord], name = "z", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = cat_, column_split = grp_top[ord], column_title = NULL, row_title = NULL,
  left_annotation = rowAnnotation(Category = anno_block(gp = gpar(fill = cc6, col = NA), labels = row_labels, labels_gp = gpar(col = "white", fontsize = 7.5, fontface = "bold", lineheight = 0.85), labels_rot = 90, width = unit(7, "mm")), show_annotation_name = FALSE),
  top_annotation = HeatmapAnnotation(Group = anno_block(gp = gpar(fill = grpcol[c("EGM", "Control")], col = NA), labels = c("EGM", "Control"), labels_gp = gpar(col = "white", fontsize = 8.5, fontface = "bold"), height = unit(4.5, "mm")), show_annotation_name = FALSE),
  row_names_gp = gpar(fontsize = 9, fontface = "italic"), row_names_side = "left", column_names_gp = gpar(fontsize = 8), column_names_rot = 45,
  rect_gp = gpar(col = "white", lwd = 0.8), row_gap = GAP, column_gap = GAP,
  width = unit(6 * 5.5, "mm"), height = unit(length(g) * cs, "mm"), heatmap_legend_param = leg)
png(file.path(FIG, "marker_heatmap_4blocks_EGMtop_genes_in_columns.png"), width = 7.2, height = 3.3, units = "in", res = 300); draw(c1, heatmap_legend_side = "bottom"); dev.off()
pdf(file.path(FIG, "marker_heatmap_4blocks_EGMtop_genes_in_columns.pdf"), width = 7.2, height = 3.3); draw(c1, heatmap_legend_side = "bottom"); dev.off()
png(file.path(FIG, "marker_heatmap_4blocks_EGMtop_genes_in_rows.png"), width = 3.4, height = 6.8, units = "in", res = 300); draw(c2, heatmap_legend_side = "bottom"); dev.off()
pdf(file.path(FIG, "marker_heatmap_4blocks_EGMtop_genes_in_rows.pdf"), width = 3.4, height = 6.8); draw(c2, heatmap_legend_side = "bottom"); dev.off()
cat("row order now:", paste(rownames(t(Zc)[ord, ]), collapse = ", "), "\n"); draw(c1, heatmap_legend_side = "bottom")
