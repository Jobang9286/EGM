## ---- setup, KEGG over-representation ----------------------------------------
suppressMessages({
  library(DESeq2); library(clusterProfiler); library(org.Mm.eg.db)
})

obj <- readRDS("analysis/deseq2_objects.rds")
dds <- obj$dds
sig <- obj$sig

cat("dds samples:", paste(colnames(dds), collapse=", "), "\n")
cat("sig rows:", nrow(sig), " padj<0.05:", sum(sig$padj < 0.05, na.rm=TRUE), "\n")

map <- bitr(rownames(dds), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
cat("mapped genes:", nrow(map), "of", nrow(dds), "\n")

sig_genes  <- sig$id[!is.na(sig$padj) & sig$padj < 0.05]
sig_entrez <- na.omit(map$ENTREZID[match(sig_genes, map$SYMBOL)])
universe   <- na.omit(map$ENTREZID)

ekegg <- enrichKEGG(gene = sig_entrez, organism = "mmu", universe = universe,
                     pvalueCutoff = 1, qvalueCutoff = 1)
ek <- as.data.frame(ekegg)
targets <- c("FoxO signaling pathway","Signaling pathways regulating pluripotency of stem cells",
             "Cell cycle","IL-17 signaling pathway","TNF signaling pathway","Pathways in cancer")
hit <- ek[ek$Description %in% targets, c("Description","GeneRatio","pvalue","p.adjust","Count")]
print(hit)

## ---- vst matrix (the 25-gene panelA/panelB built here are a first
## ---- layout, overwritten below; kept because `grp` is taken from their columns) -
suppressMessages(library(ComplexHeatmap))

get_pathway_genes <- function(kegg_id) {
  entrez <- strsplit(ek$geneID[ek$ID == kegg_id], "/")[[1]]
  sym <- map$SYMBOL[match(entrez, map$ENTREZID)]
  keep <- !is.na(sym)
  sym <- sym[keep]
  lfc <- sig$lfc[match(sym, sig$id)]
  data.frame(gene = sym, lfc = lfc)[order(-abs(lfc)), ]
}

vsd <- DESeq2::vst(dds, blind = FALSE)
mat_all <- SummarizedExperiment::assay(vsd)
samp_order <- c("CP-1","CP-2","CP-3","4FP-1","4FP-2","4FP-3")
mat_all <- mat_all[, samp_order]

CAP <- 25
build_block <- function(name, kegg_id, cap = CAP) {
  df <- get_pathway_genes(kegg_id)
  df <- df[df$gene %in% rownames(mat_all), ]
  df <- head(df, cap)
  data.frame(pathway = name, gene = df$gene, lfc = df$lfc, stringsAsFactors = FALSE)
}

panelA <- rbind(
  build_block("FoxO signaling pathway", "mmu04068"),
  build_block("Pluripotency of stem cells", "mmu04550"),
  build_block("Cell cycle", "mmu04110")
)
panelB <- rbind(
  build_block("IL-17 signaling pathway", "mmu04657"),
  build_block("TNF signaling pathway", "mmu04668"),
  build_block("Pathways in cancer", "mmu05200")
)

## ---- column groups ----------------------------------------------------------
suppressMessages({library(ComplexHeatmap); library(circlize); library(grid)})

zmat <- function(panel) {
  m <- mat_all[panel$gene, , drop = FALSE]
  z <- t(scale(t(m)))
  rownames(z) <- make.unique(panel$gene)
  z
}
zA <- zmat(panelA)

grp <- factor(ifelse(grepl("^CP", colnames(zA)), "Ctrl-BP", "P-BP"), levels = c("Ctrl-BP","P-BP"))

## ---- direction-aware gene picker --------------------------------------------
get_pathway_genes_dir <- function(kegg_id, direction = c("up","down")) {
  direction <- match.arg(direction)
  entrez <- strsplit(ek$geneID[ek$ID == kegg_id], "/")[[1]]
  sym <- map$SYMBOL[match(entrez, map$ENTREZID)]
  sym <- sym[!is.na(sym) & sym %in% rownames(mat_all)]
  lfc <- sig$lfc[match(sym, sig$id)]
  df  <- data.frame(gene = sym, lfc = lfc)
  df  <- if (direction == "up") df[df$lfc > 0, ] else df[df$lfc < 0, ]
  df[order(-abs(df$lfc)), ]
}

## ---- column annotation ------------------------------------------------------
top_ann <- HeatmapAnnotation(Group = grp,
  col = list(Group = c("Ctrl-BP" = "grey65", "P-BP" = "#3B7DD8")),
  show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(2.5, "mm"))

## ---- colour scale -1.5 .. 1.5 -----------------------------------------------
col_fun <- colorRamp2(seq(-1.5, 1.5, length.out = 7), rev(RColorBrewer::brewer.pal(7, "RdYlBu")))

## ---- pathway specs ----------------------------------------------------------
specsA <- list(
  "FoxO signaling pathway"      = list(id = "mmu04068", dir = "up"),
  "Pluripotency of stem cells"  = list(id = "mmu04550", dir = "up"),
  "Cell cycle"                  = list(id = "mmu04110", dir = "up")
)
specsB <- list(
  "IL-17 signaling pathway" = list(id = "mmu04657", dir = "down"),
  "TNF signaling pathway"   = list(id = "mmu04668", dir = "down"),
  "Pathways in cancer"      = list(id = "mmu05200", dir = "down")
)

## ---- gene selection, no duplicates across blocks, panel A -------------------
build_panel_dedup_v2 <- function(spec_list, caps) {
  used <- character(0)
  out <- list()
  for (nm in names(spec_list)) {
    id  <- spec_list[[nm]]$id
    dir <- spec_list[[nm]]$dir
    cap <- caps[[nm]]
    df <- get_pathway_genes_dir(id, dir)
    df <- df[!(df$gene %in% used), ]
    df <- head(df, cap)
    used <- c(used, df$gene)
    out[[nm]] <- data.frame(pathway = nm, gene = df$gene, lfc = df$lfc, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

capsA <- list("FoxO signaling pathway" = 14, "Pluripotency of stem cells" = 14, "Cell cycle" = 10)
panelA <- build_panel_dedup_v2(specsA, capsA)
cat("panelA total:", nrow(panelA), "\n"); print(table(panelA$pathway))

## ---- panel B and z-scores ---------------------------------------------------
capsB <- list("IL-17 signaling pathway" = 10, "TNF signaling pathway" = 9, "Pathways in cancer" = 19)
panelB <- build_panel_dedup_v2(specsB, capsB)
cat("panelB total:", nrow(panelB), "\n"); print(table(panelB$pathway))

zA <- mat_all[panelA$gene, , drop = FALSE]; zA <- t(scale(t(zA))); rownames(zA) <- panelA$gene
zB <- mat_all[panelB$gene, , drop = FALSE]; zB <- t(scale(t(zB))); rownames(zB) <- panelB$gene
row_split_A <- factor(panelA$pathway, levels = unique(panelA$pathway))
row_split_B <- factor(panelB$pathway, levels = unique(panelB$pathway))

## ---- heatmap builder with a left block-colour bar ---------------------------
mk_final_bar <- function(z, split, bar_color) {
  left_ann <- rowAnnotation(
    block = anno_block(gp = gpar(fill = bar_color, col = NA), width = unit(2.5, "mm")),
    show_annotation_name = FALSE
  )
  Heatmap(z, name = "z-score", col = col_fun, top_annotation = top_ann, left_annotation = left_ann,
    cluster_rows = TRUE, cluster_row_slices = FALSE, show_row_dend = FALSE,
    cluster_columns = FALSE, show_column_names = FALSE,
    row_labels = gsub("\\.[0-9]+$", "", rownames(z)),
    column_split = grp, column_title = c("Ctrl-BP","P-BP"), column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    column_gap = unit(0, "mm"),
    row_split = split, row_title = NULL,
    row_gap = unit(0, "mm"),
    border = FALSE,
    row_names_gp = gpar(fontsize = 8), row_names_side = "right",
    heatmap_legend_param = list(title = "z-score", legend_height = unit(2.5, "cm"),
                                 at = c(-1.5,-1,-0.5,0,0.5,1,1.5)))
}

## ---- writer -----------------------------------------------------------------
dir.create("analysis/figures/pathway_heatmap_redraw", showWarnings = FALSE)
save_seamless <- function(fname_base, ht) {
  png(file.path("analysis/figures/pathway_heatmap_redraw", paste0(fname_base, ".png")),
      width = 6.5, height = 5.2, units = "in", res = 400)
  draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
  dev.off()
  pdf(file.path("analysis/figures/pathway_heatmap_redraw", paste0(fname_base, ".pdf")),
      width = 6.5, height = 5.2)
  draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
  dev.off()
}

## ---- final block colours ----------------------------------------------------
suppressMessages({library(ComplexHeatmap); library(circlize); library(grid)})

colsA <- c("#FFB703", "#FB8500", "#E63946")   # FoxO, Pluripotency, Cell cycle - vivid amber/orange/red
colsB <- c("#48CAE4", "#00B4D8", "#0077B6")   # IL-17, TNF, Cancer - vivid cyan/blue

ht_up   <- mk_final_bar(zA, row_split_A, colsA)
ht_down <- mk_final_bar(zB, row_split_B, colsB)

## ---- final save -------------------------------------------------------------
save_seamless("panelA_up_foxo_pluripotency_cellcycle", ht_up)
save_seamless("panelB_down_il17_tnf_cancer", ht_down)
