# OncoKB oncogenes in HGPS mice: EGM vs Control (Brain / Kidney / Liver)
# - Oncogene list : OncoKB cancerGeneList.tsv, Gene Type ONCOGENE + ONCOGENE_AND_TSG
# - Human -> mouse: Ensembl BioMart orthologs (saved in ensembl_mouse_human_orthologs.csv),
#                   fallback = case-converted symbol match
# - DE            : DESeq2 per organ, HGPS_EGM (n=3) vs HGPS_control (n=3); WT aging samples excluded
#                   filter: >=10 counts in >=3 samples; display LFC = apeglm shrinkage
# - Enrichment    : % oncogenes up (padj<0.05) vs expression-matched non-oncogene background
#                   (10 baseMean bins, 2000 draws); one-sided permutation p

suppressMessages({
  library(DESeq2); library(dplyr); library(ComplexHeatmap); library(circlize)
})
set.seed(1)
root   <- "EGM_tissue"
mbase  <- file.path(root, "oncogene_analysis/Mouse")
outdir <- file.path(mbase, "oncogene_heatmap")
ht_opt$message <- FALSE

# ---- oncogene list & mapping ----
okb <- read.delim(file.path(root, "cancerGeneList.tsv"), check.names = FALSE)
onc <- subset(okb, `Gene Type` %in% c("ONCOGENE", "ONCOGENE_AND_TSG"))
orth_mm <- read.csv(file.path(outdir, "ensembl_mouse_human_orthologs.csv"), colClasses = "character")

files <- c(Brain = "brain/counts_brain.txt", Kidney = "kidney/counts_kidney.txt", Liver = "liver/counts_liver.csv")
cnt_mm <- lapply(files, function(f) {
  x <- read.table(file.path(mbase, f), header = TRUE, sep = ifelse(grepl("csv$", f), ",", "\t"),
                  check.names = FALSE, row.names = 1, comment.char = "")
  colnames(x) <- sub("\\.bam$", "", colnames(x)); x
})
genes_mm <- unique(unlist(lapply(cnt_mm, rownames)))
onc_sym  <- onc$`Hugo Symbol`
om <- orth_mm %>% filter(hs_symbol %in% onc_sym, mm_symbol %in% genes_mm) %>% distinct(mm_symbol, hs_symbol, orth_type)
fb <- tibble(hs_symbol = setdiff(onc_sym, om$hs_symbol)) %>%
  mutate(mm_symbol = paste0(substr(hs_symbol, 1, 1), tolower(substring(hs_symbol, 2)))) %>%
  filter(mm_symbol %in% genes_mm) %>% mutate(orth_type = "case_match")
om2 <- bind_rows(om, fb) %>% left_join(onc %>% select(hs_symbol = `Hugo Symbol`, type = `Gene Type`), by = "hs_symbol")

# ---- DESeq2 per organ ----
dds_m <- list(); res_m <- list(); vsd_m <- list(); shr <- list()
for (t in names(cnt_mm)) {
  x <- cnt_mm[[t]]
  cols <- grep("_HGPS_(control|EGM)_", colnames(x), value = TRUE)
  meta <- data.frame(group = factor(ifelse(grepl("_EGM_", cols), "EGM", "Control"), levels = c("Control", "EGM")), row.names = cols)
  dds <- DESeqDataSetFromMatrix(round(as.matrix(x[, cols])), meta, ~group)
  dds <- DESeq(dds[rowSums(counts(dds) >= 10) >= 3, ], quiet = TRUE)
  r <- as.data.frame(results(dds, contrast = c("group", "EGM", "Control"))); r$mm_symbol <- rownames(r)
  s <- as.data.frame(lfcShrink(dds, coef = "group_EGM_vs_Control", type = "apeglm", quiet = TRUE))
  dds_m[[t]] <- dds; res_m[[t]] <- r; vsd_m[[t]] <- vst(dds, blind = FALSE)
  shr[[t]] <- data.frame(mm_symbol = rownames(s), log2FC_shrunk = s$log2FoldChange)
}

# ---- oncogene results + expression-matched enrichment ----
summ <- list(); onc_m <- list()
for (t in names(res_m)) {
  r <- res_m[[t]] %>% filter(!is.na(padj)) %>%
    mutate(onc = mm_symbol %in% om2$mm_symbol,
           dir = ifelse(padj < 0.05 & log2FoldChange > 0, "Up", ifelse(padj < 0.05 & log2FoldChange < 0, "Down", "NS")),
           bin = cut(log10(baseMean + 1), quantile(log10(baseMean + 1), 0:10 / 10), include.lowest = TRUE))
  # one mouse gene per human oncogene: identical symbol first, else highest baseMean
  ro <- r %>% filter(onc) %>% inner_join(om2, by = "mm_symbol") %>%
    mutate(same = toupper(mm_symbol) == hs_symbol) %>% group_by(hs_symbol) %>%
    arrange(desc(same), desc(baseMean), .by_group = TRUE) %>% slice(1) %>% ungroup()
  bg <- r %>% filter(!onc); tab <- table(ro$bin)
  samp <- function() unlist(lapply(names(tab), function(b) sample(bg$mm_symbol[bg$bin == b], tab[b])))
  perm <- replicate(2000, { s <- match(samp(), r$mm_symbol); c(up = mean(r$dir[s] == "Up"), med = median(r$log2FoldChange[s])) })
  summ[[t]] <- data.frame(organ = t, onc_expressed = nrow(ro), onc_up = sum(ro$dir == "Up"), onc_down = sum(ro$dir == "Down"),
    onc_ns = sum(ro$dir == "NS"), onc_up_pct = round(100 * mean(ro$dir == "Up"), 1),
    matched_bg_up_pct = round(100 * mean(perm["up", ]), 1), perm_p_up = mean(perm["up", ] >= mean(ro$dir == "Up")),
    onc_down_pct = round(100 * mean(ro$dir == "Down"), 1),
    onc_medLFC = round(median(ro$log2FoldChange), 3), matched_bg_medLFC = round(mean(perm["med", ]), 3),
    perm_p_medLFC = mean(perm["med", ] >= median(ro$log2FoldChange)),
    onc_up_2fold = sum(ro$padj < 0.05 & ro$log2FoldChange > 1),
    bg_up2_pct = round(100 * mean(bg$padj < 0.05 & bg$log2FoldChange > 1), 1))
  onc_m[[t]] <- ro %>% mutate(organ = t)
}
summ_df <- bind_rows(summ)
long <- bind_rows(lapply(names(onc_m), function(t) onc_m[[t]] %>%
  select(organ, hs_symbol, mm_symbol, type, baseMean, log2FoldChange, lfcSE, pvalue, padj, dir) %>%
  left_join(shr[[t]], by = "mm_symbol")))
write.csv(om2, file.path(outdir, "oncokb_oncogene_mouse_mapping.csv"), row.names = FALSE)
write.csv(long %>% arrange(organ, padj), file.path(outdir, "oncogene_DE_EGM_vs_Control_by_organ.csv"), row.names = FALSE)
write.csv(summ_df, file.path(outdir, "oncogene_enrichment_summary.csv"), row.names = FALSE)

# ---- shared styling ----
organs  <- names(vsd_m)
hs_rows <- sort(unique(long$hs_symbol))
div     <- colorRamp2(c(-2, 0, 2), c("#1c5cab", "#f0efec", "#e34948"))
col_grp <- c(Control = "#d1d2d7", EGM = "#6d98f8")
save_both <- function(ht, f, w, h, title, ...) for (dev in c("pdf", "png")) {
  if (dev == "pdf") pdf(paste0(f, ".pdf"), width = w, height = h) else png(paste0(f, ".png"), width = w, height = h, units = "in", res = 300)
  draw(ht, column_title = title, column_title_gp = gpar(fontsize = 12), ...); dev.off()
}

# ---- Fig 1: sample-level z-score (genes expressed in all organs) ----
col_order <- c(paste0("control_", 1:3), paste0("EGM_", 1:3))
Z <- do.call(cbind, lapply(organs, function(t) {
  v <- assay(vsd_m[[t]]); colnames(v) <- sub(".*_HGPS_", "", colnames(v)); v <- v[, col_order]
  ids <- onc_m[[t]] %>% select(hs_symbol, mm_symbol)
  z <- t(scale(t(v[ids$mm_symbol, , drop = FALSE]))); rownames(z) <- ids$hs_symbol
  m <- matrix(NA, length(hs_rows), 6, dimnames = list(hs_rows, paste(t, col_order, sep = "_"))); m[rownames(z), ] <- z; m
}))
Zc <- Z[rowSums(is.na(Z)) == 0, ]
grp <- factor(rep(rep(c("Control", "EGM"), each = 3), 3), levels = c("Control", "EGM"))
org <- factor(rep(organs, each = 6), levels = organs)
lab1 <- sapply(organs, function(t) { d <- table(factor(onc_m[[t]]$dir[onc_m[[t]]$hs_symbol %in% rownames(Zc)], levels = c("Up", "NS", "Down")))
  sprintf("Up %d · NS %d · Down %d", d["Up"], d["NS"], d["Down"]) })
h1 <- Heatmap(Zc, name = "Row z-score\n(VST, within organ)", col = div,
  cluster_rows = hclust(dist(Zc), method = "ward.D2"), show_row_dend = FALSE, cluster_columns = FALSE,
  column_split = org, column_gap = unit(2.5, "mm"), show_row_names = FALSE, show_column_names = FALSE,
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  top_annotation = HeatmapAnnotation(Group = grp, col = list(Group = col_grp), annotation_name_side = "left", simple_anno_size = unit(3.5, "mm")),
  bottom_annotation = HeatmapAnnotation(summary = anno_block(labels = lab1, labels_gp = gpar(fontsize = 7.5), gp = gpar(fill = NA, col = NA)), show_annotation_name = FALSE),
  use_raster = FALSE, heatmap_legend_param = list(at = -2:2))
save_both(h1, file.path(outdir, "Fig_oncogene_heatmap_zscore"), 6.5, 8,
  sprintf("OncoKB oncogenes in HGPS mice: EGM vs Control\n(%d genes expressed in all three organs)", nrow(Zc)), merge_legend = TRUE)

# ---- Fig 2: log2FC (apeglm), non-significant = 0, not expressed = white ----
L <- P <- matrix(NA, length(hs_rows), length(organs), dimnames = list(hs_rows, organs))
for (t in organs) { d <- long %>% filter(organ == t); L[d$hs_symbol, t] <- d$log2FC_shrunk; P[d$hs_symbol, t] <- d$padj }
L[!is.na(P) & P >= 0.05] <- 0
ord <- order(rowMeans(L, na.rm = TRUE), decreasing = TRUE); L <- L[ord, ]
mm_lab  <- (long %>% distinct(hs_symbol, mm_symbol) %>% group_by(hs_symbol) %>% slice(1))$mm_symbol
names(mm_lab) <- (long %>% distinct(hs_symbol, mm_symbol) %>% group_by(hs_symbol) %>% slice(1))$hs_symbol
lab_idx <- which(apply(L, 1, function(x) any(x > 1, na.rm = TRUE)))
lab2 <- sapply(organs, function(t) { d <- table(factor(onc_m[[t]]$dir, levels = c("Up", "NS", "Down")))
  sprintf("Up %d\nNS %d\nDown %d", d["Up"], d["NS"], d["Down"]) })
h2 <- Heatmap(L, name = "log2FC\nEGM / Control\n(padj<0.05 only)", col = div, na_col = "white",
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, width = unit(36, "mm"),
  column_names_rot = 0, column_names_centered = TRUE, column_names_side = "top", column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  bottom_annotation = HeatmapAnnotation(n = anno_text(lab2, rot = 0, just = "center", location = 0.5, gp = gpar(fontsize = 7.5), height = unit(14, "mm")), show_annotation_name = FALSE),
  right_annotation = rowAnnotation(lab = anno_mark(at = lab_idx, labels = mm_lab[rownames(L)][lab_idx], labels_gp = gpar(fontsize = 6.5), link_width = unit(4, "mm"))),
  heatmap_legend_param = list(at = -2:2), use_raster = FALSE)
save_both(h2, file.path(outdir, "Fig_oncogene_heatmap_log2FC"), 4, 8,
  sprintf("OncoKB oncogenes (%d)\nHGPS mice, EGM vs Control", nrow(L)))

print(summ_df)

# ---- Fig 3: core oncogenes only (pre-specified, result-independent criterion) ----
# ONCOGENE listed in all 7 OncoKB-aggregated resources
# (OncoKB, MSK-IMPACT, MSK-HEME, FoundationOne, FoundationOne Heme, Vogelstein, COSMIC CGC)
core7 <- okb$`Hugo Symbol`[okb$`Gene Type` == "ONCOGENE" & okb$`# of occurrence within resources (Column K-P)` == 7]
core  <- long %>% filter(hs_symbol %in% core7)
write.csv(core %>% arrange(hs_symbol, organ), file.path(outdir, "core_oncogenes_DE_by_organ.csv"), row.names = FALSE)
core_hs <- sort(core7)
org3 <- c("Brain", "Liver", "Kidney")   # panel order
mm_of <- sapply(core_hs, function(h) { x <- unique(core$mm_symbol[core$hs_symbol == h]); if (length(x)) x[1] else NA })
R <- do.call(cbind, lapply(org3, function(t) {           # VST minus control mean, within organ
  v <- assay(vsd_m[[t]]); colnames(v) <- sub(".*_HGPS_", "", colnames(v)); v <- v[, col_order]
  m <- matrix(NA, length(core_hs), 6, dimnames = list(core_hs, paste(t, col_order, sep = "_")))
  ok <- !is.na(mm_of) & core_hs %in% core$hs_symbol[core$organ == t]
  vv <- v[mm_of[ok], , drop = FALSE]; m[ok, ] <- vv - rowMeans(vv[, 1:3]); m
}))
rownames(R) <- mm_of
R <- R[rowSums(is.na(R)) == 0, ]   # keep genes expressed in all three organs
lab3 <- sapply(org3, function(t) { d <- table(factor((core %>% filter(organ == t, mm_symbol %in% rownames(R)))$dir, levels = c("Up", "NS", "Down")))
  sprintf("Up %d · Down %d\nNS %d", d["Up"], d["Down"], d["NS"]) })
vir <- colorRamp2(seq(-2, 2, length.out = 256), viridisLite::viridis(256))
hA <- Heatmap(R, name = "log2 ratio\nvs. Control mean", col = vir, na_col = "white",
  cluster_rows = FALSE, cluster_columns = FALSE, column_split = factor(rep(org3, each = 6), levels = org3), column_gap = unit(2, "mm"),
  top_annotation = HeatmapAnnotation(Group = grp, col = list(Group = col_grp), show_annotation_name = FALSE,
    simple_anno_size = unit(3, "mm"), annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = list(Group = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)))),
  bottom_annotation = HeatmapAnnotation(summary = anno_block(labels = lab3, labels_gp = gpar(fontsize = 7),
    gp = gpar(fill = NA, col = NA), height = unit(8, "mm")), show_annotation_name = FALSE),
  show_column_names = FALSE, row_names_side = "right", row_names_gp = gpar(fontsize = 9.5),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"), width = unit(18 * 4.2, "mm"), height = unit(nrow(R) * 3.6, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.4),
  heatmap_legend_param = list(at = -2:2, title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)))
save_both(hA, file.path(outdir, "Fig_core_oncogenes_heatmap"), 5.0, 6.8,
  sprintf("Core oncogenes (OncoKB consensus, n = %d)\nHGPS mice: EGM vs Control", nrow(R)),
  merge_legend = TRUE, padding = unit(c(2, 4, 2, 6), "mm"))
