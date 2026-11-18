# OncoKB core oncogenes in cynomolgus monkeys: EGM vs Control (Brain / Liver / Kidney / Spleen)
# Same approach as Mouse/oncogene_heatmap/oncogene_heatmap.R
# - Counts        : RSEM (tximport) gene counts, rounded; EGM (n=3) vs control (n=3) per organ
# - Human mapping : Ensembl BioMart M. fascicularis -> human orthologs (ensembl_mfascicularis_human_orthologs.csv),
#                   fallback = existing gene_symbol annotation
# - DE            : DESeq2, filter >=10 counts in >=3 samples
# - Enrichment    : % oncogenes up (padj<0.05) vs expression-matched non-oncogene background (10 bins, 2000 draws)
# - Core set      : OncoKB ONCOGENE listed in all 7 aggregated resources (n = 36);
#                   genes not expressed in >=1 organ are dropped from the heatmap

suppressMessages({
  library(DESeq2); library(dplyr); library(ComplexHeatmap); library(circlize)
})
set.seed(1)
root  <- "EGM_tissue"
kbase <- file.path(root, "oncogene_analysis/Monkey")
kout  <- file.path(kbase, "oncogene_heatmap")
ht_opt$message <- FALSE

# ---- oncogene lists & mapping ----
okb   <- read.delim(file.path(root, "cancerGeneList.tsv"), check.names = FALSE)
onc   <- subset(okb, `Gene Type` %in% c("ONCOGENE", "ONCOGENE_AND_TSG"))
core7 <- okb$`Hugo Symbol`[okb$`Gene Type` == "ONCOGENE" & okb$`# of occurrence within resources (Column K-P)` == 7]
onc_sym <- onc$`Hugo Symbol`
orth <- read.csv(file.path(kout, "ensembl_mfascicularis_human_orthologs.csv"), colClasses = "character")
korg <- c("Brain", "Kidney", "Liver", "Spleen")
ann  <- bind_rows(lapply(korg, function(t)
  read.csv(file.path(kbase, t, "DESeq2_results_annotated.csv"))[, c("gene_id", "gene_symbol")])) %>% distinct()
okm <- orth %>% filter(hs_symbol %in% onc_sym) %>% distinct(gene_id, hs_symbol, orth_type)
fbk <- ann %>% filter(gene_symbol %in% setdiff(onc_sym, okm$hs_symbol)) %>%
  transmute(gene_id, hs_symbol = gene_symbol, orth_type = "annotation_symbol")
omk <- bind_rows(okm, fbk) %>% left_join(onc %>% select(hs_symbol = `Hugo Symbol`, type = `Gene Type`), by = "hs_symbol")

# ---- DESeq2 per organ ----
res_k <- list(); vsd_k <- list(); shr_k <- list()
for (t in korg) {
  x <- read.delim(file.path(kbase, t, "txi.rsem.counts.txt"), row.names = 1, check.names = FALSE)
  cols <- c(paste0("control_", t, "_", 1:3), paste0("EGM_", t, "_", 1:3))
  meta <- data.frame(group = factor(rep(c("Control", "EGM"), each = 3), levels = c("Control", "EGM")), row.names = cols)
  dds <- DESeqDataSetFromMatrix(round(as.matrix(x[, cols])), meta, ~group)
  dds <- DESeq(dds[rowSums(counts(dds) >= 10) >= 3, ], quiet = TRUE)
  r <- as.data.frame(results(dds, contrast = c("group", "EGM", "Control"))); r$gene_id <- rownames(r)
  s <- as.data.frame(lfcShrink(dds, coef = "group_EGM_vs_Control", type = "apeglm", quiet = TRUE))
  res_k[[t]] <- r; vsd_k[[t]] <- vst(dds, blind = FALSE)
  shr_k[[t]] <- data.frame(gene_id = rownames(s), log2FC_shrunk = s$log2FoldChange)
}

# ---- oncogene results + expression-matched enrichment (all oncogenes and core set) ----
summ_k <- list(); onc_k <- list()
for (t in korg) {
  r <- res_k[[t]] %>% filter(!is.na(padj)) %>%
    mutate(onc = gene_id %in% omk$gene_id,
           dir = ifelse(padj < 0.05 & log2FoldChange > 0, "Up", ifelse(padj < 0.05 & log2FoldChange < 0, "Down", "NS")),
           bin = cut(log10(baseMean + 1), quantile(log10(baseMean + 1), 0:10 / 10), include.lowest = TRUE))
  ro <- r %>% filter(onc) %>% inner_join(omk, by = "gene_id") %>% group_by(hs_symbol) %>%
    slice_max(baseMean, n = 1, with_ties = FALSE) %>% ungroup()
  bg <- r %>% filter(!onc)
  permtest <- function(ids) {
    tab <- table(r$bin[match(ids, r$gene_id)])
    p <- replicate(2000, { s <- match(unlist(lapply(names(tab), function(b) sample(bg$gene_id[bg$bin == b], tab[b]))), r$gene_id)
                           mean(r$dir[s] == "Up") })
    obs <- mean(r$dir[match(ids, r$gene_id)] == "Up"); c(obs = 100 * obs, bg = 100 * mean(p), p = mean(p >= obs))
  }
  cc <- ro %>% filter(hs_symbol %in% core7); a <- permtest(ro$gene_id); b <- permtest(cc$gene_id)
  summ_k[[t]] <- data.frame(organ = t,
    all_n = nrow(ro), all_up = sum(ro$dir == "Up"), all_down = sum(ro$dir == "Down"),
    all_up_pct = round(a["obs"], 1), all_bg_up_pct = round(a["bg"], 1), all_perm_p = a["p"],
    core_n = nrow(cc), core_up = sum(cc$dir == "Up"), core_down = sum(cc$dir == "Down"), core_ns = sum(cc$dir == "NS"),
    core_up_pct = round(b["obs"], 1), core_bg_up_pct = round(b["bg"], 1), core_perm_p = b["p"], row.names = NULL)
  onc_k[[t]] <- ro %>% mutate(organ = t)
}
summ_k_df <- bind_rows(summ_k)
long_k <- bind_rows(lapply(korg, function(t) onc_k[[t]] %>%
  select(organ, hs_symbol, gene_id, type, baseMean, log2FoldChange, lfcSE, pvalue, padj, dir) %>%
  left_join(shr_k[[t]], by = "gene_id")))
core_k <- long_k %>% filter(hs_symbol %in% core7)
write.csv(omk, file.path(kout, "oncokb_oncogene_monkey_mapping.csv"), row.names = FALSE)
write.csv(long_k %>% arrange(organ, padj), file.path(kout, "oncogene_DE_EGM_vs_Control_by_organ.csv"), row.names = FALSE)
write.csv(core_k %>% arrange(hs_symbol, organ), file.path(kout, "core_oncogenes_DE_by_organ.csv"), row.names = FALSE)
write.csv(summ_k_df, file.path(kout, "oncogene_enrichment_summary.csv"), row.names = FALSE)

# ---- core oncogene heatmap (same style as mouse figure) ----
org4 <- c("Brain", "Liver", "Kidney", "Spleen")   # panel order
core_hs <- sort(core7)
R <- do.call(cbind, lapply(org4, function(t) {     # VST minus control mean, within organ
  v <- assay(vsd_k[[t]]); colnames(v) <- sub(paste0("_", t), "", colnames(v))
  v <- v[, c(paste0("control_", 1:3), paste0("EGM_", 1:3))]
  d <- core_k %>% filter(organ == t)
  m <- matrix(NA, length(core_hs), 6, dimnames = list(core_hs, paste(t, colnames(v), sep = "_")))
  vv <- v[d$gene_id, , drop = FALSE]; m[d$hs_symbol, ] <- vv - rowMeans(vv[, 1:3]); m
}))
R <- R[rowSums(is.na(R)) == 0, ]                  # keep genes expressed in all organs
lab <- sapply(org4, function(t) { d <- table(factor((core_k %>% filter(organ == t, hs_symbol %in% rownames(R)))$dir, levels = c("Up", "NS", "Down")))
  sprintf("Up %d · Down %d\nNS %d", d["Up"], d["Down"], d["NS"]) })
grp <- factor(rep(rep(c("Control", "EGM"), each = 3), length(org4)), levels = c("Control", "EGM"))
vir <- colorRamp2(seq(-2, 2, length.out = 256), viridisLite::viridis(256))
hK <- Heatmap(R, name = "log2 ratio\nvs. Control mean", col = vir,
  cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = factor(rep(org4, each = 6), levels = org4), column_gap = unit(2, "mm"),
  top_annotation = HeatmapAnnotation(Group = grp, col = list(Group = c(Control = "#d1d2d7", EGM = "#6d98f8")),
    show_annotation_name = FALSE, simple_anno_size = unit(3, "mm"),
    annotation_legend_param = list(Group = list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)))),
  bottom_annotation = HeatmapAnnotation(summary = anno_block(labels = lab, labels_gp = gpar(fontsize = 7),
    gp = gpar(fill = NA, col = NA), height = unit(8, "mm")), show_annotation_name = FALSE),
  show_column_names = FALSE, row_names_side = "right", row_names_gp = gpar(fontsize = 9.5),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  width = unit(length(org4) * 6 * 4.2, "mm"), height = unit(nrow(R) * 3.6, "mm"),
  rect_gp = gpar(col = "white", lwd = 0.4),
  heatmap_legend_param = list(at = -2:2, title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)))
f <- file.path(kout, "Fig_core_oncogenes_heatmap_monkey")
for (dev in c("pdf", "png")) {
  if (dev == "pdf") pdf(paste0(f, ".pdf"), width = 6.2, height = 6.8) else png(paste0(f, ".png"), width = 6.2, height = 6.8, units = "in", res = 300)
  draw(hK, column_title = sprintf("Core oncogenes (OncoKB consensus, n = %d)\nCynomolgus monkeys: EGM vs Control", nrow(R)),
       column_title_gp = gpar(fontsize = 12), merge_legend = TRUE, padding = unit(c(2, 4, 2, 6), "mm"))
  dev.off()
}
print(summ_k_df)
