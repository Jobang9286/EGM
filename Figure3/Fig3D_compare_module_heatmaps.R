# =============================================================================
# EGM hdWGCNA module heatmap: current modules vs. the corrected merge, drawn with identical code
#
# The May 2026 hdWGCNA had 13 modules before merging. py_hdWGCNA's merge looked up
# module eigengenes one row off (grey counted as a row), so it merged M4 -> M3 (ME r 0.63) and M12 -> M2
# (r -0.01) instead of the pairs above the 0.2 cut, M9 -> M5 (r 0.88) and M13 -> M6 (r 0.92).
#   current   = May's final modules (M1-M3, M5-M11, M13)
#   corrected = May's pre-merge modules with M9 -> M5 and M13 -> M6 (M1-M8, M10-M12)
#   unmerged  = all 13 pre-merge modules, no merging (M1-M13)
# Cell hMEs for both were recomputed from May's saved object with the same function and Harmony
# settings (corrected_hmes.py); its output, won_hpc_corrected_hmes_6505649.out.txt,
# holds the per-sample mean hME, the Mann-Whitney p vs A1 and the pre-merge M4 / M12 gene lists.
# z = row z-score of the sample means, stars = BH padj x |z| cut-offs,
# hallmarks = one-sided Fisher vs msigdbr sets (same msigdbr for both).
# Output: heatmap_current.pdf/png, heatmap_corrected.pdf/png, comparison tables (*.csv)
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(msigdbr); library(RColorBrewer)
  library(ComplexHeatmap); library(circlize); library(grid)
})
ROOT <- "4FP_plasma/EGM_hdWGCNA_code"
DIR <- file.path(ROOT, "03_figures_R/output")
ORG <- file.path(ROOT, "data/may2026_11modules")         # = egm_final/for_rstudio
L   <- unlist(lapply(c("won_hpc_corrected_hmes_6505649.out.txt", "won_hpc_unmerged_hmes_6505654.out.txt"),
                     function(f) readLines(file.path(ROOT, "02_13modules_won_hpc/output", f))))
blk <- function(title) {
  i <- grep(title, L, fixed = TRUE); h <- grep("^###", L); j <- min(c(h[h > i], length(L) + 1))
  read.csv(text = L[(i + 1):(j - 1)], row.names = 1, check.names = FALSE)
}
premerge <- function(m) strsplit(sub(".*: ", "", grep(sprintf("^### premerge %s genes:", m), L, value = TRUE)), ",")[[1]]

# ---- module assignments -----------------------------------------------------
ma  <- read.csv(file.path(ORG, "all_modules_assignments.csv"), check.names = FALSE)
bg  <- unique(ma$gene_name)
cur <- setNames(ma$module, ma$gene_name)
cor_ <- cur
cor_[cur == "M9"]  <- "M5"; cor_[cur == "M13"] <- "M6"
stopifnot(all(cur[premerge("M4")] == "M3"), all(cur[premerge("M12")] == "M2"))
cor_[premerge("M4")] <- "M4"; cor_[premerge("M12")] <- "M12"
unm <- cur; unm[premerge("M4")] <- "M4"; unm[premerge("M12")] <- "M12"
assign_ <- list(current = cur, corrected = cor_, unmerged = unm)

# ---- gene sets ----------------------------------------------------------------
msig <- msigdbr(species = "Mus musculus")
hm_aging <- list(
  Senescence = c("REACTOME_CELLULAR_SENESCENCE","GOBP_CELLULAR_SENESCENCE","FRIDMAN_SENESCENCE_UP","SAUL_SEN_MAYO"),
  SASP       = c("REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
                 "WP_SENESCENCEASSOCIATED_SECRETORY_PHENOTYPE_SASP",
                 "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_INFLAMMATORY_RESPONSE"),
  ECM        = c("NABA_MATRISOME","GOBP_EXTRACELLULAR_MATRIX_ORGANIZATION",
                 "GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION"),   # GO:0045229, the set ECM organization was merged into
  IFN        = c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE"),
  `DNA repair`= c("HALLMARK_DNA_REPAIR","GOBP_DNA_REPAIR"),
  Telomere   = c("GOBP_TELOMERE_ORGANIZATION","REACTOME_EXTENSION_OF_TELOMERES",
                 "GOBP_POSITIVE_REGULATION_OF_TELOMERE_MAINTENANCE"),
  Chromatin  = c("GOBP_CHROMATIN_ORGANIZATION","GOBP_CHROMATIN_REMODELING"),
  UPR        = c("HALLMARK_UNFOLDED_PROTEIN_RESPONSE","GOBP_CELLULAR_RESPONSE_TO_UNFOLDED_PROTEIN"),
  Autophagy  = c("GOBP_MACROAUTOPHAGY"),
  mTOR       = c("HALLMARK_MTORC1_SIGNALING","HALLMARK_PI3K_AKT_MTOR_SIGNALING"),
  OxPhos     = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
  `Stem cell exh.` = c("HALLMARK_WNT_BETA_CATENIN_SIGNALING","HALLMARK_NOTCH_SIGNALING","HALLMARK_HEDGEHOG_SIGNALING",
                       "GOBP_MESENCHYMAL_STEM_CELL_PROLIFERATION","GOBP_STEM_CELL_DIVISION"))
hm_id <- list(`ESC core` = c("WONG_EMBRYONIC_STEM_CELL_CORE","BHATTACHARYA_EMBRYONIC_STEM_CELL",
                             "BENPORATH_ES_CORE_NINE_CORRELATED"),
              `PRC2-repressed somatic` = c("BENPORATH_PRC2_TARGETS","BENPORATH_SUZ12_TARGETS",
                                           "MIKKELSEN_MEF_HCP_WITH_H3K27ME3","MIKKELSEN_DEDIFFERENTIATED_STATE_DN"))
sets_of <- function(Lst) lapply(Lst, function(s) intersect(msig %>% filter(gs_name %in% s) %>% pull(gene_symbol) %>% unique(), bg))
S_aging <- sets_of(hm_aging); S_id <- sets_of(hm_id)

enrich <- function(mods) {   # one-sided Fisher, background = the 9,559 network genes; BH over both blocks
  N <- length(bg)
  r <- do.call(rbind, lapply(names(mods), function(m) do.call(rbind, lapply(c(names(S_aging), names(S_id)), function(h) {
    sg <- c(S_aging, S_id)[[h]]; mg <- intersect(mods[[m]], bg)
    a <- length(intersect(mg, sg)); b <- length(mg) - a; cc <- length(sg) - a; d <- N - a - b - cc
    data.frame(module = m, hallmark = h, block = if (h %in% names(S_aging)) "aging" else "id",
               p = fisher.test(matrix(c(a, cc, b, d), 2), alternative = "greater")$p.value)
  }))))
  r$padj <- p.adjust(r$p, "BH"); r$nlogp <- pmin(-log10(r$padj), 10)
  r$lab <- ifelse(r$padj < 1e-4, "***", ifelse(r$padj < 1e-2, "**", ifelse(r$padj < 0.05, "*", "")))
  r
}

samples <- c("A1","Y1","R0","R1","R2","R3","R4","R5","R6","R7","R8","R9")
mk_star <- function(p, z) {
  if (is.na(p) || is.na(z)) return("")
  if (p < 1e-20 && abs(z) > 1.2) return("***"); if (p < 1e-5 && abs(z) > 0.8) return("**")
  if (p < 1e-2 && abs(z) > 0.3) return("*"); ""
}
build <- function(tag) {
  mean_ <- blk(sprintf("### %s mean hME", tag)); pm <- as.matrix(blk(sprintf("### %s Mann-Whitney p", tag)))
  mods_order <- colnames(mean_)[order(as.integer(sub("M", "", colnames(mean_))))]
  z <- t(scale(as.matrix(mean_[samples, mods_order])))                       # row z-score across the 12 samples
  pm <- pm[mods_order, samples]
  padj <- matrix(p.adjust(pm, "BH"), nrow(pm), dimnames = dimnames(pm))
  sig_z <- matrix(mapply(mk_star, padj, z), nrow(z), dimnames = dimnames(z))
  a <- assign_[[tag]]; mods <- split(names(a), a); mods <- mods[mods_order]
  e <- enrich(mods)
  wide <- function(col) { w <- e %>% filter(block == "aging") %>% select(module, hallmark, all_of(col)) %>%
    pivot_wider(names_from = hallmark, values_from = all_of(col)) %>% as.data.frame()
    rownames(w) <- w$module; as.matrix(w[mods_order, names(hm_aging)]) }
  list(z = z, sig_z = sig_z, mat_a = wide("nlogp"), sig_a = wide("lab"), sizes = lengths(mods), enrich = e)
}
res <- lapply(c(current = "current", corrected = "corrected", unmerged = "unmerged"), build)

# ---- check: recomputed current vs the figure's own z-scores --------------------
z_org <- as.matrix(read.csv(file.path(ORG, "sample_module_zscore.csv"), row.names = 1, check.names = FALSE))
cat(sprintf("current recomputed vs original z-scores: max |diff| %.3f\n",
            max(abs(res$current$z - z_org[rownames(res$current$z), samples]))))

# ---- draw (same layout as the figure, without the identity block) -------------
pal_z <- rev(brewer.pal(11, "RdBu"))
col_z <- colorRamp2(c(-2.4, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.4), pal_z)
col_aging <- colorRamp2(c(0, 1.3, 3, 6, 10), c("white", "#fee5d9", "#fcae91", "#fb6a4a", "#cb181d"))
CELL <- 5
draw_one <- function(r, title, file) {
  n <- nrow(r$z)
  rl <- sprintf("%s (%d)", rownames(r$z), r$sizes[rownames(r$z)])
  ht_l <- Heatmap(r$z, name = "z-score", col = col_z, cluster_rows = FALSE, cluster_columns = FALSE,
    row_labels = rl, row_names_side = "left", row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7),
    column_names_rot = 0, column_names_centered = TRUE,
    column_split = factor(c("A1", "Y1", rep("R", 10)), levels = c("A1", "Y1", "R")), column_title = NULL,
    column_gap = unit(1, "mm"), width = unit(CELL * 12, "mm"), height = unit(CELL * n, "mm"),
    rect_gp = gpar(col = "white", lwd = 0.6),
    cell_fun = function(j, i, x, y, w, h, fill) grid.text(r$sig_z[i, j], x, y, gp = gpar(fontsize = 5)),
    heatmap_legend_param = list(at = c(-2, -1, 0, 1, 2), direction = "horizontal", legend_width = unit(2, "cm"),
                                labels_gp = gpar(fontsize = 5.5), title_gp = gpar(fontsize = 6.5)))
  ht_r <- Heatmap(r$mat_a, name = "-log10 padj", col = col_aging, cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, column_names_gp = gpar(fontsize = 7), column_names_rot = 45,
    width = unit(CELL * ncol(r$mat_a), "mm"), height = unit(CELL * n, "mm"), rect_gp = gpar(col = "white", lwd = 0.6),
    cell_fun = function(j, i, x, y, w, h, fill) grid.text(r$sig_a[i, j], x, y, gp = gpar(fontsize = 5)),
    heatmap_legend_param = list(at = c(0, 1.3, 3, 6, 10), labels = c("0","1.3","3","6",">=10"), direction = "horizontal",
                                legend_width = unit(2, "cm"), labels_gp = gpar(fontsize = 5.5), title_gp = gpar(fontsize = 6.5)))
  for (dev in c("pdf", "png")) {
    f <- file.path(DIR, paste0(file, ".", dev))
    H <- 3.6 + (n - 11) * CELL / 25.4      # 3.6 in fits 11 rows of 5-mm cells
    if (dev == "pdf") pdf(f, width = 7, height = H) else png(f, width = 7, height = H, units = "in", res = 300)
    draw(ht_l + ht_r, heatmap_legend_side = "bottom", column_title = title, column_title_gp = gpar(fontsize = 9, fontface = "bold"))
    invisible(dev.off())
  }
}
draw_one(res$current,   "Current modules (May: M4 -> M3, M12 -> M2)", "heatmap_current")
draw_one(res$corrected, "Corrected merge (M9 -> M5, M13 -> M6; M4, M12 kept)", "heatmap_corrected")
draw_one(res$unmerged,  "All 13 modules before merging (M1-M13)", "heatmap_unmerged")

# ---- what changed -------------------------------------------------------------
compare <- function(A, B) do.call(rbind, lapply(union(rownames(A$z), rownames(B$z)), function(m) {
  inC <- m %in% rownames(A$z); inK <- m %in% rownames(B$z)
  data.frame(module = m,
    genes_current = if (inC) A$sizes[[m]] else NA, genes_other = if (inK) B$sizes[[m]] else NA,
    z_profile_r = if (inC && inK) round(cor(A$z[m, ], B$z[m, ]), 3) else NA,
    z_stars_changed = if (inC && inK) sum(A$sig_z[m, ] != B$sig_z[m, ]) else NA,
    hallmarks_current = if (inC) paste(names(which(A$sig_a[m, ] != "")), A$sig_a[m, A$sig_a[m, ] != ""], collapse = "; ") else "",
    hallmarks_other = if (inK) paste(names(which(B$sig_a[m, ] != "")), B$sig_a[m, B$sig_a[m, ] != ""], collapse = "; ") else "")
}))
for (o in c("corrected", "unmerged")) {
  cmp <- compare(res$current, res[[o]]); cmp <- cmp[order(as.integer(sub("M", "", cmp$module))), ]
  write.csv(cmp, file.path(DIR, sprintf("module_comparison_current_vs_%s.csv", o)), row.names = FALSE)
  cat("\n== current vs", o, "\n"); print(cmp[, 1:5], row.names = FALSE)
}
