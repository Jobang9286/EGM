# ------------------------------------------------------------------
# Fig 3K check: are the lipid and telomere pathways the Q2-Q3 and Q4-Q5 parts of pseudotime?
#
# Same method as Fig3K_window_GSEA_lollipops.R:
#   1,200 pseudotime-dynamic genes (trajectory_env.rds$z_mat), ranked by
#   window score = mean z inside the window - mean z outside, then clusterProfiler::gseGO (BP).
# Windows are pseudotime quintiles of cells. The heatmap matrix has 80 equal-cell bins with 2 bins
# trimmed at each end (76 columns = original bins 3-78), so quintile q = original bins 16(q-1)+1 .. 16q:
#   Q1 = columns 1-14 (bins 3-16), Q2 = 15-30, Q3 = 31-46, Q4 = 47-62, Q5 = 63-76 (bins 65-78).
# Q2-Q3 is exactly the "transition" window of Fig 3K (columns 15-46); the Fig 3K "late" window was
# columns 46-76, i.e. Q4-Q5 plus the last bin of Q3.
# Outputs: quintile_gsea_<window>.csv, quintile_terms_summary.csv
# ------------------------------------------------------------------
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Mm.eg.db) })

z   <- readRDS("EGM_oracle/data/trajectory_env.rds")$z_mat
OUT <- "4FP_plasma/EGM_Fig3K_pathways"
stopifnot(ncol(z) == 76)
q_of_col <- ceiling((seq_len(ncol(z)) + 2) / 16)          # column j = original bin j + 2
print(table(quintile = q_of_col))

win_rank <- function(cols) sort(rowMeans(z[, cols]) - rowMeans(z[, -cols]), decreasing = TRUE)
run_gsea <- function(cols) {
  set.seed(42)
  gseGO(geneList = win_rank(cols), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
        minGSSize = 10, maxGSSize = 2000, pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE)@result
}

windows <- list("Q2-Q3" = which(q_of_col %in% 2:3), "Q4-Q5" = which(q_of_col %in% 4:5),
                Q1 = which(q_of_col == 1), Q2 = which(q_of_col == 2), Q3 = which(q_of_col == 3),
                Q4 = which(q_of_col == 4), Q5 = which(q_of_col == 5))
TERMS <- c("lipid metabolic process", "fatty acid metabolic process",
           "telomere maintenance", "telomere organization")

res <- lapply(names(windows), function(w) {
  g <- run_gsea(windows[[w]])
  write.csv(g[, c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust")],
            file.path(OUT, sprintf("quintile_gsea_%s.csv", w)), row.names = FALSE)
  sig <- g[g$NES > 0 & g$p.adjust < 0.05, ]; sig <- sig[order(-sig$NES), ]
  rows <- lapply(TERMS, function(t) {
    h <- g[g$Description == t, ]
    data.frame(window = w, term = t, NES = round(h$NES, 2), FDR = signif(h$p.adjust, 2),
               rank_among_sig_pos = if (t %in% sig$Description) match(t, sig$Description) else NA,
               n_sig_pos = nrow(sig))
  })
  cat(sprintf("\n[%s] columns %d-%d | significant positive terms: %d\n  top 8: %s\n", w,
              min(windows[[w]]), max(windows[[w]]), nrow(sig), paste(head(sig$Description, 8), collapse = " | ")))
  list(sum = do.call(rbind, rows), sig = sig$ID)
})
names(res) <- names(windows)
summ <- do.call(rbind, lapply(res, `[[`, "sum"))
write.csv(summ, file.path(OUT, "quintile_terms_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)

# overlap with the Fig 3K pathway sets (249 transition / 109 late, NES > 0 & FDR < 0.05)
fig_t <- read.csv(file.path(OUT, "Fig3K_transition_window_249pathways.csv"))$GO_ID
fig_l <- read.csv(file.path(OUT, "Fig3K_late_window_109pathways.csv"))$GO_ID
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
cat(sprintf("\nQ2-Q3 vs Fig 3K transition (249): shared %d, Jaccard %.2f\n",
            length(intersect(res[["Q2-Q3"]]$sig, fig_t)), jac(res[["Q2-Q3"]]$sig, fig_t)))
cat(sprintf("Q4-Q5 vs Fig 3K late (109):       shared %d, Jaccard %.2f\n",
            length(intersect(res[["Q4-Q5"]]$sig, fig_l)), jac(res[["Q4-Q5"]]$sig, fig_l)))
