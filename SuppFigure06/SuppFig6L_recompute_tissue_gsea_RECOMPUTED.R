# ------------------------------------------------------------------------------
# NES bars (Liver / Kidney / Brain) — GSEA on the HGPS tissue bulk RNA-seq
#
# Counts: EGM_project_for_deconvolution/{liver,kidney,brain} (featureCounts),
#         HGPS_control_1-3 vs HGPS_EGM_1-3 per tissue.
# Rankings (all gseGO BP, SYMBOL, minGSSize 10, maxGSSize 5000, seeded):
#   lfc_hgps  : DESeq2 on the 6 HGPS samples, ranked by log2FoldChange (EGM vs control)
#   stat_hgps : same DESeq2, ranked by the Wald statistic
#   lfc_all   : DESeq2 on every sample of the tissue (~group), contrast HGPS_EGM vs HGPS_control, log2FC
# Output: tissue_gsea_<tissue>_<ranking>.csv, nes_match_table.csv
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({ library(DESeq2); library(clusterProfiler); library(org.Mm.eg.db) })
OUT  <- "4FP_plasma/EGM_SuppFig7_NES"
ROOT <- "EGM_project_for_deconvolution"
files <- c(Liver = "liver/counts_liver.csv", Kidney = "kidney/kidney_bulk_data/counts_kidney.txt", Brain = "brain/counts_brain.txt")

# reference NES of the figure bars (NES_ppt), label as printed -> GO term
TARGET <- tibble::tribble(
  ~tissue,  ~label,                              ~term,                                   ~NES_ppt,
  "Liver",  "B cell signaling pathway",          "B cell receptor signaling pathway",     1.661,
  "Liver",  "Cell division",                     "cell division",                         1.379,
  "Liver",  "DNA repair",                        "DNA repair",                            1.270,
  "Liver",  "Cell cycle",                        "cell cycle",                            1.216,
  "Liver",  "Cell development",                  "cell development",                      1.204,
  "Kidney", "Calcium ion transport",             "calcium ion transport",                 1.194,
  "Kidney", "Cellular homeostasis",              "cellular homeostasis",                  1.172,
  "Kidney", "Ion transmembrane transport",       "monoatomic ion transmembrane transport",1.159,
  "Kidney", "Cell-cell adhesion",                "cell-cell adhesion",                    1.127,
  "Brain",  "Response to lectin",                "response to lectin",                    1.415,
  "Brain",  "Neural crest cell migration",       "neural crest cell migration",           1.320,
  "Brain",  "Neural crest cell differentiation", "neural crest cell differentiation",     1.314,
  "Brain",  "TLR signaling pathway",             "toll-like receptor signaling pathway", -1.372)

read_counts <- function(f) {
  x <- read.delim(file.path(ROOT, f), sep = if (grepl("\\.csv$", f)) "," else "\t", check.names = FALSE)
  x <- x[!duplicated(x$Geneid), ]; m <- as.matrix(x[, -1]); rownames(m) <- x$Geneid
  colnames(m) <- sub("\\.bam$", "", colnames(m)); m
}
group_of <- function(s) sub("_[0-9]+$", "", sub("^[A-Za-z0-9]+_", "", s))   # Liver_HGPS_EGM_1 -> HGPS_EGM

gsea <- function(r) {
  set.seed(42)
  gseGO(sort(r, decreasing = TRUE), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
        minGSSize = 10, maxGSSize = 5000, pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE)@result
}

res <- list()
for (tis in names(files)) {
  m <- read_counts(files[[tis]]); g <- group_of(colnames(m))
  cat(tis, ":", paste(names(table(g)), table(g), collapse = ", "), "\n")
  h <- g %in% c("HGPS_control", "HGPS_EGM")
  dds <- DESeqDataSetFromMatrix(m[, h], data.frame(group = factor(g[h], c("HGPS_control", "HGPS_EGM"))), ~group)
  dds <- DESeq(dds[rowSums(counts(dds)) > 0, ], quiet = TRUE)
  r6  <- as.data.frame(results(dds, contrast = c("group", "HGPS_EGM", "HGPS_control")))
  r6  <- r6[!is.na(r6$stat), ]
  dda <- DESeqDataSetFromMatrix(m, data.frame(group = factor(g)), ~group)
  dda <- DESeq(dda[rowSums(counts(dda)) > 0, ], quiet = TRUE)
  ra  <- as.data.frame(results(dda, contrast = c("group", "HGPS_EGM", "HGPS_control")))
  ra  <- ra[!is.na(ra$log2FoldChange), ]
  write.csv(r6, file.path(OUT, sprintf("DESeq2_%s_HGPS_EGM_vs_control.csv", tis)))
  ranks <- list(lfc_hgps  = setNames(r6$log2FoldChange, rownames(r6)),
                stat_hgps = setNames(r6$stat,           rownames(r6)),
                lfc_all   = setNames(ra$log2FoldChange, rownames(ra)))
  for (k in names(ranks)) {
    gs <- gsea(ranks[[k]])
    write.csv(gs[, c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust")],
              file.path(OUT, sprintf("tissue_gsea_%s_%s.csv", tis, k)), row.names = FALSE)
    res[[paste(tis, k)]] <- data.frame(tissue = tis, ranking = k, gs[, c("Description", "setSize", "NES", "p.adjust")])
  }
}
all <- do.call(rbind, res)
tab <- merge(TARGET, all, by.x = c("tissue", "term"), by.y = c("tissue", "Description"), all.x = TRUE)
tab <- tab[order(tab$tissue, tab$ranking, -tab$NES_ppt), ]
tab$NES <- round(tab$NES, 3); tab$p.adjust <- signif(tab$p.adjust, 2)
write.csv(tab, file.path(OUT, "nes_match_table.csv"), row.names = FALSE)
print(tab[, c("tissue", "ranking", "label", "setSize", "NES_ppt", "NES", "p.adjust")], row.names = FALSE)
cat("\nmean |NES - NES_ppt| by ranking:\n")
print(aggregate(abs(NES - NES_ppt) ~ ranking, tab, mean))
