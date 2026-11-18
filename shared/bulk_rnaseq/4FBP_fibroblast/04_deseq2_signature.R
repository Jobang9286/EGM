# 4FP vs CP differential expression -> ranked signature for RedRibbon
suppressMessages({library(DESeq2); library(apeglm)})

cnt  <- read.csv("4FP_CP_counts_symbol.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("4FP_CP_meta.csv", row.names = 1)
meta <- meta[colnames(cnt), , drop = FALSE]
stopifnot(identical(rownames(meta), colnames(cnt)))

cnt <- round(as.matrix(cnt))                      # RSEM expected_counts -> integers
keep <- rowSums(cnt >= 10) >= 3
cat("genes:", nrow(cnt), "-> kept:", sum(keep), "\n")
cnt <- cnt[keep, ]

meta$group <- factor(meta$group, levels = c("CP", "4FP"))   # CP = reference
dds <- DESeqDataSetFromMatrix(cnt, meta, ~ group)
dds <- DESeq(dds, quiet = TRUE)

res <- results(dds, contrast = c("group", "4FP", "CP"))
shr <- lfcShrink(dds, coef = "group_4FP_vs_CP", type = "apeglm", quiet = TRUE)

sig <- data.frame(
  id        = rownames(res),
  baseMean  = res$baseMean,
  lfc_raw   = res$log2FoldChange,
  lfc       = shr$log2FoldChange,          # apeglm-shrunk (3v3 -> essential)
  pvalue    = res$pvalue,
  padj      = res$padj,
  stat      = res$stat,
  row.names = NULL
)
sig <- sig[!is.na(sig$pvalue) & !is.na(sig$lfc), ]
# RedRibbon / RRHO ranking metric: signed significance
sig$rank_metric <- -log10(sig$pvalue) * sign(sig$lfc)
sig <- sig[order(-sig$rank_metric), ]

write.csv(sig, "analysis/4FP_vs_CP_signature.csv", row.names = FALSE)

cat("\n== summary ==\n")
cat("tested genes      :", nrow(sig), "\n")
cat("padj < 0.05       :", sum(sig$padj < 0.05, na.rm = TRUE), "\n")
cat("padj<0.05 & |lfc|>0.5:", sum(sig$padj < 0.05 & abs(sig$lfc) > 0.5, na.rm = TRUE), "\n")
cat("  up  :", sum(sig$padj < 0.05 & sig$lfc >  0.5, na.rm = TRUE), "\n")
cat("  down:", sum(sig$padj < 0.05 & sig$lfc < -0.5, na.rm = TRUE), "\n")
cat("\n-- top 15 UP in 4FP --\n");   print(head(sig[, c("id","lfc","padj")], 15), digits = 3)
cat("\n-- top 15 DOWN in 4FP --\n"); print(head(sig[order(sig$rank_metric), c("id","lfc","padj")], 15), digits = 3)

saveRDS(list(dds = dds, sig = sig), "analysis/deseq2_objects.rds")
