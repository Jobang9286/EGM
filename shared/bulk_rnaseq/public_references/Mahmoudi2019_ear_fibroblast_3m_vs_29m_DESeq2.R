# ------------------------------------------------------------------------------
# Mahmoudi et al. 2018 ear fibroblast, 3 month vs 29 month DEG — from RAW COUNTS
#
# Input  : 18 htseq count files
#            3 month  (n =  8): 3m1 3m2 3m4 3M3c4 3M6c4 3M7c4 3m3c10R 3m5c11L
#           29 month  (n = 10): 29m2 29m6 29m7 29M2c4 29M4c4 29M8c4
#                               29M3C10L 29M4C11R 29M5C10L 29M5C11R
#          DESeq2_all_genes.csv: the original DE table (source of
#          fb_young_vs_old_for_RRHO.csv), used for comparison
# ------------------------------------------------------------------------------
library(DESeq2)

root <- "reproduce"
dir  <- file.path(root, "data", "ear_fib_counts")

files <- list.files(dir, "_htseq_count\\.txt$", full.names = TRUE)
lst <- lapply(files, function(p) {
  x <- read.delim(p, header = FALSE, row.names = 1)
  colnames(x) <- sub("^Fib_|_htseq_count\\.txt$", "", basename(p)); x
})
stopifnot(all(sapply(lst, function(x) identical(rownames(x), rownames(lst[[1]])))))
m <- do.call(cbind, lst)
m <- as.matrix(m[!grepl("^__", rownames(m)), ])            # 23,366 genes x 18

cd <- data.frame(age = factor(ifelse(grepl("^29", colnames(m)), "m29", "m3"),
                              levels = c("m29", "m3")), row.names = colnames(m))
dds <- DESeq(DESeqDataSetFromMatrix(m, cd, ~ age), quiet = TRUE)
res <- results(dds, contrast = c("age", "m3", "m29"))       # + = up in 3 month (young)

tgt <- read.csv(file.path(dir, "DESeq2_all_genes.csv"), row.names = 1)   # 14,310 genes
g   <- intersect(rownames(tgt), rownames(res))
cat(sprintf("shared %d | cor(LFC) %.4f | baseMean ratio median %.4f\n",
            length(g), cor(tgt[g, "log2FoldChange"], res[g, "log2FoldChange"]),
            median(tgt[g, "baseMean"] / res[g, "baseMean"])))

write.csv(as.data.frame(res), file.path(root, "out", "earfib_3m_vs_29m_REPRODUCED.csv"))
