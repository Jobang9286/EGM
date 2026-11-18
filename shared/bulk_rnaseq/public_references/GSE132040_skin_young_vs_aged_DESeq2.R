# ------------------------------------------------------------------------------
# Skin young (1,3 mo) vs aged (21,24,27 mo) DEG — from RAW COUNTS
#
# Input  : matrix_all.csv  54,344 genes x 933 samples, integer counts
#                          (GSE132040 / Tabula Muris Senis bulk, 17 tissues)
#          meta_info.csv   Sample name, source name, tissue, age, range, pair_1, pair_2
# ------------------------------------------------------------------------------
library(data.table); library(DESeq2)

root <- "reproduce"
meta <- fread(file.path(root, "data", "meta_info.csv"),  data.table = FALSE)
cnt  <- fread(file.path(root, "data", "matrix_all.csv"), data.table = FALSE)

# `range` is the age binning:
#   aaa =  1, 3 mo           (n = 12)  -> "Young Skin (1, 3 months)"
#   bbb =  6, 9, 12, 15, 18  (n = 29)
#   ccc = 21, 24, 27 mo      (n = 10)  -> "Aged Skin (21, 24, 27 months)"
# `pair_1` splits aaa into age 1 (abc) vs age 3 (cba)   -> sk_young1_vs_young3.csv
# `pair_2` splits ccc into age 21 (bcd) vs age 27 (dcb) -> sk_old27_vs_old21.csv
sk  <- meta[meta$tissue == "Skin", ]                     # 51 samples
ids <- sk$`Sample name`

m  <- as.matrix(cnt[, ids]); rownames(m) <- cnt$gene
cd <- data.frame(range = factor(sk$range[match(ids, sk$`Sample name`)]), row.names = ids)

# NOTE: all 51 skin samples go into the model (bbb included) and the contrast is
# pulled afterwards.
dds <- DESeq(DESeqDataSetFromMatrix(m, cd, ~ range), quiet = TRUE)
res <- results(dds, contrast = c("range", "aaa", "ccc"))  # + = up in YOUNG

deg <- as.data.frame(res[which(res$pvalue < 0.05), ])     # 15,052 genes
rrho_in <- data.frame(DDE = deg$log2FoldChange, row.names = rownames(deg))  # DDE = log2FC

## ---- checks against the original result tables -----------------------------------
tgt <- read.csv(file.path(root, "data", "skin_young_vs_old.csv"),         row.names = 1)
inp <- read.csv(file.path(root, "data", "skin_young_vs_old_for_RRHO.csv"), row.names = 1)
g <- rownames(tgt)
stopifnot(
  setequal(rownames(deg), g),
  max(abs(deg[g, "log2FoldChange"] - tgt[g, "log2FoldChange"])) < 1e-3,
  max(abs(deg[g, "baseMean"]       - tgt[g, "baseMean"]))       < 1e-6,
  max(abs(rrho_in[g, "DDE"]        - inp[g, "DDE"]))            < 1e-3
)

write.csv(deg,     file.path(root, "out", "skin_young_vs_old_REPRODUCED.csv"))
write.csv(rrho_in, file.path(root, "out", "skin_young_vs_old_for_RRHO_REPRODUCED.csv"))
