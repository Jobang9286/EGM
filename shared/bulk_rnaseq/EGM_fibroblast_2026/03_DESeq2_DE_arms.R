## ---- DESeq2 model --------------------------------------------------------------
suppressPackageStartupMessages(library(DESeq2))

C  <- read.csv("egm_blood/counts_symbol.csv", row.names = 1, check.names = FALSE)
mp <- read.delim("egm_blood/mapfile.tsv")
mp$group <- factor(mp$TREATMENT, levels = c("untreated","egm","eno1","gclc","map1s"))
rownames(mp) <- mp$sample
C <- as.matrix(C[, mp$sample])
cat("counts:", nrow(C), "genes x", ncol(C), "samples\n")
print(table(mp$group))

keep <- rowSums(C >= 10) >= 3
C <- C[keep, ]
cat("after filter (count>=10 in >=3 samples):", nrow(C), "genes\n")

# one model over all 9 samples: the n=1 arms borrow dispersion from the
# replicated GFP and EGM groups, which is the only way to get usable stats there
dds <- DESeqDataSetFromMatrix(round(C), mp, ~ group)
dds <- DESeq(dds, quiet = TRUE)
cat("\ncoefficients:", paste(resultsNames(dds), collapse = ", "), "\n")

## ---- each arm vs untreated -----------------------------------------------------
OUT <- "RRHO_new/reproduced"; dir.create(OUT, showWarnings = FALSE)

arms <- c(EGM = "egm", Eno1 = "eno1", Gclc = "gclc", Map1s = "map1s")
DE <- list()
for (nm in names(arms)) {
  r <- results(dds, contrast = c("group", arms[[nm]], "untreated"))
  d <- data.frame(gene = rownames(r), DDE = r$log2FoldChange,
                  pvalue = r$pvalue, padj = r$padj)
  d <- d[!is.na(d$DDE) & !is.na(d$pvalue), ]
  DE[[nm]] <- d
  write.csv(d[, c("gene","DDE")], file.path(OUT, sprintf("%s_vs_nEGM_for_RRHO.csv", nm)),
            row.names = FALSE)
  cat(sprintf("%-6s vs nEGM : %5d genes tested | padj<0.05 %5d | |lfc|>1 %5d | DDE %.2f .. %.2f\n",
              nm, nrow(d), sum(d$padj < .05, na.rm = TRUE),
              sum(abs(d$DDE) > 1), min(d$DDE), max(d$DDE)))
}
saveRDS(DE, file.path(OUT, "DE_arms.rds"))
