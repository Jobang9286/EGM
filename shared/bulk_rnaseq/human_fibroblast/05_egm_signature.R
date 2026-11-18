# Human EGM vs Control (3 vs 3, primary human fibroblasts) -> EGM signature.
# Same recipe as ../4FBP_fibroblast/04_deseq2_signature.R so the two signatures
# are built identically: RSEM expected counts, count>=10 in >=3 samples, DESeq2,
# apeglm-shrunk LFC (3v3 -> essential), padj < 0.05 & |lfc| > 0.5 for the modules.
suppressMessages({library(DESeq2); library(apeglm)})

find_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  cands <- c(if (length(f)) dirname(dirname(normalizePath(f[1]))), getwd(), dirname(getwd()))
  hit <- cands[dir.exists(file.path(cands, "data"))]
  if (!length(hit)) stop("cannot locate the data/ directory; run from the project root")
  hit[1]
}
ROOT <- find_root(); D <- file.path(ROOT, "data")
SRC  <- file.path(dirname(ROOT), "human_egm_results")   # share bundle written by 04_make_share_bundle.sh

cnt  <- read.delim(file.path(SRC, "human_EGM_counts.expected_count.tsv"),
                   row.names = 1, check.names = FALSE)
ss   <- read.delim(file.path(SRC, "human_EGM_samplesheet.tsv"))
meta <- data.frame(row.names = ss$sample,
                   group = factor(ifelse(ss$TREATMENT == "untreated", "Control", "EGM"),
                                  levels = c("Control", "EGM")))        # Control = reference
cnt  <- round(as.matrix(cnt[, rownames(meta)]))
keep <- rowSums(cnt >= 10) >= 3
cat("genes:", nrow(cnt), "-> kept:", sum(keep), "\n")
cnt  <- cnt[keep, ]

dds <- DESeqDataSetFromMatrix(cnt, meta, ~ group)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("group", "EGM", "Control"))
shr <- lfcShrink(dds, coef = "group_EGM_vs_Control", type = "apeglm", quiet = TRUE)

map <- read.delim(file.path(SRC, "human_gene_id_to_symbol.tsv"))
sym <- setNames(map$Gene.Name, map$Gene.ID)
sig <- data.frame(
  gene_id  = rownames(res),
  id       = unname(sym[rownames(res)]),
  baseMean = res$baseMean,
  lfc_raw  = res$log2FoldChange,
  lfc      = shr$log2FoldChange,          # apeglm-shrunk
  pvalue   = res$pvalue,
  padj     = res$padj,
  stat     = res$stat,
  row.names = NULL)
miss <- is.na(sig$id) | sig$id == ""
sig$id[miss] <- sig$gene_id[miss]
sig$named <- !grepl("^ENSG\\d+$", sig$id)   # FALSE = no gene symbol in Ensembl 110
sig <- sig[!is.na(sig$pvalue) & !is.na(sig$lfc), ]
# one row per symbol (71 symbols are carried by >1 Ensembl id): keep the most expressed
sig <- sig[order(sig$id, -sig$baseMean), ]
sig <- sig[!duplicated(sig$id), ]
# RRHO / RedRibbon ranking metric; 74 genes have pvalue == 0, so cap before the log
sig$rank_metric <- -log10(pmax(sig$pvalue, .Machine$double.xmin)) * sign(sig$lfc)
sig <- sig[order(-sig$rank_metric), ]

write.csv(sig, file.path(D, "EGM_vs_Control_signature.csv"), row.names = FALSE)
saveRDS(list(dds = dds, sig = sig), file.path(D, "deseq2_objects.rds"))

cat("\n== summary (named genes only) ==\n")
s <- sig[sig$named, ]
cat("tested genes            :", nrow(s), "\n")
cat("padj < 0.05             :", sum(s$padj < 0.05, na.rm = TRUE), "\n")
cat("padj<0.05 & |lfc|>0.5   :", sum(s$padj < 0.05 & abs(s$lfc) > 0.5, na.rm = TRUE),
    " (up", sum(s$padj < 0.05 & s$lfc >  0.5, na.rm = TRUE),
    "/ down", sum(s$padj < 0.05 & s$lfc < -0.5, na.rm = TRUE), ")\n")
cat("padj<0.05 & |lfc|>1     :", sum(s$padj < 0.05 & abs(s$lfc) > 1, na.rm = TRUE),
    " (up", sum(s$padj < 0.05 & s$lfc >  1, na.rm = TRUE),
    "/ down", sum(s$padj < 0.05 & s$lfc < -1, na.rm = TRUE), ")\n")
cat("shrunk vs raw lfc  r    :", round(cor(s$lfc, s$lfc_raw), 3), "\n")
cat("\n-- top 15 UP in EGM --\n");   print(head(s[, c("id","lfc","padj")], 15), digits = 3)
cat("\n-- top 15 DOWN in EGM --\n"); print(head(s[order(s$rank_metric), c("id","lfc","padj")], 15), digits = 3)
