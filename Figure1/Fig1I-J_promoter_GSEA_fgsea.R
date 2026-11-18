## ---- H3K27me3 ---------------------------------------------------------------
suppressMessages({library(tidyverse); library(fgsea); library(data.table)})
D <- "chromatin/rank"

gmt <- gmtPathways(file.path(D,"go_bp_mouse.gmt"))
# gmtPathways keys on the first column (GO id); keep a lookup to term names
gmt_lines <- read_lines(file.path(D,"go_bp_mouse.gmt"))
go_name <- tibble(go_id = str_extract(gmt_lines, "^[^\t]+"),
                  term  = str_split_i(gmt_lines, "\t", 2))
cat("gene sets:", length(gmt), "  median size:", median(lengths(gmt)), "\n")

sig <- read_tsv(file.path(D,"h3k27me3_promoter_gene.tsv"), show_col_types=FALSE)
cat("genes:", nrow(sig), "\n")
r <- setNames(sig$log2FC, sig$gene)
r <- r[is.finite(r)]
r <- r - median(r)                      # centre; ranks are unchanged, scale reads cleaner
r <- sort(r, decreasing = TRUE)
cat("ranking: ", sprintf("max %+.2f  median %+.2f  min %+.2f\n", max(r), median(r), min(r)))
cat("top 6 (4FP-high):", paste(names(head(r,6)), collapse=", "), "\n")
cat("bottom 6 (CP-high):", paste(names(tail(r,6)), collapse=", "), "\n\n")

set.seed(1)
fg <- fgsea(pathways = gmt, stats = r, minSize = 15, maxSize = 400, nPermSimple = 10000)
fg <- as_tibble(fg) |> left_join(go_name, by = c("pathway"="go_id"))
cat("pathways tested:", nrow(fg), " | padj < 0.05:", sum(fg$padj < 0.05, na.rm=TRUE), "\n")
cat("  NES > 0 (4FP-specific):", sum(fg$padj<0.05 & fg$NES>0, na.rm=TRUE),
    " NES < 0 (CP-specific):",  sum(fg$padj<0.05 & fg$NES<0, na.rm=TRUE), "\n\n")
fg |> filter(padj<0.05) |> arrange(desc(NES)) |> select(term, NES, padj, size) |> head(10) |> print()
fg |> filter(padj<0.05) |> arrange(NES)       |> select(term, NES, padj, size) |> head(10) |> print()
saveRDS(fg, file.path(D,"fgsea_h3k27me3.rds"))

## ---- ATAC (and H3K27ac) -----------------------------------------------------
suppressMessages({library(tidyverse); library(fgsea)})
D <- "chromatin/rank"
gmt <- gmtPathways(file.path(D,"go_bp_mouse.gmt"))
gl  <- read_lines(file.path(D,"go_bp_mouse.gmt"))
nm  <- tibble(pathway = str_extract(gl,"^[^\t]+"), term = str_split_i(gl,"\t",2))

run <- function(f, tag){
  s <- read_tsv(file.path(D,f), show_col_types=FALSE, progress=FALSE)
  r <- setNames(s$log2FC, s$gene); r <- r[is.finite(r)]; r <- sort(r - median(r), decreasing=TRUE)
  set.seed(1)
  fg <- fgsea(gmt, r, minSize=15, maxSize=400, nPermSimple=10000) |>
    as_tibble() |> left_join(nm, by="pathway")
  saveRDS(fg, file.path(D, paste0("fgsea_", tag, ".rds")))
  cat(sprintf("\n%-9s tested %d | padj<0.05: %d  (NES>0 %d, NES<0 %d)\n",
              tag, nrow(fg), sum(fg$padj<0.05, na.rm=TRUE),
              sum(fg$padj<0.05 & fg$NES>0, na.rm=TRUE),
              sum(fg$padj<0.05 & fg$NES<0, na.rm=TRUE)))
  fg
}
at <- run("atac_promoter_gene.tsv",    "atac")
ac <- run("h3k27ac_promoter_gene.tsv", "h3k27ac")

cat("\n=== ATAC, top P-BP-specific (NES > 0) ===\n")
at |> filter(pval<0.05, NES>0) |> arrange(desc(NES)) |> head(15) |>
  mutate(NES=round(NES,2), pval=signif(pval,2), padj=signif(padj,2)) |>
  select(term, NES, pval, padj, size) |> as.data.frame() |> print(right=FALSE)
cat("\n=== H3K27ac, top P-BP-specific (NES > 0) ===\n")
ac |> filter(pval<0.05, NES>0) |> arrange(desc(NES)) |> head(15) |>
  mutate(NES=round(NES,2), pval=signif(pval,2), padj=signif(padj,2)) |>
  select(term, NES, pval, padj, size) |> as.data.frame() |> print(right=FALSE)
