#!/usr/bin/env Rscript
# Module-wise enrichment (GO BP/CC/MF + KEGG + Reactome) for EGM hdWGCNA modules
#   nohup Rscript run_module_enrichment.R > enrichment_run.log 2>&1 &
# 결과 CSV 를 로컬 tables/ 로 내려받아 plot_module_enrichment.R 로 그림을 그림.
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ReactomePA)
})

base   <- "EGM_multiflex/1_module_heatmap_trajectory_script"
indir  <- file.path(base, "module_gene_lists")
outdir <- file.path(base, "GO_results", "tables")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

modules <- c("M1","M2","M3","M5","M6","M7","M8","M9","M10","M11","M13")
P_CUT <- 0.05
Q_CUT <- 0.05

log_line <- function(...) { cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = ""); flush(stdout()) }

n_terms <- function(res) if (is.null(res)) NA_integer_ else nrow(as.data.frame(res))

write_res <- function(res, path, label) {
  n <- n_terms(res)
  if (is.na(n) || n == 0) { log_line(label, ": 0 terms (파일 생성 안 함)"); return(invisible(n)) }
  write.csv(as.data.frame(res), path, row.names = FALSE)
  log_line(label, ": ", n, " terms -> ", basename(path))
  invisible(n)
}

summary_rows <- list()
add_row <- function(m, db, n) summary_rows[[length(summary_rows) + 1]] <<- data.frame(module = m, db = db, n_terms = n)

for (m in modules) {
  f <- file.path(indir, paste0(m, "_genes.csv"))
  if (!file.exists(f)) { log_line(m, ": gene list 없음 -> skip"); next }
  genes <- read.csv(f, stringsAsFactors = FALSE)$gene_name
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  log_line("=== ", m, " : ", length(genes), " genes ===")

  for (ont in c("BP", "CC", "MF")) {
    ego <- tryCatch(
      enrichGO(gene = genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = ont,
               pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT),
      error = function(e) { log_line(m, " GO_", ont, " ERROR: ", conditionMessage(e)); NULL })
    n <- write_res(ego, file.path(outdir, paste0(m, "_GO_", ont, ".csv")), paste0(m, " GO_", ont))
    add_row(m, paste0("GO_", ont), n)
  }

  map <- tryCatch(suppressWarnings(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)),
                  error = function(e) NULL)
  eg <- unique(map$ENTREZID)
  log_line(m, ": ENTREZ 매핑 ", length(eg), "/", length(genes))

  ekegg <- tryCatch(
    enrichKEGG(gene = eg, organism = "mmu", keyType = "kegg",
               pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT),
    error = function(e) { log_line(m, " KEGG ERROR: ", conditionMessage(e)); NULL })
  if (!is.null(ekegg) && n_terms(ekegg) > 0)
    ekegg <- setReadable(ekegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  n <- write_res(ekegg, file.path(outdir, paste0(m, "_KEGG.csv")), paste0(m, " KEGG"))
  add_row(m, "KEGG", n)

  erea <- tryCatch(
    enrichPathway(gene = eg, organism = "mouse", pAdjustMethod = "BH",
                  pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT, readable = TRUE),
    error = function(e) { log_line(m, " Reactome ERROR: ", conditionMessage(e)); NULL })
  n <- write_res(erea, file.path(outdir, paste0(m, "_Reactome.csv")), paste0(m, " Reactome"))
  add_row(m, "Reactome", n)
}

summ <- do.call(rbind, summary_rows)
write.csv(summ, file.path(outdir, "enrichment_summary_counts.csv"), row.names = FALSE)
log_line("=== 완료 ===")
print(summ)
