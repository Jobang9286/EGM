## ---- paths -------------------------------------------------------------------------------------------
root <- "monkey"
outdir <- file.path(root,"functional_study"); dir.create(outdir, showWarnings=FALSE)

## DESeq2 results per organ
## (res_4organ.rds: list Brain / Liver / Kidney / Spleen)
res_all <- readRDS(file.path(outdir,"res_4organ.rds")); resB <- res_all$Brain; resL <- res_all$Liver

## ---- ortholog Entrez map, go_up(), Brain GO-BP -------------------------------------------------------
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(data.table); library(ggplot2)})
om2 <- fread(file.path(root,"ortholog_map.csv"))
om2 <- om2[!is.na(entrez) & entrez!=""][!duplicated(gene_id), .(gene_id, entrez=as.character(entrez))]
go_up <- function(res){
  d <- merge(as.data.table(res), om2, by="gene_id")
  uni <- unique(d[!is.na(padj)]$entrez)
  up  <- unique(d[!is.na(padj) & padj<0.05 & log2FoldChange>0]$entrez)
  enrichGO(up, universe=uni, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="BP",
           pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.1,
           minGSSize=10, maxGSSize=500, readable=TRUE)
}
egB <- go_up(resB)
saveRDS(egB, file.path(outdir,"GOup_Brain.rds"))
cat("Brain EGM-UP enriched BP terms:", nrow(as.data.frame(egB)), "\n\n")
print(head(as.data.frame(egB)[,c("Description","Count","GeneRatio","p.adjust")], 35))

## ---- Liver GO-BP --------------------------------------------------------------------------------------
egL <- go_up(resL)
saveRDS(egL, file.path(outdir,"GOup_Liver.rds"))
dfL <- as.data.frame(egL)
cat("Liver EGM-UP enriched BP terms:", nrow(dfL), "\n\n")

## ---- Kidney GO-BP (setup repeated) -------------------------------------------------------------------
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(data.table); library(ggplot2)})
root <- "monkey"; outdir <- file.path(root,"functional_study")
res_all <- readRDS(file.path(outdir,"res_4organ.rds"))   # Brain/Liver/Kidney/Spleen DESeq2
resK <- res_all$Kidney
om2 <- fread(file.path(root,"ortholog_map.csv"))
om2 <- om2[!is.na(entrez) & entrez!=""][!duplicated(gene_id), .(gene_id, entrez=as.character(entrez))]
go_up <- function(res){
  d <- merge(as.data.table(res), om2, by="gene_id")
  uni <- unique(d[!is.na(padj)]$entrez)
  up  <- unique(d[!is.na(padj) & padj<0.05 & log2FoldChange>0]$entrez)
  enrichGO(up, universe=uni, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="BP",
           pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.1,
           minGSSize=10, maxGSSize=500, readable=TRUE)
}
egK <- go_up(resK); saveRDS(egK, file.path(outdir,"GOup_Kidney.rds"))
dfK <- as.data.frame(egK)
cat("Kidney EGM-UP enriched BP terms:", nrow(dfK), "\n\n")

## ---- Kidney KEGG -----------------------------------------------------------------------------------
suppressMessages(library(clusterProfiler))
d <- merge(as.data.table(resK), om2, by="gene_id")
up  <- unique(d[!is.na(padj) & padj<0.05 & log2FoldChange>0]$entrez)
uni <- unique(d[!is.na(padj)]$entrez)
cat("KEGG: up genes", length(up), " universe", length(uni), "\n")
kk <- tryCatch(enrichKEGG(gene=up, organism="hsa", universe=uni,
                          pvalueCutoff=0.1, qvalueCutoff=0.2, minGSSize=8, maxGSSize=500),
               error=function(e) {cat("KEGG error:", conditionMessage(e), "\n"); NULL})
if(!is.null(kk) && nrow(as.data.frame(kk))>0){
  dk <- as.data.frame(kk)
  saveRDS(kk, file.path(outdir,"KEGGup_Kidney.rds"))
  cat("\nKEGG enriched:", nrow(dk), "pathways\n\n")
  print(head(dk[order(dk$p.adjust), c("Description","Count","GeneRatio","p.adjust")], 30), row.names=FALSE)
} else cat("no KEGG result (likely offline)\n")
cat("\nReactomePA installed:", requireNamespace("ReactomePA", quietly=TRUE), "\n")

## ---- Spleen GO-BP + KEGG, EGM-up (Supp Fig 12L) ----------------------------------------------------
resS <- res_all$Spleen
egS <- go_up(resS); saveRDS(egS, file.path(outdir,"GOup_Spleen.rds"))
dfS <- as.data.frame(egS); dfS$ratio <- sapply(strsplit(dfS$GeneRatio,"/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
cat("Spleen EGM-UP GO-BP terms:", nrow(dfS), "\n=== TOP 25 ===\n")
print(head(dfS[order(dfS$p.adjust), c("Description","Count","p.adjust")], 25), row.names=FALSE)

# KEGG
dS <- merge(as.data.table(resS), om2, by="gene_id")
upS  <- unique(dS[!is.na(padj) & padj<0.05 & log2FoldChange>0]$entrez)
uniS <- unique(dS[!is.na(padj)]$entrez)
kkS <- enrichKEGG(gene=upS, organism="hsa", universe=uniS, pvalueCutoff=1, qvalueCutoff=1, minGSSize=5, maxGSSize=600)
dallS <- as.data.frame(kkS); saveRDS(kkS, file.path(outdir,"KEGGup_Spleen.rds"))
cat("\n=== KEGG top 20 (Spleen EGM-UP) ===\n")
print(head(dallS[order(dallS$p.adjust), c("Description","Count","GeneRatio","p.adjust")], 20), row.names=FALSE)

## ---- Spleen GO-BP + KEGG, EGM-down (Spleen_down figure) ---------------------------------------------
d <- merge(as.data.table(resS), om2, by="gene_id")
dn  <- unique(d[!is.na(padj) & padj<0.05 & log2FoldChange<0]$entrez)
uni <- unique(d[!is.na(padj)]$entrez)
cat("Spleen DOWN genes:", length(dn), " universe:", length(uni), "\n")
egSd <- enrichGO(dn, universe=uni, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="BP",
                 pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.1, minGSSize=10, maxGSSize=500, readable=TRUE)
saveRDS(egSd, file.path(outdir,"GOdown_Spleen.rds"))
dfSd <- as.data.frame(egSd); dfSd$ratio <- sapply(strsplit(dfSd$GeneRatio,"/"),function(x) as.numeric(x[1])/as.numeric(x[2]))
cat("\n=== GO-BP DOWN top 22 ===\n")
print(head(dfSd[order(dfSd$p.adjust), c("Description","Count","p.adjust")], 22), row.names=FALSE)
kkSd <- enrichKEGG(dn, organism="hsa", universe=uni, pvalueCutoff=0.1, qvalueCutoff=0.2, minGSSize=8, maxGSSize=500)
saveRDS(kkSd, file.path(outdir,"KEGGdown_Spleen.rds"))
dallSd <- as.data.frame(kkSd)
cat("\n=== KEGG DOWN (padj<0.1):", nrow(dallSd), "===\n")
print(head(dallSd[order(dallSd$p.adjust), c("Description","Count","GeneRatio","p.adjust")], 20), row.names=FALSE)
