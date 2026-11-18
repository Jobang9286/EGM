#!/usr/bin/env bash
export ROOT="$PWD"

## ---- export full GO / KEGG tables -----------------------------------------
cd "$ROOT/monkey/functional_study" && mkdir -p GO_tables && Rscript -e '
suppressWarnings({
setwd(file.path(Sys.getenv("ROOT"), "monkey/functional_study"))
OUT <- "GO_tables"
grab <- function(f){ r <- readRDS(f)@result; r$ratio <- sapply(strsplit(r$GeneRatio,"/"), function(x) as.numeric(x[1])/as.numeric(x[2])); r }

## full tables
full <- list(GOup_Brain="GOup_Brain.rds", GOup_Liver="GOup_Liver.rds", GOup_Kidney="GOup_Kidney.rds",
             GOup_Spleen="GOup_Spleen.rds", GOdown_Spleen="GOdown_Spleen.rds",
             KEGGup_Kidney="KEGGup_Kidney.rds", KEGGup_Spleen="KEGGup_Spleen.rds", KEGGdown_Spleen="KEGGdown_Spleen.rds")
for(n in names(full)) write.csv(grab(full[[n]]), file.path(OUT, paste0("FULL_", n, ".csv")), row.names=FALSE)
})' 2>&1 | grep -v "^Warning\|^In addition\|package .* was built"

## ---- Brain / Liver gene sets (universe, EGM-up Entrez) ---------------------
cd "$ROOT/scratch" && Rscript -e '
suppressMessages({library(clusterProfiler); library(data.table)})
cat("clusterProfiler:", as.character(packageVersion("clusterProfiler")), "\n")
res_all <- readRDS(file.path(Sys.getenv("ROOT"), "monkey/functional_study/res_4organ.rds"))
om <- fread(file.path(Sys.getenv("ROOT"), "monkey/ortholog_map.csv"))
om2 <- om[!is.na(entrez) & entrez!=""][!duplicated(gene_id), .(gene_id, entrez=as.character(entrez))]
for(org in c("Brain","Liver")){
  d <- merge(as.data.table(res_all[[org]]), om2, by="gene_id")
  uni <- unique(d[!is.na(padj)]$entrez); g <- unique(d[!is.na(padj)&padj<0.05&log2FoldChange>0]$entrez)
  cat(sprintf("%s: universe=%d  UP genes(entrez)=%d\n", org, length(uni), length(g)))
  saveRDS(list(uni=uni,g=g), paste0(org,"_sets.rds"))
}' 2>&1 | grep -v "^Warning\|built under\|^In addition"

## ---- Brain KEGG ORA (EGM-up) -----------------------------------------------
cd "$ROOT/scratch" && Rscript -e '
suppressMessages(library(clusterProfiler))
s <- readRDS("Brain_sets.rds")
kk <- enrichKEGG(s$g, organism="hsa", universe=s$uni, pvalueCutoff=1, qvalueCutoff=1, minGSSize=5, maxGSSize=600)
if(is.null(kk)){ cat("NULL returned (KEGG fetch failed)\n") } else {
saveRDS(kk, "KEGGup_Brain.rds"); r <- kk@result
cat("총 pathway:", nrow(r), "  p.adj<0.05:", sum(r$p.adjust<0.05,na.rm=TRUE), "\n\n")
r <- r[order(r$p.adjust),]; s2 <- r[r$p.adjust<0.05,]
cat("--- Brain ORA-KEGG (UP), p.adj<0.05 전체 ---\n")
print(data.frame(ID=s2$ID, Description=substr(s2$Description,1,44), GeneRatio=s2$GeneRatio,
                 Count=s2$Count, FE=round(s2$FoldEnrichment,2), p.adjust=signif(s2$p.adjust,3)), row.names=FALSE, right=FALSE)
cat("\n--- 시냅스/신경 KEGG term (유의성 무관) ---\n")
n <- r[grepl("synap|Synap|Long-term|Neuroact|Axon|Dopamin|GABA|Glutamat|Cholinerg|Seroton|Circadian|Calcium sig", r$Description),]
print(data.frame(Description=substr(n$Description,1,44), GeneRatio=n$GeneRatio, Count=n$Count,
                 FE=round(n$FoldEnrichment,2), pvalue=signif(n$pvalue,3), p.adjust=signif(n$p.adjust,3)), row.names=FALSE, right=FALSE)
}' 2>&1 | grep -v "^Warning\|built under\|^In addition\|경고메시지"

## ---- Liver KEGG ORA (EGM-up) -----------------------------------------------
cd "$ROOT/scratch" && Rscript -e '
suppressMessages(library(clusterProfiler))
s <- readRDS("Liver_sets.rds")
kk <- enrichKEGG(s$g, organism="hsa", universe=s$uni, pvalueCutoff=1, qvalueCutoff=1, minGSSize=5, maxGSSize=600)
saveRDS(kk, "KEGGup_Liver.rds"); r <- kk@result[order(kk@result$p.adjust),]
cat("Liver ORA-KEGG (UP): 총", nrow(r), " p.adj<0.05:", sum(r$p.adjust<0.05,na.rm=TRUE), "\n\n")
s2 <- r[r$p.adjust<0.05,][1:25,]
print(data.frame(ID=s2$ID, Description=substr(s2$Description,1,44), GeneRatio=s2$GeneRatio,
                 Count=s2$Count, FE=round(s2$FoldEnrichment,2), p.adjust=signif(s2$p.adjust,3)), row.names=FALSE, right=FALSE)
' 2>&1 | grep -v "^Warning\|built under\|^In addition\|경고메시지\|Reading KEGG"

## ---- copy KEGG objects, write FULL_KEGGup_Brain / Liver tables ---------------
SP="$ROOT/scratch"; FS="$ROOT/monkey/functional_study"
cp "$SP/KEGGup_Brain.rds" "$SP/KEGGup_Liver.rds" "$FS/" && Rscript -e '
FS <- file.path(Sys.getenv("ROOT"), "monkey/functional_study")
for(n in c("KEGGup_Brain","KEGGup_Liver")){
  r <- readRDS(file.path(FS, paste0(n,".rds")))@result
  r$ratio <- sapply(strsplit(r$GeneRatio,"/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  r <- r[order(r$p.adjust),]
  write.csv(r, file.path(FS,"GO_tables", paste0("FULL_",n,".csv")), row.names=FALSE)
  cat(sprintf("FULL_%s.csv  %d rows, %d sig\n", n, nrow(r), sum(r$p.adjust<0.05,na.rm=TRUE)))
}' 2>&1 | grep -v "^Warning\|built under" && echo "" && ls -la "$FS"/KEGGup_*.rds && echo "" && ls "$FS/GO_tables/"
