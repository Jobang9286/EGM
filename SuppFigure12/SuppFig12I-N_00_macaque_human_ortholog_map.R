suppressMessages({library(gprofiler2); library(org.Hs.eg.db); library(AnnotationDbi)})
base <- "monkey"
uni <- read.csv(file.path(base,"Brain","DESeq2_results_annotated.csv"))
gid_all <- unique(uni$gene_id); cat("universe gene_ids:", length(gid_all), "\n")
# chunked gorth
chunks <- split(gid_all, ceiling(seq_along(gid_all)/8000))
t0 <- Sys.time()
gl <- lapply(chunks, function(cc)
  tryCatch(gorth(cc,"mfascicularis","hsapiens",mthreshold=1,filter_na=TRUE),
           error=function(e) NULL))
g <- do.call(rbind, gl)
cat("gorth elapsed:", round(difftime(Sys.time(),t0,units="secs"),1),"s  rows:", nrow(g), "\n")
conv <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, unique(g$ortholog_ensg),c("ENTREZID","SYMBOL"),"ENSEMBL"))
g2 <- merge(g, conv, by.x="ortholog_ensg", by.y="ENSEMBL")
master <- unique(data.frame(gene_id=g2$input, human_ensg=g2$ortholog_ensg,
                            human_symbol=g2$SYMBOL, entrez=g2$ENTREZID))
master <- master[!is.na(master$entrez),]
write.csv(master, file.path(base,"ortholog_map.csv"), row.names=FALSE)
cat("master map: ", nrow(master), "rows,", length(unique(master$gene_id)),
    "unique macaque genes mapped (", round(100*length(unique(master$gene_id))/length(gid_all),1), "% of universe )\n")
