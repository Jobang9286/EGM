# ---- packages ----
suppressMessages({library(DESeq2); library(dplyr)})
# ---- BioMart REST query ----
mbase <- "EGM_tissue/oncogene_analysis/Mouse"
xml_mm <- '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6"><Dataset name="mmusculus_gene_ensembl" interface="default"><Attribute name="ensembl_gene_id"/><Attribute name="external_gene_name"/><Attribute name="hsapiens_homolog_associated_gene_name"/><Attribute name="hsapiens_homolog_orthology_type"/></Dataset></Query>'
orth_mm <- NULL
for (h in c("asia.ensembl.org","www.ensembl.org","useast.ensembl.org")) {
  uu <- paste0("https://", h, "/biomart/martservice?query=", URLencode(xml_mm, reserved=TRUE))
  x <- tryCatch(read.delim(url(uu), stringsAsFactors=FALSE, colClasses="character"), error=function(e) NULL)
  if (!is.null(x) && ncol(x)==4) { orth_mm <- x; cat("from", h, ":", nrow(x), "rows\n"); break }
}
names(orth_mm) <- c("ens","mm_symbol","hs_symbol","orth_type")
# ---- output folder ----
outdir <- file.path(mbase, "oncogene_heatmap"); dir.create(outdir, showWarnings=FALSE)
# ---- save ----
write.csv(orth_mm %>% filter(hs_symbol != ""), file.path(outdir, "ensembl_mouse_human_orthologs.csv"), row.names=FALSE)
