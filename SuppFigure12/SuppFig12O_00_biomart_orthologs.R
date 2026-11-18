suppressMessages(library(dplyr))

## ---- BioMart query (M. fascicularis -> human orthologs) ---------------------------------------------
xml <- '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6"><Dataset name="mfascicularis_gene_ensembl" interface="default"><Attribute name="ensembl_gene_id"/><Attribute name="external_gene_name"/><Attribute name="hsapiens_homolog_associated_gene_name"/><Attribute name="hsapiens_homolog_orthology_type"/></Dataset></Query>'

## ---- query the Ensembl mirrors -----------------------------------------------------------------------
get_orth <- function(host) {
  uu <- paste0("https://", host, "/biomart/martservice?query=", URLencode(xml, reserved=TRUE))
  tryCatch({
    x <- read.delim(url(uu), stringsAsFactors=FALSE, colClasses="character")
    if (ncol(x) != 4) stop("bad response: ", names(x)[1]) else x
  }, error=function(e) conditionMessage(e))
}
for (h in c("asia.ensembl.org","useast.ensembl.org","www.ensembl.org")) {
  orth <- get_orth(h); cat(h, "->", if (is.data.frame(orth)) nrow(orth) else substr(orth,1,120), "\n")
  if (is.data.frame(orth)) break
}

## ---- column names ------------------------------------------------------------------------------------
names(orth) <- c("gene_id","mfa_symbol","hs_symbol","orth_type")

## ---- save --------------------------------------------------------------------------------------------
kbase <- "EGM_tissue/oncogene_analysis/Monkey"
kout <- file.path(kbase, "oncogene_heatmap"); dir.create(kout, showWarnings=FALSE)
stopifnot(exists("orth"), all(c("gene_id","mfa_symbol","hs_symbol","orth_type") %in% names(orth)))
write.csv(orth %>% filter(hs_symbol != ""), file.path(kout, "ensembl_mfascicularis_human_orthologs.csv"), row.names=FALSE)
