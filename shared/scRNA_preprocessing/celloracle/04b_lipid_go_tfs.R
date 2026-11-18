suppressPackageStartupMessages({library(org.Mm.eg.db); library(GO.db); library(AnnotationDbi)})
roots <- c(lipid_metabolic="GO:0006629", reg_lipid_metabolic="GO:0019216",
           lipid_localization="GO:0010876", reg_lipid_localization="GO:1905952",
           lipid_homeostasis="GO:0055088", lipid_storage="GO:0019915", reg_lipid_storage="GO:0010883",
           fat_cell_diff="GO:0045444", reg_fat_cell_diff="GO:0045598")
res <- AnnotationDbi::select(org.Mm.eg.db, keys=unname(roots), keytype="GOALL", columns=c("SYMBOL"))
res$category <- names(roots)[match(res$GOALL, roots)]
agg <- aggregate(category ~ SYMBOL, data=unique(res[,c("SYMBOL","category")]), FUN=function(x) paste(sort(unique(x)), collapse="|"))
write.csv(agg, "lipid_go_genes.csv", row.names=FALSE)
cat("genes with lipid GO:", nrow(agg), "\n")
