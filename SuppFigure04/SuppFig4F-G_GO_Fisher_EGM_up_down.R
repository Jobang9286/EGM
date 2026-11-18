library(limma)
mm.kegg <- getGeneKEGGLinks(species = "Mm")

getGeneKEGG

load("39_fibroblasts/RNA_seq/htseq/KEGG_2016_disease_aging_activation_lower_case_for_PAGODA.RData")


aa <- kegg_aging.env$`2-Oxocarboxylic acid metabolism`

# Make it into a data frame
genes<- as.data.frame(resSig)

summary(resSig)

genes

genes1 <- toupper(rownames(genes))

rownames(genes) <- genes1

genes

# Upload kegg pathways (the .csv version)
gene.sets <- read.table("39_fibroblasts/RNA_seq/htseq/m5.go.bp.v2024.1.Mm.symbols.txt", sep="\t", header = F, quote="",row.names = 1)
gene.set.data <- as.data.frame(t(gene.sets))
gene.sets

gene.sets <- read.table("41_RNA_seq_test/gene_set_library_crisp.txt", sep="\t", header = F, quote="",row.names = 1)
gene.set.data <- as.data.frame(t(gene.sets))
gene.sets
########################################################################
# Extract genes that are DOWN in Old
deg.down <- genes[which(genes$pvalue < 0.05 & genes$log2FoldChange > 0.5),]

deg.down
# Make a data frame with 3 columns for p-value, FDR-corrected p-values , and the string of genes within the pathway
my.data.down <- as.data.frame(matrix(0, ncol = 3, nrow = ncol(gene.set.data)))
colnames(my.data.down) <- c("p.values", "FDR.BH", "genes")
rownames(my.data.down) <- colnames(gene.set.data)

# In a loop go though each and one of the KEGG pathways, and determine
# 1) how many genes in the KEGG pathway are found in my gene list (background) (kegg_genes)
# 2) how many genes in the KEGG pathway are found DEG in my data set (deg_kegg_genes)
# 3) how many genes in the KEGG pathway are found in my gene list but are not DEG in my data set (deg_non_kegg_genes)
# 4) how many DEG genes are not found in the KEGG pathway (non_deg_kegg_genes)
# 5) how many genes are not DEG nor found in the KEGG pathway (non_deg_non_kegg_genes)
# Use the numbers to perform fischer's t-test, and save the p-value in the data frame
# Finally correct for FDR using benjimini_hochberg

for(i in 1:ncol(gene.set.data)){
  kegg_genes <- nrow(genes[which(rownames(genes) %in% gene.set.data[,i]),])
  deg_kegg_genes <- nrow(genes[which(rownames(deg.down) %in% gene.set.data[,i]),])
  deg_kegg_genes.names <- as.data.frame(row.names(deg.down[which(rownames(deg.down) %in% gene.set.data[,i]),]))
  deg_non_kegg_genes <- nrow(deg.down) - deg_kegg_genes
  non_deg_kegg_genes <- kegg_genes - deg_kegg_genes
  non_deg_non_kegg_genes <- (nrow(genes) -nrow(deg.down)) - kegg_genes
  ftest <- fisher.test(matrix(c(deg_kegg_genes,deg_non_kegg_genes,non_deg_kegg_genes,non_deg_non_kegg_genes),nrow=2,ncol=2),alternative="greater")
  my.data.down$p.values[i] <- ftest$p.value
  my.data.down$genes[i] <- paste(deg_kegg_genes.names[1:nrow(deg_kegg_genes.names),], collapse=", ")
}
my.data.down$FDR.BH <- p.adjust(my.data.down$p.values, "BH")
# Order by FDR adjusted p-values
my.data.downOrdered <- my.data.down[order(my.data.down$FDR.BH),]
my.data.downOrdered
# Write a table with the final results
write.csv(my.data.downOrdered, file="39_fibroblasts/RNA_seq/final_EGM_vs_no_EGM/GO_EGM_UP.csv")

#####################################################################################################
# Extract genes that are UP in old
deg.up <- genes[which(genes$pvalue < 0.05 & genes$log2FoldChange < -0.5),] #633 UP in Old
deg.up

write.csv(deg.up, "./deg_up_3f.csv")
# Make a data frame with 3 columns for p-value, FDR-corrected p-values , and the string of genes within the pathway
my.data.up <- as.data.frame(matrix(0, ncol = 3, nrow = ncol(gene.set.data)))
colnames(my.data.up) <- c("p.values", "FDR.BH", "genes")
rownames(my.data.up) <- colnames(gene.set.data)

# In a loop go though each and one of the KEGG pathways, and determine
# 1) how many genes in the KEGG pathway are found in my gene list (background) (kegg_genes)
# 2) how many genes in the KEGG pathway are found DEG in my data set (deg_kegg_genes)
# 3) how many genes in the KEGG pathway are found in my gene list but are not DEG in my data set (deg_non_kegg_genes)
# 4) how many DEG genes are not found in the KEGG pathway (non_deg_kegg_genes)
# 5) how many genes are not DEG nor found in the KEGG pathway (non_deg_non_kegg_genes)
# Use the numbers to perform fischer's t-test, and save the p-value in the data frame
# Finally correct for FDR using benjimini_hochberg

for(i in 1:ncol(gene.set.data)){
  kegg_genes <- nrow(genes[which(rownames(genes) %in% gene.set.data[,i]),])
  deg_kegg_genes <- nrow(genes[which(rownames(deg.up) %in% gene.set.data[,i]),])
  deg_kegg_genes.names <- as.data.frame(row.names(deg.up[which(rownames(deg.up) %in% gene.set.data[,i]),]))
  deg_non_kegg_genes <- nrow(deg.up) - deg_kegg_genes
  non_deg_kegg_genes <- kegg_genes - deg_kegg_genes
  non_deg_non_kegg_genes <- (nrow(genes) -nrow(deg.up)) - kegg_genes
  ftest <- fisher.test(matrix(c(deg_kegg_genes,deg_non_kegg_genes,non_deg_kegg_genes,non_deg_non_kegg_genes),nrow=2,ncol=2),alternative="greater")
  my.data.up$p.values[i] <- ftest$p.value
  my.data.up$genes[i] <- paste(deg_kegg_genes.names[1:nrow(deg_kegg_genes.names),], collapse=", ")
}
my.data.up$FDR.BH <- p.adjust(my.data.up$p.values, "BH")
# Order by FDR adjusted p-values
my.data.upOrdered <- my.data.up[order(my.data.up$FDR.BH),]
my.data.upOrdered
# Write a table with the final results
write.csv(my.data.upOrdered, file="39_fibroblasts/RNA_seq/final_EGM_vs_no_EGM/GO_EGM_DOWN.csv")
