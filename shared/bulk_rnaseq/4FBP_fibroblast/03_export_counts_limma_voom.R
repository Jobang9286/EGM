suppressMessages({library(data.table); library(limma); library(edgeR)})
TG <- "fib_target"
sam <- c("4FP-1","4FP-2","4FP-3","CP-1","CP-2","CP-3")
rd <- function(s) fread(file.path(TG, paste0(s,".genes.results")), select=c("gene_id","expected_count","TPM"))
L <- lapply(sam, rd); g <- L[[1]]$gene_id
stopifnot(all(sapply(L, function(x) identical(x$gene_id, g))))
C <- sapply(L, function(x) round(x$expected_count)); T <- sapply(L, function(x) x$TPM)
colnames(C) <- colnames(T) <- sam; rownames(C) <- rownames(T) <- g

fwrite(data.table(gene=g, as.data.table(C)), "4FP_CP_counts_symbol.csv")
fwrite(data.table(gene=g, as.data.table(round(T,4))), "4FP_CP_tpm_symbol.csv")
meta <- data.table(sample=sam, group=ifelse(grepl("^4FP",sam),"4FP","CP"),
  flowcell="HWLLKDSX7", lane=1L, platform="NovaSeq", read_len=101L, layout="paired",
  strandedness="reverse (dUTP)", forward_prob=0,
  reference="GRCm39 + GENCODE vM38, RSEM symbol-level",
  pipeline="trim_galore --poly_g -> STAR --quantMode TranscriptomeSAM -> RSEM --forward-prob 0")
fwrite(meta, "4FP_CP_meta.csv")

grp <- factor(meta$group, levels=c("CP","4FP"))
d <- DGEList(C); keep <- filterByExpr(d, group=grp); d <- d[keep,,keep.lib.sizes=FALSE]; d <- calcNormFactors(d)
des <- model.matrix(~grp); v <- voom(d, des); f <- eBayes(lmFit(v, des))
tt <- as.data.table(topTable(f, coef=2, number=Inf, sort.by="P"), keep.rownames="gene")
fwrite(tt, "4FP_CP_DE_limma_voom.csv")
fwrite(data.table(gene=rownames(v$E), as.data.table(round(v$E,4))), "4FP_CP_logCPM_voom.csv")

cat("genes (all):", nrow(C), "  genes after filterByExpr:", nrow(d), "\n")
cat("FDR<0.05:", sum(tt$adj.P.Val<0.05), "  library sizes (M):", paste(round(colSums(C)/1e6,1), collapse=" "), "\n")
