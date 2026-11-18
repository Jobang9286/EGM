############################################################################
# EGM functional study — reproduce  (crab-eating macaque, 4 organs)
# ------------------------------------------------------------------------
# RSEM raw counts  ->  DESeq2 (EGM vs Saline, within-organ)  ->
#   Phase 1: SASP / senescence / inflammaging marker heatmaps (VST z-score)
#   Phase 2: GO-BP / KEGG enrichment dot plots (EGM-restored functions)
#
# Run from THIS folder (setwd to the bundle dir). KEGG steps need internet;
# pre-computed enrichment objects are provided in data/cache/ to skip them.
############################################################################
suppressMessages({
  library(DESeq2); library(data.table); library(ComplexHeatmap); library(circlize)
  library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)
})

DATA   <- "data"                          # holds <Organ>/txi.rsem.counts.txt + mapfile.tsv, ortholog_map.csv
OUT    <- "outputs"; dir.create(OUT, showWarnings = FALSE)
organs <- c("Brain","Liver","Kidney","Spleen")

## ======================================================================
## 1. DESeq2 per organ (EGM vs Saline) + VST
## ======================================================================
om     <- fread(file.path(DATA,"ortholog_map.csv"))          # ENSMFAG -> human_symbol / entrez
om_sym <- om[!duplicated(gene_id), .(gene_id, sym = human_symbol)]
om2    <- om[!is.na(entrez) & entrez!=""][!duplicated(gene_id), .(gene_id, entrez = as.character(entrez))]

run_organ <- function(org){
  cts <- as.data.frame(fread(file.path(DATA,org,"txi.rsem.counts.txt")))
  rownames(cts) <- cts$gene_id; cts$gene_id <- NULL; cts <- round(as.matrix(cts))   # RSEM expected counts -> int
  map <- as.data.frame(fread(file.path(DATA,org,"mapfile.tsv")))
  rownames(map) <- map$sample; map <- map[colnames(cts),,drop=FALSE]
  map$grp <- factor(ifelse(map$TREATMENT=="treated","EGM","CTRL"), levels=c("CTRL","EGM"))
  dds <- DESeqDataSetFromMatrix(cts, map, ~grp)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]                # expressed filter
  dds <- DESeq(dds, quiet=TRUE)
  res <- as.data.frame(results(dds, contrast=c("grp","EGM","CTRL")))  # log2FC>0 = up by EGM
  res$gene_id <- rownames(res); res$tissue <- org
  list(res = merge(res, as.data.frame(om_sym), by="gene_id", all.x=TRUE),
       vst = assay(vst(dds, blind=FALSE)), map = map)
}
R       <- lapply(setNames(organs,organs), run_organ)
res_all <- lapply(R, `[[`, "res")
vst_all <- lapply(R, function(x) list(vst=x$vst, map=x$map))
saveRDS(res_all, file.path(OUT,"res_4organ.rds"))
saveRDS(vst_all, file.path(OUT,"vst_4organ.rds"))

## ======================================================================
## 2. Phase 1 — SASP / senescence / inflammaging marker heatmaps
## ======================================================================
panel2 <- list(
 Senescence  = c("CDKN1A","CDKN2A","CDKN2B","TP53","RB1","SERPINE1","GLB1","E2F1"),
 SASP_cyt    = c("IL6","IL1A","IL1B","IL18","IL15"),
 SASP_chemo  = c("CXCL1","CXCL2","CXCL8","CXCL9","CXCL10","CXCL12","CCL2","CCL5","CCL20"),
 SASP_prot   = c("MMP1","MMP2","MMP3","MMP9","MMP12","MMP13","MMP14","TIMP1","PLAU","PLAUR"),
 SASP_gf     = c("IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP7","GDF15","VEGFA","HGF","AREG","EREG","CSF1"),
 NFkB        = c("NFKB1","NFKB2","RELA","RELB","REL","NFKBIZ","NFKBIE","NFKBIB",
                 "TNF","TNFAIP3","TNFAIP6","TNFAIP2","TNFRSF1A","TNFRSF1B","LTB","PTGS2","NOS2","PTX3"),
 Adhesion    = c("ICAM1","VCAM1","SELE","SELP"),
 Inflammasome= c("NLRP3","CASP1","CASP4","PYCARD","AIM2","GSDMD"),
 Complement  = c("C1QA","C1QB","C1QC","C1R","C1S","C2","C3","C4A","CFB","CFH","SERPING1","C3AR1","C5AR1"),
 Myeloid     = c("AIF1","TYROBP","TREM2","CD68","CD74","GFAP","CTSS","HEXB","FCER1G","LAPTM5","ITGAM","FCGR3A","LY86"),
 IFN         = c("STAT1","NFKBIA","IFI6","IFI27","IFI44","IFI44L","ISG15","MX1","MX2","OAS1","OAS2","OASL",
                 "IFIT1","IFIT2","IFIT3","IFIH1","IRF7","DDX58","RSAD2","USP18","XAF1","BST2"),
 Aging       = c("LGALS3","LGALS3BP","B2M","VIM","CD44","S100B","S100A8","S100A9","HMGB1","EDA2R"))
mk <- data.table(sym=unlist(panel2), cat=rep(names(panel2), lengths(panel2)))
mk[cat %in% c("SASP_cyt","SASP_chemo","SASP_prot","SASP_gf"), cat := "SASP"]      # merge 4 SASP sub-axes
clev <- c("Senescence","SASP","NFkB","Adhesion","Inflammasome","Complement","Myeloid","IFN","Aging")
clab <- c("Senescence","SASP","NF-kB / TNF","Adhesion","Inflammasome","Complement",
          "Myeloid / Macrophage","Type 1 IFN / ISG","Aging / galectin")

# heatmap of markers significantly REDUCED by EGM (padj<.05 & log2FC<0), VST z-score, EGM|Saline split
make_hm <- function(org, keep=NULL, aging_label=NULL, rowgap=0.6){
  vs <- vst_all[[org]]$vst; map <- vst_all[[org]]$map; res <- res_all[[org]]
  d <- as.data.table(res)[!is.na(sym) & sym %in% mk$sym & padj<0.05 & log2FoldChange<0]
  d <- d[order(padj)][!duplicated(sym)]; d <- merge(d, mk, by="sym")
  if(!is.null(keep)) d <- d[cat %in% keep]
  d <- d[gene_id %in% rownames(vs)]
  lab <- clab; if(!is.null(aging_label)) lab[clev=="Aging"] <- aging_label
  d[, cat := factor(cat, levels=clev, labels=lab)]; d <- droplevels(d)
  d <- d[order(cat, log2FoldChange)]
  m <- vs[d$gene_id,,drop=FALSE]; rownames(m) <- d$sym; z <- t(scale(t(m)))
  eg <- rownames(map)[map$grp=="EGM"]; ct <- rownames(map)[map$grp=="CTRL"]
  z <- z[, c(eg,ct)]; colnames(z) <- c(paste0("EGM ",seq_along(eg)), paste0("Saline ",seq_along(ct)))
  csplit <- factor(c(rep("EGM",length(eg)), rep("Saline",length(ct))), levels=c("EGM","Saline"))
  Heatmap(z, name="Z-score", col=colorRamp2(c(-1.5,0,1.5), c("#2166AC","white","#B2182B")),
    cluster_columns=FALSE, cluster_rows=FALSE, column_split=csplit, column_gap=unit(0,"mm"),
    row_split=d$cat, row_gap=unit(rowgap,"mm"), row_title_rot=0,
    row_title_gp=gpar(fontsize=10,fontface="bold"),
    column_title=org, column_title_gp=gpar(fontsize=13,fontface="bold"),
    row_names_gp=gpar(fontsize=9,fontface="italic"), column_names_gp=gpar(fontsize=9.5),
    column_names_rot=45, rect_gp=gpar(col="white",lwd=0.6),
    heatmap_legend_param=list(title="Z-score", at=c(-2,0,2)))
}
save_hm <- function(ht,f,w,h){ png(file.path(OUT,f),width=w,height=h,units="in",res=200); draw(ht,heatmap_legend_side="right"); dev.off() }

# per-organ CORE categories (Brain = Type-I IFN/ISG + Reactive gliosis[=Aging relabelled])
save_hm(make_hm("Brain",  keep=c("IFN","Aging"), aging_label="Reactive Gliosis", rowgap=0.7), "heatmap_expr_Brain.png",  5.6, 4.6)
save_hm(make_hm("Liver",  keep=c("Senescence","SASP","Complement","Myeloid")),                "heatmap_expr_Liver.png",  6.0, 7.4)
save_hm(make_hm("Kidney", keep=c("Senescence","SASP","Inflammasome","Complement","Myeloid")), "heatmap_expr_Kidney.png", 6.0, 7.4)
save_hm(make_hm("Spleen", keep=c("Senescence","Inflammasome","Complement","Myeloid","IFN")),  "heatmap_expr_Spleen.png", 6.0, 8.8)

## ======================================================================
## 3. Phase 2 — GO-BP / KEGG enrichment dot plots
## ======================================================================
go_dir <- function(res, dir=c("up","down")){
  dir <- match.arg(dir); d <- merge(as.data.table(res), om2, by="gene_id")
  uni <- unique(d[!is.na(padj)]$entrez)
  g <- if(dir=="up") unique(d[!is.na(padj)&padj<0.05&log2FoldChange>0]$entrez)
       else          unique(d[!is.na(padj)&padj<0.05&log2FoldChange<0]$entrez)
  enrichGO(g, universe=uni, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="BP",
           pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.1, minGSSize=10, maxGSSize=500, readable=TRUE)
}
kegg_dir <- function(res, dir=c("up","down"), cutoff=1){
  dir <- match.arg(dir); d <- merge(as.data.table(res), om2, by="gene_id")
  uni <- unique(d[!is.na(padj)]$entrez)
  g <- if(dir=="up") unique(d[!is.na(padj)&padj<0.05&log2FoldChange>0]$entrez)
       else          unique(d[!is.na(padj)&padj<0.05&log2FoldChange<0]$entrez)
  enrichKEGG(g, organism="hsa", universe=uni, pvalueCutoff=cutoff, qvalueCutoff=cutoff, minGSSize=5, maxGSSize=600)
}
# reusable dot-plot: x=GeneRatio, size=Count, color=p.adjust (yellow->red default; blue for down)
mk_dot <- function(dd, title, sbreaks, cols=c("#B2182B","#EF6548","#FDAE61","#FFFFBF")){
  dd$ratio <- sapply(strsplit(dd$GeneRatio,"/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  dd <- dd[order(dd$ratio),]; dd$Description <- factor(dd$Description, levels=dd$Description)
  ggplot(dd, aes(ratio, Description, size=Count, color=p.adjust)) + geom_point() +
    scale_color_gradientn(colours=cols, name="p.adjust") +
    scale_size_continuous(range=c(2,6.5), breaks=sbreaks, name="Gene count") +
    guides(size=guide_legend(order=1), color=guide_colorbar(order=2)) +
    labs(x="Gene ratio", y=NULL, title=title) + theme_bw(base_size=12) +
    theme(panel.grid=element_blank(), axis.text.y=element_text(size=10.5,color="black"),
          plot.title=element_text(face="bold", size=12.5))
}
egB  <- go_dir(res_all$Brain ,"up");   dfB   <- as.data.frame(egB)
egL  <- go_dir(res_all$Liver ,"up");   dfL   <- as.data.frame(egL)
egK  <- go_dir(res_all$Kidney,"up");   dfK   <- as.data.frame(egK)
kkK  <- kegg_dir(res_all$Kidney,"up"); dallK <- as.data.frame(kkK)
egSd <- go_dir(res_all$Spleen,"down"); dfSd  <- as.data.frame(egSd)
kkSd <- kegg_dir(res_all$Spleen,"down",cutoff=0.1); dallSd <- as.data.frame(kkSd)

## --- Brain (UP): synaptic / neuronal / cognitive recovery
termsB <- c("cell morphogenesis involved in neuron differentiation","synapse organization",
 "modulation of chemical synaptic transmission","axonogenesis","regulation of neuron projection development",
 "regulation of membrane potential","cognition","synapse assembly","learning or memory",
 "neurotransmitter transport","positive regulation of neurogenesis","action potential",
 "neuron migration","synaptic transmission, glutamatergic","dendritic spine development")
ggsave(file.path(OUT,"dotplot_GOup_Brain.png"),
  mk_dot(dfB[dfB$Description %in% termsB,], "Brain - aging-declined functions restored by EGM", c(50,100,150)),
  width=8.0, height=5.4, dpi=200)

## --- Liver (UP): diversified hepatocyte functions
termsL <- c("fatty acid beta-oxidation","cholesterol metabolic process","bile acid biosynthetic process",
 "cellular detoxification","xenobiotic metabolic process","cellular oxidant detoxification","urea cycle",
 "gluconeogenesis","glucose metabolic process","plasma lipoprotein particle organization","cholesterol transport",
 "heme metabolic process","ketone metabolic process","amino acid catabolic process")
ggsave(file.path(OUT,"dotplot_GOup_Liver.png"),
  mk_dot(dfL[dfL$Description %in% termsL,], "Liver - aging-declined functions restored by EGM", c(10,30,60)),
  width=8.4, height=5.2, dpi=200)

## --- Kidney (UP): GO + KEGG mix, kidney-specific favourable (P450 removed)
goKp <- c("sodium ion transmembrane transport","proton transmembrane transport","regulation of blood pressure","cellular respiration")
keKp <- c("Aldosterone-regulated sodium reabsorption","Proximal tubule bicarbonate reclamation","Mineral absorption",
          "Endocrine and other factor-regulated calcium reabsorption","Oxidative phosphorylation","Citrate cycle (TCA cycle)","PPAR signaling pathway")
g <- dfK[dfK$Description %in% goKp, c("Description","GeneRatio","Count","p.adjust")]; g$src <- "GO"
k <- dallK[dallK$Description %in% keKp, c("Description","GeneRatio","Count","p.adjust")]; k$src <- "KEGG"
comb <- rbind(g,k)
comb$Description <- sub("Endocrine and other factor-regulated calcium reabsorption","Endocrine calcium reabsorption",comb$Description)
comb$Description <- paste0(comb$Description, "  (", comb$src, ")")
ggsave(file.path(OUT,"dotplot_Kidney_GO_KEGG.png"),
  mk_dot(comb, "Kidney - aging-declined renal functions restored by EGM", c(15,45,80)), width=8.6, height=5.0, dpi=200)

## --- Spleen (DOWN): immune / inflammatory activation reduced (BLUE)
goSp <- c("adaptive immune response","positive regulation of cytokine production","T cell activation",
          "leukocyte mediated immunity","regulation of innate immune response","cytokine-mediated signaling pathway","leukocyte migration")
keSp <- c("Antigen processing and presentation","Phagocytosis","NOD-like receptor signaling pathway",
          "Natural killer cell mediated cytotoxicity","Th17 cell differentiation")
gS <- dfSd[dfSd$Description %in% goSp, c("Description","GeneRatio","Count","p.adjust")]; gS$src <- "GO"
kS <- dallSd[dallSd$Description %in% keSp, c("Description","GeneRatio","Count","p.adjust")]; kS$src <- "KEGG"
combS <- rbind(gS,kS); combS$Description <- paste0(combS$Description, "  (", combS$src, ")")
ggsave(file.path(OUT,"dotplot_Spleen_DOWN_GO_KEGG.png"),
  mk_dot(combS, "Spleen - aging immune/inflammatory activation reduced by EGM", c(50,110,180),
         cols=c("#08306B","#2171B5","#6BAED6","#C6DBEF")), width=8.9, height=5.0, dpi=200)

cat("DONE. Figures + rds written to", normalizePath(OUT), "\n")
