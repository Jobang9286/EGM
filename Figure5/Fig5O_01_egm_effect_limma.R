WORK_DIR    <- "public_monkey"   # holds merged/
setwd(WORK_DIR)

suppressMessages({library(data.table);library(limma);library(edgeR);library(AnnotationDbi)})
orgdb<-loadDb("../monkey/orgdb_build/org.Mfascicularis.eg.db/inst/extdata/org.Mfascicularis.eg.sqlite")
markers<-c("CDKN1A","CDKN2A","TP53","IL6","TNF","LMNB1","GLB1","SERPINE1","CDKN1B","MKI67",
 "COL1A1","COL3A1","COL1A2","IGF1","FOXO3","SIRT1","NFKB1","B2M","GDF15","TERT",
 "IL1B","CXCL8","MMP3","MMP9","IGFBP3","IGFBP7","TIMP1","HMGB1","NAMPT","KLOTHO","KL")
tis<-c("Brain","Liver","Kidney","Spleen")
allres<-list()
for(ti in tis){
  cnt<-fread(sprintf("merged/%s_counts.csv",ti));meta<-fread(sprintf("merged/%s_meta.csv",ti))
  keep<-meta$group %in% c("EGM","CTRL");m<-meta[keep];
  M<-as.matrix(cnt[,m$sample,with=FALSE]);rownames(M)<-cnt$gene_id
  grp<-factor(m$group,levels=c("CTRL","EGM"))
  d<-DGEList(M);d<-d[filterByExpr(d,group=grp),,keep.lib.sizes=FALSE];d<-calcNormFactors(d)
  des<-model.matrix(~grp);v<-voom(d,des);fit<-eBayes(lmFit(v,des))
  tt<-topTable(fit,coef=2,number=Inf,sort.by="none")
  tt$gene_id<-rownames(tt);tt$tissue<-ti
  ens<-sub("\\..*","",tt$gene_id)
  tt$symbol<-mapIds(orgdb,keys=ens,column="SYMBOL",keytype="ENSEMBL",multiVals="first")
  allres[[ti]]<-as.data.table(tt)
  cat(sprintf("== %s: nGenes=%d, adjP<.05=%d\n",ti,nrow(tt),sum(tt$adj.P.Val<0.05,na.rm=TRUE)))
}
res<-rbindlist(allres)
fwrite(res,"merged/egm_effect_all.csv")
fwrite(res[symbol %in% markers],"merged/egm_effect_markers.csv")
# fig4g plotdata: top 14 by |t| per tissue, symbol-named
pd<-res[!is.na(symbol)][,head(.SD[order(-abs(t))],14),by=tissue]
fwrite(pd,"merged/egm_fig4g_plotdata.csv")
cat("\n== MARKERS (EGM vs CTRL, logFC + adjP) ==\n")
mk<-res[symbol %in% markers][order(tissue,adj.P.Val)]
print(mk[,.(tissue,symbol,logFC=round(logFC,3),P=signif(P.Value,2),adjP=signif(adj.P.Val,2))])
