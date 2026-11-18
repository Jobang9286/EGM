Sys.setenv(RETICULATE_PYTHON = file.path(getwd(), ".venv/bin/python"))
suppressMessages({library(tAge); library(Biobase)})
CK <- "clocks"

load_tissue <- function(twf, ogf){
  tw <- read.table(twf, header=TRUE, row.names=1, sep=if(grepl("\\.csv$",twf)) "," else "\t", check.names=FALSE)
  og <- read.table(ogf, header=TRUE, row.names=1, sep=if(grepl("\\.csv$",ogf)) "," else "\t", check.names=FALSE)
  colnames(tw) <- sub("\\.bam$","",colnames(tw)); colnames(og) <- sub("\\.bam$","",colnames(og))
  hg <- og[, grep("HGPS", colnames(og)), drop=FALSE]; g <- intersect(rownames(tw), rownames(og))
  cbind(tw[g,,drop=FALSE], hg[g,,drop=FALSE])
}
predict_tissue <- function(tissue, twf, ogf, model, scaling="scaled_diff"){
  M   <- load_tissue(twf, ogf)
  age <- as.numeric(sub(".*WT_(\\d+)[Mm]_.*","\\1", colnames(M)))
  trt <- ifelse(grepl("HGPS_control",colnames(M)),"HGPS_ctrl",
         ifelse(grepl("HGPS_EGM",colnames(M)),"HGPS_EGM",paste0("WT_",age,"M")))
  meta <- data.frame(sample=colnames(M), treatment=trt, age=age, row.names=colnames(M))
  eset <- make_ExpressionSet(as.data.frame(M), meta)
  pe   <- tAge_preprocessing(eset, species="mouse", gene_mapping_type="Gene.Symbol",
                             control_group_column="treatment", control_group_label="WT_1M", verbose=FALSE)
  pr   <- as.data.frame(predict_tAge(pe, setNames(list(model), scaling), species="mouse", mode="EN"))
  col  <- grep("_tAge$", colnames(pr), value=TRUE)[1]
  data.frame(sample=rownames(pr), tissue=tissue, treatment=meta[rownames(pr),"treatment"],
             age=meta[rownames(pr),"age"], tAge=pr[[col]], row.names=NULL)
}
cfg <- list(
  Brain =c("data/counts_brain_time.csv","data/counts_brain.txt", file.path(CK,"EN_Mortality_Rodents_Multitissue_scaleddiff.pkl")),
  Kidney=c("data/counts_kidney_time.csv","data/counts_kidney.txt",file.path(CK,"EN_Mortality_Rodents_Multitissue_scaleddiff.pkl")),
  Liver =c("data/counts_liver_for_time.csv","data/counts_liver.csv",file.path(CK,"EN_Mortality_Rodents_Liver_scaleddiff.pkl")))
out <- do.call(rbind, lapply(names(cfg), function(t) predict_tissue(t, cfg[[t]][1], cfg[[t]][2], cfg[[t]][3])))
dir.create("outputs", showWarnings=FALSE)
write.csv(out, "outputs/tage_mortality_persample.csv", row.names=FALSE)
lev <- c(paste0("WT_",c(1,3,6,12,18,24,27),"M"),"HGPS_ctrl","HGPS_EGM")
for (t in names(cfg)){ dd <- out[out$tissue==t,]; wt <- dd[!is.na(dd$age),]
  cat("\n=====",t,"(Mortality tAge log10HR) =====\n")
  agg <- aggregate(tAge~treatment, dd, function(x) round(mean(x),3)); print(agg[order(match(agg$treatment,lev)),], row.names=FALSE)
  cat("WT cor(mortAge,age)=", round(cor(wt$tAge,wt$age),3),
      "| Saline", round(mean(dd$tAge[dd$treatment=="HGPS_ctrl"]),3),
      "EGM", round(mean(dd$tAge[dd$treatment=="HGPS_EGM"]),3),
      "delta", round(mean(dd$tAge[dd$treatment=="HGPS_EGM"])-mean(dd$tAge[dd$treatment=="HGPS_ctrl"]),3),
      "p", signif(t.test(dd$tAge[dd$treatment=="HGPS_EGM"], dd$tAge[dd$treatment=="HGPS_ctrl"])$p.value,3),"\n") }
cat("\nsaved outputs/tage_mortality_persample.csv\n")
