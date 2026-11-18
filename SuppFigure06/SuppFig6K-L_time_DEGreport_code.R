root <- getwd()

library(SummarizedExperiment)
library(ggplot2)
library(DESeq2)
library(DEGreport)
library(ggsci)
library(dplyr)


setwd(file.path(root, "41_RNA_seq_test/"))

liv_aging_mat <- as.data.frame(read.csv(file.path(root, "39_fibroblasts/2_RNA_seq_tissue/Liver_final/counts_liver_for_time.csv")))
row.names(liv_aging_mat) <- liv_aging_mat[,"Geneid"]
liv_aging_mat <- as.data.frame(liv_aging_mat[,-1])

liv_aging_col1 <- colnames(liv_aging_mat)
liv_aging_col1

condition = c(rep('bbb',7),  rep('aaa',6), rep('ccc',7), rep('aaa',6), rep('bbb',4))
condition


exp.data <- data.frame(
  row.names=colnames(liv_aging_mat),
  condition = condition)

exp.data

liver_aging_dds <- DESeqDataSetFromMatrix(countData = liv_aging_mat,
                                    colData = exp.data,
                                    design = ~ condition)

liver_aging_dds <- DESeq(liver_aging_dds)

#Ift70a1, 

genelist_liver <- c("Sccpdh", "Coq8a", "Spc25", "Klhl23", "Ift70a1", "Capn3", "Abhd18", "Postn", "Pdgfc", "Glt28d2", "Tdrkh", "Cgn", "Hsd3b5", "Slc16a1", "Prok1", "Sort1", "Gm11832", "Gm12679", "Mup3", "Crot", "Spon2", "Gm32304", "Ugt2b1", "Arhgap24", "Oas3", "Gm34333", "Smo", "Fabp1", "Mug1", "Gm10319", "Gm44478", "Slco1a1", "Ndn", "Hsd3b7", "Eps8l2", "Gm2814", "4833411C07Rik", "Trmt9b", "Cyp4v3", "Acsl1", "Sh2d4a", "Ces2b", "Ces2c", "Ces3a", "Cd276", "Gm32017", "Gclc", "Cyp8b1", "Gask1a", "Cdk1", "Rab36", "Tjp3", "Gm36041", "Cpeb4", "Gm16170", "Mrnip", "Ulk2", "Slc16a11", "Abca8a", "Lrrn3", "Serpina9", "Gm36264", "Gm32063", "Gm40935", "S1pr3", "Appl1", "Fermt2", "Chrna2", "Ccn4", "Adck5", "Slc39a4", "Bik", "Gm38456", "Cldn1", "Gm19522", "Tmem45a", "Pde9a", "Zfp101", "Rnf125", "Gm51400", "Gm36629", "Gm4952", "Keg1", "C330002G04Rik", "Loxl4", "Elovl3", "Xk", "Slc9a7", "Mbnl3", "Bhlhb9", "Tceal1")

genelist_liver2 <- c("Sccpdh", "Coq8a", "Spc25", "Klhl23", "Ift70a1", "Capn3", "Abhd18", "Postn", "Pdgfc", "Glt28d2", "Tdrkh", "Cgn", "Hsd3b5", "Slc16a1", "Prok1", "Sort1", "Gm11832", "Gm12679", "Mup3", "Crot", "Spon2", "Gm32304", "Ugt2b1", "Arhgap24", "Oas3", "Gm34333", "Smo", "Fabp1", "Mug1", "Gm10319", "Gm44478", "Slco1a1", "Ndn", "Hsd3b7", "Eps8l2", "Gm2814", "4833411C07Rik", "Trmt9b", "Cyp4v3", "Acsl1", "Sh2d4a", "Ces2b", "Ces2c", "Ces3a", "Cd276", "Gm32017", "Gclc", "Cyp8b1", "Gask1a", "Cdk1", "Rab36", "Tjp3", "Gm36041", "Cpeb4", "Gm16170", "Mrnip", "Ulk2", "Slc16a11", "Abca8a", "Lrrn3", "Serpina9", "Gm36264", "Gm32063", "Gm40935", "S1pr3", "Appl1", "Fermt2", "Chrna2", "Ccn4", "Adck5", "Slc39a4", "Bik", "Gm38456", "Cldn1", "Gm19522", "Tmem45a", "Pde9a", "Zfp101", "Rnf125", "Gm51400", "Gm36629", "Gm4952", "Keg1", "C330002G04Rik", "Loxl4", "Elovl3", "Xk", "Slc9a7", "Mbnl3", "Bhlhb9", "Tceal1")


rld_liver <- rlogTransformation(liver_aging_dds)

rld_mat_liver <- assay(rld_liver)
rld_mat_liver

cluster_liver_rlog <- rld_mat_liver[genelist_liver2, ]

clusters <- degPatterns(cluster_liver_rlog, metadata = colData(liver_aging_dds), time = "condition", minc = 0, plot = T, cutoff = 0.7)

library(ggplot2)
library(viridis)

clusters <- degPatterns(cluster_liver_rlog, metadata = colData(liver_dds), time = "condition", minc = 0, plot = T, cutoff = 0.7) +
  scale_color_viridis_d(option = "viridis") + # "viridis", "magma", "plasma" 중 선택
  theme_minimal() +
  theme(legend.position = "bottom") # 범례 위치 조정

clusters <- degPatterns(cluster_liver_rlog, metadata = colData(liver_dds), time = "condition", minc = 0, plot = T, cutoff = 0.7)
filtclusters<-clusters$normalized
filtclusters$cluster
upc = c("1")
downc = c("2")

filtclusters1<-filtclusters[filtclusters$cluster %in% upc,]
filtclusters1[filtclusters1$cluster %in% upc,"cluster"] = 1
p1 <- degPlotCluster(filtclusters1,time="condition",color="condition")+geom_line(aes_string(group="genes"),alpha=0.05)+scale_color_nejm()+theme_bw()

filtclusters2<-filtclusters[filtclusters$cluster %in% downc,]
filtclusters2[filtclusters2$cluster %in% downc,"cluster"] = 2
p2 <- degPlotCluster(filtclusters2,time="condition",color="condition", min_genes = 0)+geom_line(aes_string(group="genes"),alpha=0.05)+scale_color_nejm()+theme_bw()

p <- plot_grid(p1, p2)
p

setwd(file.path(root, "39_fibroblasts/2_RNA_seq_tissue/final/"))

ggsave("DEGreport_liver_1216.pdf", plot = p, height = 2, width = 12, units = "in", dpi = 300)

############################################################################################################################

setwd(file.path(root, "41_RNA_seq_test/"))

br_mat <- as.data.frame(read.csv(file.path(root, "39_fibroblasts/2_RNA_seq_tissue/Brain_final/bam/counts_brain_time.csv")))
row.names(br_mat) <- br_mat[,"Geneid"]
br_mat <- as.data.frame(br_mat[,-1])

br_col1 <- colnames(br_mat)
br_col1

condition = c(rep('bbb',8),  rep('aaa',6), rep('ccc',8), rep('aaa',6), rep('bbb',4))
condition


exp.data <- data.frame(
  row.names=colnames(br_mat),
  condition = condition)

exp.data

br_dds <- DESeqDataSetFromMatrix(countData = br_mat,
                                    colData = exp.data,
                                    design = ~ condition)

br_dds <- DESeq(br_dds)

#Ift70a1, 

genelist_brain <- c("Smap1", "Ptp4a1", "Lman2l", "Tmeff2", "Pgap1", "Dsel", "Lhx9", "Kcnt2", "B3galt2", "Mpzl1", "Fam107b", "Frmd4a", "Sec61a2", "Galnt13", "Nr4a2", "Scn3a", "B3galt1", "Gm52555", "Agps", "Ube2e3", "Mir670hg", "Apip", "Ryr3", "Pdia3", "Il1bos", "Dnaaf9", "Pak5", "Rgs19", "Cyp7b1", "Hltf", "Skil", "Ttc14", "C230034O21Rik", "Syt11", "Hax1", "Them4", "Nras", "Capza1", "St7l", "Eeig2", "Plppr4", "Plppr5", "Fnbp1l", "Ugt8a", "Elovl6", "Pkn2", "Penk", "Ubxn2b", "Ccnc", "Tmeff1", "Tnc", "Nfib", "Tek", "Nrdc", "Slc6a9", "Marcksl1", "Crocc", "Cfap69", "Slc4a2", "Fbxl5", "Kcnip4", "Tmem33", "Gabra2", "Miat", "Slc29a4", "Nptx2", "Pds5b", "St7", "Rbm28", "Tra2a", "Chn2", "Lrrtm1", "Wnt7a", "Tafa1", "Grm7", "Marchf8", "Slco1c1", "Mia", "Akt2", "Gm45552", "Mkrn3", "Sh3gl3", "Gm36584", "Il16", "Rab30", "Rras2", "Gm39080", "Btbd16", "Agpat5", "Thap1", "Rnf122", "Enpp6", "Tll1", "Ints10", "Hhip", "Cbln1", "Kirrel3", "Crtam", "Ncam1", "AI593442", "Aph1b", "Stt3b", "Lztfl1", "Ccr5", "Rmnd1", "Ust", "Fabp7", "Gamt", "Gal3st1", "Emid1", "Vstm2a", "Ccnjl", "Gria1", "Trim17", "Gjc2", "Mettl16", "Tlcd3a", "Dhrs13", "Ypel2", "Tmem100", "Hpcal1", "Adam17", "Rnf144a", "Lamb1", "Cog5", "Arhgap5", "Prps1l3", "Dbpht2", "Akap5", "Fut8", "Eif2s1", "Eml1", "Dip2c", "Drd1", "Exoc3", "Rgs7bp", "Elovl7", "Rab3c", "Top2b", "Bmpr1a", "Gm20642", "Fermt2", "Exoc5", "Rem2", "Cmtm5", "Xkr6", "Pcdh17", "Pcdh20", "Sema5a", "Gm31462", "Desi1", "Shisal1", "Tuba1a", "Shisa9", "2510009E07Rik", "Nectin3", "Ift57", "Cadm2", "Slc5a3", "Pde10a", "Syngr3", "Srpk1", "Mdga1", "Zfp871", "Eml4", "Prkce", "Nrxn1", "Rnf138", "Dtna", "Gpr17", "Kctd16", "Atg12", "Ap3s1", "Sema6a", "Csnk1g3", "Dok6", "Tmx3", "Rfx3", "Apln", "Igsf1", "Fmr1", "Atp2b3", "Plp1", "Trpc5", "Cnksr2", "Glra2", "Kdm5d", "Uty")

genelist_brain2 <- c("Smap1","Ptp4a1","Lman2l","Tmeff2","Pgap1","Dsel","Lhx9","Kcnt2","B3galt2","Mpzl1","Fam107b","Frmd4a","Sec61a2","Galnt13","Nr4a2","Scn3a","B3galt1","Gm52555","Agps","Ube2e3","Mir670hg","Apip","Ryr3","Pdia3","Il1bos","Dnaaf9","Pak5","Rgs19","Cyp7b1","Hltf","Skil","Ttc14","C230034O21Rik","Syt11","Hax1","Them4","Nras","Capza1","St7l","Eeig2","Plppr4","Plppr5","Fnbp1l","Ugt8a","Elovl6","Pkn2","Penk","Ubxn2b","Ccnc","Tmeff1","Tnc","Nfib","Tek","Nrdc","Slc6a9","Marcksl1","Crocc","Cfap69","Slc4a2","Fbxl5","Kcnip4","Tmem33","Gabra2","Miat","Slc29a4","Nptx2","Pds5b","St7","Rbm28","Tra2a","Chn2","Lrrtm1","Wnt7a","Tafa1","Grm7","Marchf8","Slco1c1","Mia","Akt2","Gm45552","Mkrn3","Sh3gl3","Gm36584","Il16","Rab30","Rras2","Gm39080","Btbd16","Agpat5","Thap1","Rnf122","Enpp6","Tll1","Ints10","Hhip","Cbln1","Kirrel3","Crtam","Ncam1","AI593442","Aph1b","Stt3b","Lztfl1","Ccr5","Rmnd1","Ust","Fabp7","Gamt","Gal3st1","Emid1","Vstm2a","Ccnjl","Gria1","Trim17","Gjc2","Mettl16","Tlcd3a","Dhrs13","Ypel2","Tmem100","Hpcal1","Adam17","Rnf144a","Lamb1","Cog5","Arhgap5","Prps1l3","Dbpht2","Akap5","Fut8","Eif2s1","Eml1","Dip2c","Drd1","Exoc3","Rgs7bp","Elovl7","Rab3c","Top2b","Bmpr1a","Gm20642","Fermt2","Exoc5","Rem2","Cmtm5","Xkr6","Pcdh17","Pcdh20","Sema5a","Gm31462","Desi1","Shisal1","Tuba1a","Shisa9","2510009E07Rik","Nectin3","Ift57","Cadm2","Slc5a3","Pde10a","Syngr3","Srpk1","Mdga1","Zfp871","Eml4","Prkce","Nrxn1","Rnf138","Dtna","Gpr17","Kctd16","Atg12","Ap3s1","Sema6a","Csnk1g3","Dok6","Tmx3","Rfx3","Apln","Igsf1","Fmr1","Atp2b3","Plp1","Trpc5","Cnksr2","Glra2","Kdm5d","Uty")


rld_brain <- rlogTransformation(br_dds)

rld_mat_brain <- assay(rld_brain)
rld_mat_brain

cluster_brain_rlog <- rld_mat_brain[genelist_brain2, ]

clusters <- degPatterns(cluster_brain_rlog, metadata = colData(br_dds), time = "condition", minc = 0, plot = T, cutoff = 0.7)

library(ggplot2)
library(viridis)

clusters <- degPatterns(cluster_liver_rlog, metadata = colData(liver_dds), time = "condition", minc = 0, plot = T, cutoff = 0.7) +
  scale_color_viridis_d(option = "viridis") + # "viridis", "magma", "plasma" 중 선택
  theme_minimal() +
  theme(legend.position = "bottom") # 범례 위치 조정

clusters <- degPatterns(cluster_brain_rlog, metadata = colData(br_dds), time = "condition", minc = 0, plot = T, cutoff = 0.7)
filtclusters<-clusters$normalized
filtclusters$cluster
upc = c("1")
downc = c("2")

filtclusters1<-filtclusters[filtclusters$cluster %in% upc,]
filtclusters1[filtclusters1$cluster %in% upc,"cluster"] = 1
p1 <- degPlotCluster(filtclusters1,time="condition",color="condition")+geom_line(aes_string(group="genes"),alpha=0.05)+scale_color_nejm()+theme_bw()+ theme(legend.position = "none")+ theme(axis.title.y = element_blank())

filtclusters2<-filtclusters[filtclusters$cluster %in% downc,]
filtclusters2[filtclusters2$cluster %in% downc,"cluster"] = 2
p2 <- degPlotCluster(filtclusters2,time="condition",color="condition")+geom_line(aes_string(group="genes"),alpha=0.05)+scale_color_nejm()+theme_bw()+ theme(legend.position = "none")+ theme(axis.title.y = element_blank())

p <- plot_grid(p1, p2)
p

setwd(file.path(root, "39_fibroblasts/2_RNA_seq_tissue/final/"))

ggsave("DEGreport_brain_1216.pdf", plot = p, height = 2, width = 12, units = "in", dpi = 300)



############################################################################################################################

setwd(file.path(root, "41_RNA_seq_test/"))

kd_mat <- as.data.frame(read.csv(file.path(root, "39_fibroblasts/2_RNA_seq_tissue/public/kidney/bam/bam/counts_kidney_time.csv")))
row.names(kd_mat) <- kd_mat[,"Geneid"]
kd_mat <- as.data.frame(kd_mat[,-1])

kd_col1 <- colnames(kd_mat)
kd_col1

condition = c(rep('bbb',8),  rep('aaa',6), rep('ccc',7), rep('aaa',5), rep('bbb',4))
condition


exp.data <- data.frame(
  row.names=colnames(kd_mat),
  condition = condition)

exp.data

kd_dds <- DESeqDataSetFromMatrix(countData = kd_mat,
                                 colData = exp.data,
                                 design = ~ condition)

kd_dds <- DESeq(kd_dds)


genelist_kidney <- c("Lactb2", "Cnot11", "Osgepl1", "1700066M21Rik", "Spats2l", "Tmem169", "Gpc1", "Slc35f5", "Cryzl2", "Pigc", "Sccpdh", "Enah", "Capn2", "Disp1", "Nek2", "Stam", "Col5a1", "Wdsub1", "Spc25", "Ift70a1", "Cat", "E2f1", "Rprd1b", "Chrna4", "Car3", "Ccna2", "Abhd18", "Etfdh", "Tmem144", "Sh3d19", "Cgn", "Txnip", "Hsd3b2", "Slc16a1", "Capza1", "Papss1", "Npnt", "Ifi44", "Decr1", "Nbn", "Rad23b", "Cyp2j5", "C8a", "Cdc20", "Crot", "Bst1", "Lap3", "Adgra3", "Corin", "Dcun1d4", "Afm", "Lrrc8d", "Rpap2", "Fam20c", "Gna12", "Tfpi2", "Umad1", "Irf5", "Smo", "Gmcl1", "Rbsn", "Raf1", "Ift122", "Bid", "Foxm1", "Slco1a1", "Opa3", "Cebpa", "Klk1b4", "Tmem86a", "Gab2", "Plk1", "Ap2a2", "Cdc16", "Tdrp", "Mcph1", "Erlin2", "Acsl1", "Casp3", "Gpm6a", "Fbxo8", "Csgalnact1", "Phkb", "Lonp2", "Mmp15", "Cfap20", "Ces2b", "Ces2c", "Nudt7", "Zfp426", "Ncam1", "Acat1", "Dis3l", "Fam81a", "Elovl5", "Sh3bgrl2", "Trpc1", "Nmnat3", "Dag1", "Scap", "Pth1r", "Dync1li1", "Lama4", "Cd164", "Slc29a3", "Herc4", "Cdk1", "Sppl2b", "Nopchap1", "Inpp5j", "Gal3st1", "Emid1", "Gm11992", "Hmmr", "Prpsap2", "Ulk2", "Slc16a13", "Tefm", "Arsg", "Zfp750", "Lratd1", "Tmem18", "Hbp1", "Klhdc2", "Vipas39", "Tnfaip2", "Akr1c14", "Gpr137b", "AW209491", "Elovl2", "Ubqln1", "Naa35", "Zfp874a", "Arrdc3", "Ccnb1", "Ngly1", "4931414P19Rik", "Pbk", "Ghr", "Ugt3a1", "Enpp2", "Mapk15", "Racgap1", "Cers5", "Snap29", "Ehhadh", "Bdh1", "Pcyt1a", "Btnl12", "Pde9a", "Zfp101", "Cdc5l", "Ltbp1", "Ammecr1l", "Nrep", "Dctn4", "Zfp532", "Spire1", "Gm33703", "Slc22a28", "Fads3", "Keg1", "C330002G04Rik", "Ermp1", "A1cf", "Lipo3", "Loxl4", "Acsl5", "Hs6st2", "Slc10a3", "Slc16a2", "Bhlhb9", "Cldn2", "Col4a5", "Nxt2", "Fgd1", "Prps2")

genelist_kidney2 <- c("Gm4208","Pecr","Mogat1","Steap3","Fmod","Rnpep","Gm8947","Dpt","Kcnh1","Dipk1b","Pomt1","Aplnr","Rhov","Ckmt1","E2f1","Matn4","Kcng1","Gm6300","Car3","Gm33051","Insrr","Rhbg","Cct3","Ssr2","Dcst1","Smcp","H2bc21","H3c15","Phgdh","Sars1","Gm40137","Adh1","Ubap2","Stra6l","Col15a1","Mast2","Cdc20","Zmynd12","Cdca8","Ldc1","Arhgef19","Fam131c","Ubiad1","Perm1","Eif2b4","Htra3","Evc","Cckar","Gm42114","Gm10399","Castor2","Cldn3","Gper1","Smo","Tmem150a","D830050J10Rik","Efcab12","Ift122","Mfap5","Foxm1","Fkbp4","Slco1a1","Wdr62","Cebpa","Gm35724","Myh14","Tmem86a","Aldh1a3","Wdr93","Mogat2","Serpinh1","Ucp2","Clpb","Acsm3","Htra1","Chid1","Cars2","Gm15348","Gm39133","Pllp","Mmp15","Ces2b","Hsd11b2","Plcg2","Hsd17b2","Gm31774","Cldn18","Col6a6","Slc26a6","Vill","Clec3b","Mfsd4b3-ps","Cdk1","Ggt1","Gm867","Col6a2","Col18a1","Igf1","Rdh5","Nhp2","Serpinf2","Car4","Dcakd","Wnt9b","Mapt","Ace","Fasn","Dpf3","Larp4b","H4c18","Wnk2","Clptm1l","Ccnb1","D14Ertd670e","Nt5dc2","Eddm3b","4931414P19Rik","Nynrin","Atp12a","Trim35","AA536875","Osr2","Kdelr3","Rtl6","Nckap5l","Aifm3","Vps8","Fbxo40","Nhlrc4","Cryaa","Mllt1","Slc30a6","Smim3","Cndp2","Acy3","Tbx10","Cabp2","Slc22a12","Gm32245","Slc22a6","Fads3","Aldh18a1","Chst7","Apln","Fgd1","Asb9")

rld_kidney <- rlogTransformation(kd_dds)

rld_mat_kidney <- assay(rld_kidney)
rld_mat_kidney

cluster_kidney_rlog <- rld_mat_kidney[genelist_kidney2, ]

clusters <- degPatterns(cluster_kidney_rlog, metadata = colData(kd_dds), time = "condition", minc = 0, plot = T, cutoff = 0.7)


library(ggsci)
library(cowplot)
library(DEGreport)

clusters <- degPatterns(cluster_kidney_rlog, metadata = colData(kd_dds), time = "condition", minc = 0, plot = T, cutoff = 0.7)
filtclusters<-clusters$normalized
filtclusters$cluster
upc = c("1")
downc = c("2")

filtclusters1<-filtclusters[filtclusters$cluster %in% upc,]
filtclusters1[filtclusters1$cluster %in% upc,"cluster"] = 1
p1 <- degPlotCluster(filtclusters1,time="condition",color="condition")+geom_line(aes_string(group="genes"),alpha=0.05)+scale_color_nejm()+theme_bw()+ theme(legend.position = "none")+ theme(axis.title.y = element_blank())
p1

filtclusters2<-filtclusters[filtclusters$cluster %in% downc,]
filtclusters2[filtclusters2$cluster %in% downc,"cluster"] = 2
p2 <- degPlotCluster(filtclusters2,time="condition",color="condition")+geom_line(aes_string(group="genes"),alpha=0.05)+scale_color_nejm()+theme_bw()+ theme(legend.position = "none")+ theme(axis.title.y = element_blank())
p2

p <- plot_grid(p1, p2)
combined_plot <- p1 / p2
plot_grid(p1, p2, ncol = 2, align = "v")

setwd(file.path(root, "39_fibroblasts/2_RNA_seq_tissue/final/"))


ggsave("DEGreport_kidney_1204.pdf", plot = p2, height = 2, width = 6, units = "in", dpi = 300)


#######################################################################################################

egm_mat <- as.data.frame(read.csv(file.path(root, "39_fibroblasts/2_RNA_seq_tissue/only_egm_effect_in_tissue/counts_egm.csv")))
row.names(egm_mat) <- egm_mat[,"Geneid"]
egm_mat <- as.data.frame(egm_mat[,-1])

egm_col1 <- colnames(egm_mat)
egm_col1

condition = c(rep('progeria_control',3),  rep('progeria_egm',3), rep('progeria_control',3), rep('progeria_egm',3), rep('progeria_control',3),  rep('progeria_egm',3))
condition


exp.data <- data.frame(
  row.names=colnames(egm_mat),
  condition = condition)

exp.data

egm_dds <- DESeqDataSetFromMatrix(countData = egm_mat,
                                    colData = exp.data,
                                    design = ~ condition)

egm_dds <- DESeq(egm_dds)

# 1. DESeq2 결과로 ranks 생성

library(msigdbr)
library(fgsea)

res <- results(egm_dds, contrast = c("condition", "progeria_egm", "progeria_control"))


write.csv(res, file.path(root, "39_fibroblasts/2_RNA_seq_tissue/only_egm_effect_in_tissue/egm_mechanism.csv"))

ranks <- res$log2FoldChange
names(ranks) <- rownames(res)
ranks <- ranks[!is.na(ranks)]
ranks <- sort(ranks, decreasing = TRUE)

filtered_res <- res[res$pvalue < 0.05 & !is.na(res$pvalue), ]

# log2FoldChange를 ranks로 설정
ranks <- filtered_res$log2FoldChange

# gene 이름으로 ranks의 이름 설정
names(ranks) <- rownames(filtered_res)

# NA 제거 (필터링 과정에서 이미 제거되었겠지만, 안정성을 위해 추가)
ranks <- ranks[!is.na(ranks)]

# ranks를 내림차순으로 정렬
ranks <- sort(ranks, decreasing = TRUE)

# ranks 출력 확인
head(ranks)

mouse_gobp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")

filtered_gobp <- mouse_gobp %>%
  filter(grepl("CELL|CYCLE|TRANSPORT|REGENERATION", gs_name, ignore.case = TRUE))

filtered_gobp

mouse_gobp_pathways <- split(x = filtered_gobp$gene_symbol, f = filtered_gobp$gs_name)

# 결과 확인
length(mouse_gobp_pathways)  # 필터링된 경로 수

# fgsea 실행
fgsea_res <- fgsea(
  pathways = mouse_gobp,
  stats = ranks,
  nperm = 1000
)

significant_pathways <- fgsea_res[fgsea_res$pval < 0.05, ]

significant_pathways




selected_pathways <- c("GOBP_CELLULAR_HOMEOSTASIS", "GOBP_ION_TRANSMEMBRANE_TRANSPORT", "GOBP_CALCIUM_ION_TRANSPORT", "GOBP_CELL_CELL_ADHESION")

# 필터링
filtered_fgsea_res <- fgsea_res %>%
  filter(pathway %in% selected_pathways) %>%
  mutate(pathway = factor(pathway, pathway[order(NES, decreasing = TRUE)]))


ggplot(filtered_fgsea_res, aes(x = pathway, y = NES, fill = -log10(pval))) +
  geom_bar(stat = "identity", width = 0.8, color = "black") +
  coord_flip() +
  labs(
    title = "Top Enriched Pathways",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)",
    fill = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_gradient(
    low = "pink",    # 연한 빨간색
    high = "red"     # 진한 빨간색
  ) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

ggplot(filtered_fgsea_res, aes(x = pathway, y = NES, fill = -log10(pval))) +
  geom_bar(stat = "identity", width = 0.8, color = "black") +  # 바 테두리
  labs(
    title = "Top Enriched Pathways",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)",
    fill = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_gradient(
    low = "pink",    # 연한 빨간색
    high = "red"     # 진한 빨간색
  ) +
  scale_y_continuous(
    breaks = seq(0, 1.5, by = 0.5),  # y축 눈금 표시 간격 설정
    limits = c(0, 1.5)               # y축 범위 설정
  ) +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),        # y축 레이블 색상
    axis.ticks.y = element_line(color = "black", size = 0.5),     # y축 눈금선 색상
    axis.line.y = element_line(color = "black", size = 0.8),      # y축 축선 색상
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1), # x축 레이블 스타일
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # 그래프 바깥쪽 테두리
  )
