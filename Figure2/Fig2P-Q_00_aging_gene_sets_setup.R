## ---- packages, project root -------------------------------------------------
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ComplexHeatmap); library(circlize); library(ggplot2); library(patchwork)})
root    <- "EGM_project_for_deconvolution"
# ---- 2. Aging gene sets (mouse symbols) ----
aging_sets <- list(
  Senescence = c("Cdkn1a","Cdkn2a","Cdkn2b","Cdkn1b","Cdkn1c","Trp53","Glb1","Serpine1","Lmnb1","Hmgb1","Mdm2","Rb1"),
  SASP = c("Il6","Il1a","Il1b","Il18","Il1rn","Cxcl1","Cxcl2","Cxcl10","Cxcl12","Cxcl16","Ccl2","Ccl3","Ccl4","Ccl5","Ccl7","Ccl8","Ccl20",
           "Mmp3","Mmp9","Mmp12","Mmp13","Mmp14","Timp1","Igfbp2","Igfbp3","Igfbp4","Igfbp5","Igfbp7","Tnf","Vegfa","Areg","Ereg",
           "Csf1","Csf2","Fgf2","Hgf","Plaur","Lif","Osm","Il15","Icam1","Vcam1","Gdf15"),
  Inflammation = c("Nfkb1","Nfkb2","Rela","Relb","Ikbkb","Tlr2","Tlr4","Myd88","Stat1","Stat3","Irf1","Irf7","C3","C1qa","C1qb","C1qc","C4b",
                   "B2m","H2-K1","H2-D1","Cd74","Lgals3","Cd68","Adgre1","Ptprc","S100a8","S100a9","Nlrp3","Casp1","Cd44","Ifit1","Ifit3",
                   "Isg15","Oasl1","Cxcl9","Ccr2","Itgam","Trem2"),
  Fibrosis = c("Col1a1","Col1a2","Col3a1","Col4a1","Col4a2","Col5a1","Col5a2","Col6a1","Col6a3","Fn1","Acta2","Tagln","Tgfb1","Tgfb2","Tgfb3",
               "Ccn2","Ctgf","Postn","Vim","Lox","Loxl2","Sparc","Mmp2","Thbs1","Thbs2","Eln","Fbn1","Dcn","Lum","Pdgfrb","Pdgfb","Des","Fap",
               "Timp2","Serpinh1","Bgn","Tnc"),
  Injury_Damage = c("Havcr1","Lcn2","Spp1","Klf6","Atf3","Sox9","Krt20","Clu","Tnfrsf12a","Cst3","Gpnmb","Ccn1","Cyr61","Egr1","Fos","Jun","Junb",
                    "Lgals1","Anxa2","Anxa3","Cp","Tnfaip3","Nfkbia","Ier3"),
  Longevity_Nutrient = c("Kl","Klb","Sirt1","Sirt6","Sirt7","Foxo1","Foxo3","Foxo4","Mtor","Rptor","Rictor","Prkaa1","Prkaa2","Igf1","Igf1r","Insr",
                         "Irs1","Irs2","Fgf23","Fgf21","Nampt","Nmnat1","Cd38","Akt1","Tsc1","Tsc2","Ppard","Adipor1"),
  Tubular_Function = c("Umod","Egf","Lrp2","Cubn","Slc34a1","Slc5a2","Slc22a6","Slc22a8","Aqp1","Aqp2","Slc12a1","Slc12a3","Nphs1","Nphs2","Podxl",
                       "Wt1","Cdh16","Hnf1b","Hnf4a","Pax2","Pax8","Ggt1","Kap","Calb1","Pvalb","Fxyd2","Atp1a1"),
  DNA_Damage = c("Atm","Atr","Chek1","Chek2","H2ax","H2afx","Gadd45a","Gadd45b","Gadd45g","Parp1","Brca1","Brca2","Rad51","Xrcc1","Xrcc5","Xrcc6",
                 "Trp53bp1","Mre11a","Nbn","Rad50","Ddb2","Ercc1","Ercc2","Ogg1","Apex1","Lig4","Prkdc","Mgmt","Rpa1","Rad23b","Xpa","Xpc"),
  Oxidative_Stress = c("Sod1","Sod2","Sod3","Cat","Gpx1","Gpx3","Gpx4","Nfe2l2","Keap1","Hmox1","Nqo1","Txn1","Txnrd1","Prdx1","Prdx3","Prdx6",
                       "Gclc","Gclm","Gsr","Gstm1","Gsta4","Gstp1","Nox4","Cybb","Ncf1","Ncf2","Srxn1","Mt1","Mt2"),
  Mitochondrial = c("Ppargc1a","Ppargc1b","Tfam","Nrf1","Sirt3","Sirt4","Sirt5","Mfn1","Mfn2","Opa1","Dnm1l","Fis1","Cox4i1","Cox5a","Ndufa4","Ndufs1",
                    "Ndufb8","Sdha","Sdhb","Uqcrc1","Uqcrc2","Atp5a1","Atp5f1a","Cycs","Pink1","Prkn","Ucp2","Cpt1a","Acadm","Ppara","Esrra","Polg","Tomm20","Vdac1"),
  Autophagy_Proteostasis = c("Atg5","Atg7","Atg12","Atg16l1","Becn1","Map1lc3b","Map1lc3a","Sqstm1","Ulk1","Ulk2","Lamp1","Lamp2","Ctsb","Ctsd","Ctsl","Tfeb",
                             "Hspa1a","Hspa1b","Hspa5","Hspa8","Hsp90aa1","Hspb1","Hspd1","Dnajb1","Hsf1","Bnip3","Psmb5","Psmd1","Psma1","Ubb","Ubc",
                             "Atf4","Ddit3","Xbp1","Ern1","Eif2ak3","Nbr1","Wipi2","Gabarap"),
  Epigenetic = c("Dnmt1","Dnmt3a","Dnmt3b","Tet1","Tet2","Tet3","Hdac1","Hdac2","Hdac3","Hdac4","Kdm6a","Kdm6b","Kdm5b","Kdm4a","Ezh2","Kat2a","Kat6a",
                 "Setd7","Suv39h1","Cbx5","H3f3b","Brd4","Kmt2a","Ep300","Crebbp","Mecp2"),
  Nuclear_Lamina = c("Lmna","Lmnb2","Emd","Lbr","Sun1","Sun2","Zmpste24","Banf1","Tmpo","Syne1","Syne2","Nup153")
)
# 방향 주석: Pro-aging(노화 시 증가/유해) vs Protective(항노화/기능 유지) vs Ambiguous
protective <- c("Lmnb1", aging_sets$Tubular_Function, aging_sets$Mitochondrial,
  "Sod1","Sod2","Sod3","Cat","Gpx1","Gpx3","Gpx4","Nfe2l2","Hmox1","Nqo1","Txn1","Txnrd1","Prdx1","Prdx3","Prdx6","Gclc","Gclm","Gsr","Gstm1","Gsta4","Gstp1","Srxn1","Mt1","Mt2",
  "Kl","Klb","Sirt1","Sirt6","Sirt7","Foxo1","Foxo3","Foxo4","Prkaa1","Prkaa2","Nampt","Nmnat1","Tsc1","Tsc2","Ppard","Adipor1","Fgf21",
  "Atg5","Atg7","Atg12","Atg16l1","Becn1","Map1lc3b","Map1lc3a","Ulk1","Ulk2","Lamp1","Lamp2","Ctsb","Ctsd","Ctsl","Tfeb","Hspa1a","Hspa1b","Hspa5","Hspa8",
  "Hsp90aa1","Hspb1","Hspd1","Dnajb1","Hsf1","Psmb5","Psmd1","Psma1","Ubb","Ubc","Nbr1","Wipi2","Gabarap",
  "Brca1","Brca2","Rad51","Xrcc1","Xrcc5","Xrcc6","Trp53bp1","Mre11a","Nbn","Rad50","Ddb2","Ercc1","Ercc2","Ogg1","Apex1","Lig4","Prkdc","Mgmt","Rpa1","Rad23b","Xpa","Xpc",
  "Lmnb2","Emd","Lbr","Sun1","Sun2","Zmpste24","Banf1","Tmpo","Syne1","Syne2","Nup153")
ambiguous <- c(aging_sets$Epigenetic, "Atm","Atr","Chek1","Chek2","H2ax","H2afx","Gadd45a","Gadd45b","Gadd45g","Parp1","Lmna","Hmgb1","Trp53","Rb1","Mdm2","Cd38")

## grid.text()/gpar() in the heatmap cell functions need it.
library(grid)

