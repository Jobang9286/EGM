import pandas as pd
cat={"Stat1":"Immune","Stat2":"Immune","Irf7":"Immune","Rela":"Immune",
     "Foxm1":"Cell cycle","E2f6":"Cell cycle","Tead1":"Hippo","Tead2":"Hippo","Twist2":"Developmental",
     "Ctcf":"Chromatin","Mbd2":"Chromatin","Smarcc2":"Chromatin","Setdb1":"Chromatin","Trim28":"Chromatin","Hdac2":"Chromatin","Smc3":"Chromatin",
     "Ets1":"ETS","Elf1":"ETS","Ets2":"ETS","Elk3":"ETS","Etv5":"ETS","Gabpa":"ETS","Fli1":"ETS"}
d=pd.concat([pd.read_csv(f"phase3_figs/ps_{t}/ps_scores.csv") for t in ["negctrl","negctrl2","negctrl3"]]).drop_duplicates("TF")
d["category"]=d.TF.map(cat)
d["no_reversal"]=d.p_neg_early>0.05
d=d.sort_values(["category","mean_ps_early"])
d.to_csv("phase3_figs/PS_negctrl_all.csv",index=False)
print(d[["category","TF","mean_ps_early","p_neg_early","no_reversal","cell_cos_mean"]].round(4).to_string(index=False))
