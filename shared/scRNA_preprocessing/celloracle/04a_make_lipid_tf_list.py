import pandas as pd
c = pd.read_csv("phase2_centrality.csv")
reg = set(c[c.degree_centrality_out>0].gene)
groups = {
 "SREBP/cholesterol":   ["Srebf1","Nfya","Nfyb","Nfyc","Sp1"],
 "nuclear receptors":   ["Nr1h2","Nr1h3","Nr1h4","Ppara","Ppard","Pparg","Rxra","Rxrb","Rxrg","Nr2f2","Nr1d1","Nr1d2","Rora","Rorc","Thra","Thrb","Esrra","Nr5a2","Hnf4a","Nr0b2","Nr1i2","Nr1i3","Nr3c1"],
 "lipogenesis (ChREBP)":["Mlxipl","Mlx"],
 "adipogenic/CEBP":     ["Cebpa","Cebpb","Cebpd","Klf15","Klf5","Ebf1","Prdm16","Zfp423","Bhlhe40"],
 "ER stress/lipid":     ["Xbp1","Atf4","Atf6","Creb3l1","Creb3l2","Creb3l3","Ddit3"],
 "lysosome/lipophagy":  ["Tfeb","Tfe3","Mitf"],
 "oxidative/FOXO":      ["Nfe2l1","Nfe2l2","Foxo1","Foxo3"],
}
use, skip = [], []
for g, tfs in groups.items():
    ok = [t for t in tfs if t in reg]; no = [t for t in tfs if t not in reg]
    print(f"{g:22s} usable={ok}  not_regulator={no}")
    use += ok; skip += no
print("\nUSABLE:", len(use), " ".join(use))
open("lipid_tf_list.txt","w").write("\n".join(use)+"\n")
