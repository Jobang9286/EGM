import pandas as pd
c = pd.read_csv("phase2_centrality.csv")
piv = c[c.degree_centrality_out>0].pivot_table(index="gene",columns="cluster",values="degree_centrality_out").fillna(0).astype(int)
go = pd.read_csv("lipid_go_genes.csv").set_index("SYMBOL")
already = set(l.strip() for l in open("lipid_tf_list.txt")) | {"Srebf2"}
hits = [g for g in piv.index if g in go.index]
new = [g for g in hits if g not in already]
tab = piv.loc[new].copy(); tab["GO_category"] = go.loc[new,"category"]
tab["max_outdeg"] = tab[["Q1","Q2","Q3","Q4","Q5"]].max(axis=1)
tab = tab.sort_values("max_outdeg", ascending=False)
print(f"lipid-GO regulators total={len(hits)}, already running={len(hits)-len(new)}, NEW={len(new)}")
print(tab.to_string())
open("lipid2_tf_list.txt","w").write("\n".join(tab.index)+"\n")
