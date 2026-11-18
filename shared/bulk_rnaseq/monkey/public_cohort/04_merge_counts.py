import pandas as pd, os
pub="public_monkey/proj"; egm="monkey/project/quant"
dl=pd.read_csv("public_monkey/download_list.csv")
out="public_monkey/merged"; os.makedirs(out,exist_ok=True)
def col(f): return pd.read_csv(f,sep='\t',usecols=['gene_id','expected_count']).set_index('gene_id')['expected_count'].round().astype(int)
for organ in ["Brain","Liver","Kidney","Spleen"]:
    cols={}; meta=[]
    for _,r in dl[dl.organ==organ].iterrows():
        f=f"{pub}/{organ}/{r.CRR}.genes.results"
        if not os.path.exists(f): print("MISS",f); continue
        n=f"P_{r.age}_{r.CRR}"; cols[n]=col(f); meta.append((n,r.age,"public"))
    for grp,lab in [("EGM","EGM"),("control","CTRL")]:
        for i in (1,2,3):
            f=f"{egm}/{grp}_{organ}_{i}.genes.results"; n=f"{lab}_{i}"
            cols[n]=col(f); meta.append((n,lab,"egm"))
    mat=pd.DataFrame(cols); mat.index.name="gene_id"
    mat.to_csv(f"{out}/{organ}_counts.csv")
    pd.DataFrame(meta,columns=["sample","group","source"]).to_csv(f"{out}/{organ}_meta.csv",index=False)
    print(organ, "genes x samples =", mat.shape, "| groups:", pd.Series([m[1] for m in meta]).value_counts().to_dict())
