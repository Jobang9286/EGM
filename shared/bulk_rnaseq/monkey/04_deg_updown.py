import pandas as pd, os
base="monkey/project/deg"
tissues=["Brain","Kidney","Liver","Spleen"]
hdr="%-8s | %6s %6s %6s | %6s %6s %6s" % ("Tissue","sig10","up10","dn10","sig05","up05","dn05")
print(hdr); print("-"*len(hdr))
rows=[]
for t in tissues:
    df=pd.read_csv(os.path.join(base,t,"DESeq2_results_annotated.csv"))
    def c(p):
        s=df[(df.padj<p)&(df.log2FoldChange.abs()>=1)]
        return len(s), int((s.log2FoldChange>=1).sum()), int((s.log2FoldChange<=-1).sum())
    s10,u10,d10=c(0.10); s05,u05,d05=c(0.05)
    print("%-8s | %6d %6d %6d | %6d %6d %6d" % (t,s10,u10,d10,s05,u05,d05))
    rows.append([t,s10,u10,d10,s05,u05,d05])
out=pd.DataFrame(rows,columns=["tissue","sig_padj0.10","up_0.10","down_0.10","sig_padj0.05","up_0.05","down_0.05"])
out.to_csv(os.path.join(base,"DEG_updown_summary.csv"),index=False)
print("\nSaved: deg/DEG_updown_summary.csv  (|log2FC|>=1)")
