import pandas as pd, os
base="monkey/project/deg"
tissues=["Brain","Kidney","Liver","Spleen"]
print("%-8s %10s %6s %6s  (padj<0.05 & |log2FC|>=1)" % ("Tissue","sig0.05","up","down"))
for t in tissues:
    f=os.path.join(base,t,"DESeq2_results_annotated.csv")
    df=pd.read_csv(f)
    sig=df[(df.padj<0.05)&(df.log2FoldChange.abs()>=1)].copy()
    up=sig[sig.log2FoldChange>=1]; down=sig[sig.log2FoldChange<=-1]
    sig.sort_values("padj").to_csv(os.path.join(base,t,"DEG_padj0.05_lfc1.csv"),index=False)
    up.sort_values("log2FoldChange",ascending=False).to_csv(os.path.join(base,t,"DEG_up_padj0.05.csv"),index=False)
    down.sort_values("log2FoldChange").to_csv(os.path.join(base,t,"DEG_down_padj0.05.csv"),index=False)
    print("%-8s %10d %6d %6d" % (t,len(sig),len(up),len(down)))
print("\nSaved per tissue: DEG_padj0.05_lfc1.csv, DEG_up_padj0.05.csv, DEG_down_padj0.05.csv")
