tfs=["Srebf2"]+[l.strip() for l in open("lipid_tf_list.txt") if l.strip()]+[l.strip() for l in open("lipid2_tf_list.txt") if l.strip()]
tfs=list(dict.fromkeys(tfs)); print("n TFs:",len(tfs))
half=(len(tfs)+1)//2
open("v2a_tf_list.txt","w").write("\n".join(tfs[:half])+"\n"); open("v2b_tf_list.txt","w").write("\n".join(tfs[half:])+"\n")
print("v2a:",tfs[:half]); print("v2b:",tfs[half:])
