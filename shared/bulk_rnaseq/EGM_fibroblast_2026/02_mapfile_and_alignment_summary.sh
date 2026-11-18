#!/bin/bash
ROOT=$PWD

## ---- count matrices with EGM labels, mapfile -----------------------------------
cd ${ROOT}/egm_blood
# rebuild the matrices with EGM labels
python3 - <<'PY'
import pandas as pd, os
PROJ="."
S=["GFP-1","GFP-2","GFP-3","GME-1","GME-2","GME-3","Eno1-1","Gclc-1","Map1s-1"]
R={"GME-1":"EGM-1","GME-2":"EGM-2","GME-3":"EGM-3"}
cnt,tpm={},{}
for s in S:
    d=pd.read_csv(f"{PROJ}/{s}.genes.results",sep="\t",index_col=0)
    cnt[s],tpm[s]=d["expected_count"],d["TPM"]
C=pd.DataFrame(cnt).rename(columns=R); T=pd.DataFrame(tpm).rename(columns=R)
C.to_csv(f"{PROJ}/txi.rsem.counts.txt",sep="\t"); T.to_csv(f"{PROJ}/txi.rsem.abundance.TPM.txt",sep="\t")
C.round().astype(int).to_csv(f"{PROJ}/counts_symbol.csv")
print("columns:", ", ".join(C.columns))
print(f"genes={C.shape[0]}")
PY
# mapfile with EGM condition, keeping the fastq prefix for traceability
{ printf "sample\tTREATMENT\tfastq_prefix\n"
  for i in 1 2 3; do printf "GFP-%d\tuntreated\tGFP-%d\n" $i $i; done
  for i in 1 2 3; do printf "EGM-%d\tegm\tGME-%d\n" $i $i; done
  printf "Eno1-1\teno1\tEno1-1\nGclc-1\tgclc\tGclc-1\nMap1s-1\tmap1s\tMap1s-1\n"; } > mapfile.tsv
cat mapfile.tsv

## ---- alignment summary ---------------------------------------------------------
cd ${ROOT}/egm_blood
cat > /tmp/summ.py <<'PY'
import pandas as pd, re, os
PROJ = "."
S = ["GFP-1","GFP-2","GFP-3","GME-1","GME-2","GME-3","Eno1-1","Gclc-1","Map1s-1"]
R = {"GME-1":"EGM-1","GME-2":"EGM-2","GME-3":"EGM-3"}
rows = []
for s in S:
    L = dict(re.findall(r"\s*(.+?)\s*\|\s*(.+)", open(f"{PROJ}/{s}Log.final.out").read()))
    g = lambda k: L.get(k, "").strip()
    rpt = open(f"{PROJ}/trim_galore/{s}_1.fastq.gz_trimming_report.txt").read()
    raw = re.search(r"Total reads processed:\s+([\d,]+)", rpt).group(1).replace(",", "")
    d = pd.read_csv(f"{PROJ}/{s}.genes.results", sep="\t", index_col=0)
    rows.append(dict(sample=R.get(s, s),
                     raw_M=round(int(raw)/1e6, 1),
                     star_in_M=round(int(g("Number of input reads"))/1e6, 1),
                     unique=g("Uniquely mapped reads %"),
                     multi=g("% of reads mapped to multiple loci"),
                     too_short=g("% of reads unmapped: too short"),
                     counts_M=round(d["expected_count"].sum()/1e6, 1),
                     genes_det=int((d["expected_count"] >= 1).sum())))
A = pd.DataFrame(rows)
A.to_csv(f"{PROJ}/alignment_summary.csv", index=False)
print(A.to_string(index=False))
PY
python3 /tmp/summ.py
