#!/usr/bin/env python
"""Join the region annotation to the RNA-seq differential expression table.

Usage: 10_join_de.py OUTDIR DE_CSV MARK:HEATMAP_DIR [MARK:HEATMAP_DIR ...]

Both sides are GENCODE vM38 symbol-level, so gene symbols join directly.
logFC is 4FP over CP: positive means higher in 4FP.

Two tables:

  genes_de_long.tsv   one row per gene per mark per class, for grouped
                      comparisons. Genes carrying no region of a mark are
                      emitted with class "none" - without that baseline there
                      is nothing to compare a class against.
  gene_matrix.tsv     one row per gene, region counts spread across columns,
                      for scatter and correlation work.
"""
import csv
import os
import sys
from collections import defaultdict

CLASSES = ["4FP_specific", "Shared", "CP_specific"]
FILES = {"4FP_specific": "4FP_specific", "Shared": "shared", "CP_specific": "CP_specific"}


def main():
    out, de_csv = sys.argv[1], sys.argv[2]
    specs = [a.split(":", 1) for a in sys.argv[3:]]
    os.makedirs(out, exist_ok=True)

    de = {}
    with open(de_csv) as fh:
        for r in csv.DictReader(fh):
            de[r["gene"]] = r
    print(f"DE table: {len(de)} genes")

    # mark -> class -> gene -> [n_regions, n_promoter, nearest_abs_dist]
    data = {m: {c: defaultdict(lambda: [0, 0, None]) for c in CLASSES} for m, _ in specs}
    for mark, hm in specs:
        for cls in CLASSES:
            path = os.path.join(hm, f"{FILES[cls]}.annotated.tsv")
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                next(fh)
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    g, d, feat = f[3], f[6], f[7]
                    if g == ".":
                        continue
                    rec = data[mark][cls][g]
                    rec[0] += 1
                    if feat == "promoter":
                        rec[1] += 1
                    if d != "NA":
                        ad = abs(int(d))
                        if rec[2] is None or ad < rec[2]:
                            rec[2] = ad

    DECOLS = ["logFC", "AveExpr", "t", "P.Value", "adj.P.Val"]
    long_path = os.path.join(out, "genes_de_long.tsv")
    with open(long_path, "w") as fh:
        fh.write("gene\tmark\tclass\tn_regions\tn_promoter\tnearest_dist\t"
                 + "\t".join(DECOLS) + "\n")
        n = 0
        for mark, _ in specs:
            for g, r in de.items():
                dec = "\t".join(r[c] for c in DECOLS)
                hit = False
                for cls in CLASSES:
                    rec = data[mark][cls].get(g)
                    if rec:
                        hit = True
                        nd = "NA" if rec[2] is None else rec[2]
                        fh.write(f"{g}\t{mark}\t{cls}\t{rec[0]}\t{rec[1]}\t{nd}\t{dec}\n")
                        n += 1
                if not hit:
                    fh.write(f"{g}\t{mark}\tnone\t0\t0\tNA\t{dec}\n")
                    n += 1
    print(f"genes_de_long.tsv        {n} rows")

    marks = [m for m, _ in specs]
    cols = []
    for m in marks:
        for c in CLASSES:
            cols += [f"{m}_{c}", f"{m}_{c}_prom"]
    mat_path = os.path.join(out, "gene_matrix.tsv")
    with open(mat_path, "w") as fh:
        fh.write("gene\t" + "\t".join(DECOLS) + "\t" + "\t".join(cols) + "\n")
        for g, r in de.items():
            vals = []
            for m in marks:
                for c in CLASSES:
                    rec = data[m][c].get(g)
                    vals += [str(rec[0]) if rec else "0", str(rec[1]) if rec else "0"]
            fh.write(f"{g}\t" + "\t".join(r[c] for c in DECOLS) + "\t" + "\t".join(vals) + "\n")
    print(f"gene_matrix.tsv          {len(de)} rows, {len(cols)} region columns")


if __name__ == "__main__":
    main()
