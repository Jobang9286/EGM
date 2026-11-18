#!/usr/bin/env python
"""Collect regions, gene annotation and GO enrichment into tidy tables for R.

Usage: 09_export.py OUTDIR MARK:HEATMAP_DIR [MARK:HEATMAP_DIR ...]

Everything comes out long-format with mark and class as columns, so the whole
set reads in as one data frame per file and facets without reshaping.
"""
import os
import sys
from collections import defaultdict
from importlib.machinery import SourceFileLoader

BIN = os.path.dirname(os.path.abspath(__file__))
_e = SourceFileLoader("e", os.path.join(BIN, "07_go_enrich.py")).load_module()

CLASSES = [("4FP_specific", "4FP_specific"),
           ("shared", "Shared"),
           ("CP_specific", "CP_specific")]


def main():
    out = sys.argv[1]
    specs = [a.split(":", 1) for a in sys.argv[2:]]
    os.makedirs(out, exist_ok=True)
    os.makedirs(os.path.join(out, "by_mark"), exist_ok=True)

    name, ns, parents = _e.parse_obo(_e.OBO)
    sys.setrecursionlimit(10000)
    gene2go = _e.load_annotations(_e.GAF, "P", parents)
    go_ns = ns

    reg = open(os.path.join(out, "regions_annotated.tsv"), "w")
    gen = open(os.path.join(out, "genes.tsv"), "w")
    go = open(os.path.join(out, "go_bp.tsv"), "w")
    smry = open(os.path.join(out, "summary.tsv"), "w")
    reg.write("mark\tclass\tchrom\tstart\tend\twidth\tgene\tgene_type\tstrand\tdist_to_tss\tfeature\n")
    gen.write("mark\tclass\tgene\tgene_type\tn_regions\tnearest_dist\n")
    go.write("mark\tclass\tgo_id\tterm\tk\tn\texpected\tfold\tp\tfdr\n")
    smry.write("mark\tclass\tn_regions\tn_genes\tn_genes_pc\n")

    for mark, hm in specs:
        classes_pc = {}
        for fn, label in CLASSES:
            src = os.path.join(hm, f"{fn}.annotated.tsv")
            if not os.path.exists(src):
                continue
            per_gene = defaultdict(lambda: [0, None, None])
            nreg = 0
            with open(src) as fh:
                next(fh)
                lines = []
                for line in fh:
                    c, s, e, g, gt, st, d, cls = line.rstrip("\n").split("\t")
                    nreg += 1
                    row = f"{mark}\t{label}\t{c}\t{s}\t{e}\t{int(e)-int(s)}\t{g}\t{gt}\t{st}\t{d}\t{cls}\n"
                    reg.write(row); lines.append(row)
                    if g == ".":            # contig with no gene at all
                        continue
                    rec = per_gene[(g, gt)]
                    rec[0] += 1
                    if rec[1] is None or abs(int(d)) < abs(rec[1]):
                        rec[1] = int(d)
            with open(os.path.join(out, "by_mark", f"{mark}_{fn}.annotated.tsv"), "w") as fh:
                fh.write("mark\tclass\tchrom\tstart\tend\twidth\tgene\tgene_type\tstrand\tdist_to_tss\tfeature\n")
                fh.writelines(lines)
            for (g, gt), (n, d, _) in sorted(per_gene.items()):
                if g != ".":
                    gen.write(f"{mark}\t{label}\t{g}\t{gt}\t{n}\t{d}\n")
            pc = _e.read(os.path.join(hm, f"{fn}.genes_pc.txt"))
            classes_pc[label] = pc
            allg = _e.read(os.path.join(hm, f"{fn}.genes.txt"))
            smry.write(f"{mark}\t{label}\t{nreg}\t{len(allg)}\t{len(pc)}\n")

        bg = set().union(*classes_pc.values()) if classes_pc else set()
        for label, st in classes_pc.items():
            if not st or not bg:
                continue
            res, tested = _e.enrich(st, bg, gene2go, name, 10 ** 9)
            for q, p, t, k, n, f in res:
                if go_ns.get(t) != "P":
                    continue
                exp = n * len(st) / len(bg)
                nm = name.get(t, t).replace("\t", " ")
                go.write(f"{mark}\t{label}\t{t}\t{nm}\t{k}\t{n}\t{exp:.2f}\t{f:.4f}\t{p:.6g}\t{q:.6g}\n")

    for fh in (reg, gen, go, smry):
        fh.close()
    for f in ("regions_annotated.tsv", "genes.tsv", "go_bp.tsv", "summary.tsv"):
        p = os.path.join(out, f)
        print(f"{f:<26} {sum(1 for _ in open(p)) - 1:>8} rows")


if __name__ == "__main__":
    main()
