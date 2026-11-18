#!/usr/bin/env python
"""GO enrichment for the region classes 04_heatmap.sh produces.

Usage: 07_go_enrich.py HEATMAP_DIR [HEATMAP_DIR ...] [--aspect P|F|C] [--top N]

Background is every protein-coding gene carrying a peak of the SAME assay, not
the genome. Peaks sit at promoters and CpG islands far more often than chance,
so a genome background returns the housekeeping-and-development terms that any
promoter-biased set returns, and says nothing about the contrast being tested.
The question here is: given a gene that carries this mark at all, does landing
in the 4FP-specific class predict a function?

Annotations are propagated up is_a and part_of before testing, because a gene
annotated to a leaf term is by definition also annotated to its ancestors, and
skipping that silently empties every general term.
"""
import os
import re
import sys
import gzip
from collections import defaultdict

from scipy.stats import hypergeom

REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../_setup/ref")
OBO = os.path.join(REF, "go-basic.obo")
GAF = os.path.join(REF, "mgi.gaf.gz")
MIN_TERM, MAX_TERM, MIN_HIT = 5, 1500, 3


def parse_obo(path):
    name, ns, parents = {}, {}, defaultdict(set)
    cur = None
    obsolete = set()
    for line in open(path):
        line = line.rstrip("\n")
        if line == "[Term]":
            cur = None
        elif line.startswith("id: GO:"):
            cur = line[4:]
        elif cur:
            if line.startswith("name: "):
                name[cur] = line[6:]
            elif line.startswith("namespace: "):
                ns[cur] = {"biological_process": "P", "molecular_function": "F",
                           "cellular_component": "C"}.get(line[11:], "?")
            elif line.startswith("is_a: GO:"):
                parents[cur].add(line[6:].split(" ")[0])
            elif line.startswith("relationship: part_of GO:"):
                parents[cur].add(line[22:].split(" ")[0])
            elif line == "is_obsolete: true":
                obsolete.add(cur)
    for t in obsolete:
        name.pop(t, None)
    return name, ns, parents


def ancestors(term, parents, cache):
    if term in cache:
        return cache[term]
    cache[term] = set()                      # guards the cyclic case
    out = set()
    for p in parents.get(term, ()):
        out.add(p)
        out |= ancestors(p, parents, cache)
    cache[term] = out
    return out


def load_annotations(gaf, aspect, parents):
    cache = {}
    gene2go = defaultdict(set)
    with gzip.open(gaf, "rt") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.split("\t")
            if len(f) < 10 or f[8] != aspect:
                continue
            if "NOT" in f[3]:
                continue
            sym, go = f[2], f[4]
            gene2go[sym].add(go)
            gene2go[sym] |= ancestors(go, parents, cache)
    return gene2go


def bh(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        val = min(prev, pvals[i] * n / (n - rank + 1))
        q[i] = prev = val
    return q


def read(path):
    if not os.path.exists(path):
        return set()
    return {l.strip() for l in open(path) if l.strip()}


def enrich(study, bg, gene2go, name, top):
    term_bg, term_st = defaultdict(set), defaultdict(set)
    for g in bg:
        for t in gene2go.get(g, ()):
            term_bg[t].add(g)
    for g in study:
        for t in gene2go.get(g, ()):
            term_st[t].add(g)
    M, N = len(bg), len(study)
    rows = []
    for t, gs in term_bg.items():
        n = len(gs)
        k = len(term_st.get(t, ()))
        if n < MIN_TERM or n > MAX_TERM or k < MIN_HIT:
            continue
        exp = n * N / M
        if k <= exp:
            continue
        rows.append((hypergeom.sf(k - 1, M, n, N), t, k, n, k / exp))
    if not rows:
        return []
    qs = bh([r[0] for r in rows])
    out = [(q, p, t, k, n, f) for q, (p, t, k, n, f) in zip(qs, rows)]
    out.sort(key=lambda r: (r[1], -r[5]))
    return out, len(rows)


def main():
    args = [a for a in sys.argv[1:]]
    aspect, top, pcut = "P", 12, None
    if "--aspect" in args:
        i = args.index("--aspect"); aspect = args[i + 1]; del args[i:i + 2]
    if "--top" in args:
        i = args.index("--top"); top = int(args[i + 1]); del args[i:i + 2]
    if "--p" in args:
        i = args.index("--p"); pcut = float(args[i + 1]); del args[i:i + 2]

    name, ns, parents = parse_obo(OBO)
    sys.setrecursionlimit(10000)
    gene2go = load_annotations(GAF, aspect, parents)
    print(f"aspect {aspect}   annotated genes {len(gene2go)}\n")

    for hm in args:
        mark = os.path.basename(os.path.dirname(hm.rstrip("/")))
        classes = {c: read(os.path.join(hm, f"{c}.genes_pc.txt"))
                   for c in ("4FP_specific", "shared", "CP_specific")}
        bg = set().union(*classes.values())
        if not bg:
            print(f"{mark}: no annotated gene lists\n"); continue
        print("=" * 96)
        print(f"{mark.upper()}   background {len(bg)} protein-coding genes with a peak")
        for cls in ("4FP_specific", "CP_specific"):
            st = classes[cls]
            print(f"\n  {cls}  ({len(st)} genes)")
            res, tested = enrich(st, bg, gene2go, name, top)
            if not res:
                print("    nothing enriched")
                continue
            if pcut is not None:
                hits = [r for r in res if r[1] < pcut]
                exp_fp = tested * pcut
                print(f"    {len(hits)} of {tested} terms at p < {pcut}; "
                      f"{exp_fp:.0f} expected by chance, so ~{100*exp_fp/max(1,len(hits)):.0f}% of this list is noise")
                sig = sum(1 for r in res if r[0] < 0.05)
                print(f"    {sig} survive FDR 0.05")
                res = hits[:top]
            else:
                res = [r for r in res if r[0] < 0.05][:top] or res[:top]
            print(f"    {'p':>9} {'FDR':>9} {'fold':>5} {'k/n':>10}  term")
            for q, p, t, k, n, f in res:
                flag = "  " if q < 0.05 else " ."
                print(f"   {flag}{p:>8.1e} {q:>9.1e} {f:>5.2f} {k:>4}/{n:<5}  {name.get(t, t)[:56]}")
        print()


if __name__ == "__main__":
    main()
