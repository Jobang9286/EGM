#!/usr/bin/env python
"""Do the two arm-specific gene sets differ from EACH OTHER?

Usage: 08_go_contrast.py HEATMAP_DIR [HEATMAP_DIR ...] [--top N] [--p P]

07_go_enrich.py tests each class against all genes carrying the mark, which
answers "is this class enriched for X". Both classes can come back enriched for
the same X, and then the contrast between the arms is empty even though both
lists look full. This script tests the classes against each other instead:
study is the genes specific to 4FP alone, background is those plus the genes
specific to CP alone, so a term only scores when the two arms actually differ.

Genes appearing in both classes are dropped, not assigned to one side.
"""
import os
import sys
from collections import defaultdict

from scipy.stats import fisher_exact

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib.machinery import SourceFileLoader
_e = SourceFileLoader("e", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "07_go_enrich.py")).load_module()


def main():
    args = list(sys.argv[1:])
    top, pcut = 12, 0.05
    if "--top" in args:
        i = args.index("--top"); top = int(args[i + 1]); del args[i:i + 2]
    if "--p" in args:
        i = args.index("--p"); pcut = float(args[i + 1]); del args[i:i + 2]

    name, ns, parents = _e.parse_obo(_e.OBO)
    sys.setrecursionlimit(10000)
    gene2go = _e.load_annotations(_e.GAF, "P", parents)

    for hm in args:
        mark = os.path.basename(os.path.dirname(hm.rstrip("/")))
        a = _e.read(os.path.join(hm, "4FP_specific.genes_pc.txt"))
        b = _e.read(os.path.join(hm, "CP_specific.genes_pc.txt"))
        both = a & b
        A, B = a - both, b - both
        print("=" * 92)
        print(f"{mark.upper()}   4FP-only {len(A)}   CP-only {len(B)}   "
              f"in both {len(both)} (dropped)")
        if not A or not B:
            print("  one side empty\n"); continue

        ta, tb = defaultdict(int), defaultdict(int)
        for g in A:
            for t in gene2go.get(g, ()):
                ta[t] += 1
        for g in B:
            for t in gene2go.get(g, ()):
                tb[t] += 1

        rows = []
        for t in set(ta) | set(tb):
            ka, kb = ta.get(t, 0), tb.get(t, 0)
            n = ka + kb
            if n < 10:
                continue
            odds, p = fisher_exact([[ka, len(A) - ka], [kb, len(B) - kb]])
            rows.append((p, t, ka, kb, odds))
        if not rows:
            print("  no testable terms\n"); continue
        qs = _e.bh([r[0] for r in rows])
        res = sorted(zip(qs, rows), key=lambda r: r[1][0])

        nsig = sum(1 for q, r in res if r[0] < pcut)
        nfdr = sum(1 for q, r in res if q < 0.05)
        exp = len(rows) * pcut
        print(f"  {nsig} of {len(rows)} terms at p < {pcut} ({exp:.0f} expected by chance); "
              f"{nfdr} survive FDR 0.05")
        if nfdr == 0 and nsig <= exp * 1.3:
            print("  -> the two classes are not functionally distinguishable\n")
            continue
        print(f"    {'p':>9} {'FDR':>9} {'OR':>6}  {'4FP':>5} {'CP':>5}  term")
        for q, (p, t, ka, kb, odds) in res[:top]:
            side = "4FP" if odds > 1 else "CP "
            flag = "  " if q < 0.05 else " ."
            print(f"   {flag}{p:>8.1e} {q:>9.1e} {odds:>6.2f}  {ka:>5} {kb:>5}  "
                  f"[{side}] {name.get(t, t)[:48]}")
        print()


if __name__ == "__main__":
    main()
