"""Are TF genes over-represented among the 4FP-specific regions?

Background is the other peak classes of the SAME mark, not the whole genome.
Peaks land on promoters far more often than chance, so a genome background
would report TF enrichment for any promoter-biased peak set and say nothing
about 4FP. Comparing 4FP-specific against shared + CP-specific asks the actual
question: given that a region carries this mark, is being 4FP-specific
associated with sitting at a TF gene?

Usage: 06_tf_enrich.py TF_LIST HEATMAP_DIR [HEATMAP_DIR ...]
"""
import os
import sys
from math import lgamma


def log_choose(n, k):
    if k < 0 or k > n:
        return float('-inf')
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def hyper_sf(k, M, n, N):
    """P(X >= k): M population, n successes in population, N drawn."""
    total = float('-inf')
    hi = min(n, N)
    denom = log_choose(M, N)
    acc = 0.0
    for i in range(k, hi + 1):
        lp = log_choose(n, i) + log_choose(M - n, N - i) - denom
        acc += 2.718281828459045 ** lp if lp > -700 else 0.0
    return min(1.0, max(0.0, acc))


def read_genes(path):
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def main():
    tf = read_genes(sys.argv[1])
    print(f"TF reference: {len(tf)} genes\n")
    header = f"{'MARK':<10} {'CLASS':<14} {'GENES':>6} {'TF':>5} {'TF%':>6} {'BG_TF%':>7} {'FOLD':>6} {'P':>10}"
    print(header)
    print("-" * len(header))
    for hm in sys.argv[2:]:
        mark = os.path.basename(os.path.dirname(hm.rstrip('/')))
        sets = {}
        for label, fn in (('4FP_specific', '4FP_specific.genes_pc.txt'),
                          ('Shared', 'shared.genes_pc.txt'),
                          ('CP_specific', 'CP_specific.genes_pc.txt')):
            g = read_genes(os.path.join(hm, fn))
            if g:
                sets[label] = g
        if '4FP_specific' not in sets:
            print(f"{mark:<10} (not annotated yet)")
            continue
        for label, genes in sets.items():
            others = set().union(*[v for k, v in sets.items() if k != label]) - genes
            M = len(genes) + len(others)
            n = len((genes | others) & tf)
            N = len(genes)
            k = len(genes & tf)
            exp = n / M * N if M else 0
            fold = (k / exp) if exp else float('nan')
            p = hyper_sf(k, M, n, N) if k > exp else 1.0
            print(f"{mark:<10} {label:<14} {N:>6} {k:>5} {100*k/N:>5.1f}% "
                  f"{100*n/M:>6.1f}% {fold:>6.2f} {p:>10.2e}")
        print()


if __name__ == '__main__':
    main()
