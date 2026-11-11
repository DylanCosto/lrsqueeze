"""Paired block bootstrap of the F1 difference between two vcfeval runs.

Counts TP (baseline and call), FP and FN per fixed genome block from each run's vcfeval output
VCFs, resamples blocks with replacement (the same blocks for both runs), and recomputes F1 from
the pooled counts of each draw. Prints the observed difference and its percentile intervals.

usage: python -I bootstrap.py BASE_DIR CANDIDATE_DIR [BLOCK_BP] [DRAWS]
"""
import collections
import gzip
import json
import sys

import numpy as np

FILES = {'tpb': 'tp-baseline.vcf.gz', 'tpc': 'tp.vcf.gz', 'fp': 'fp.vcf.gz', 'fn': 'fn.vcf.gz'}


def counts(run, block):
    c = collections.defaultdict(lambda: np.zeros(4, np.int64))
    for i, name in enumerate(FILES.values()):
        with gzip.open(f'{run}/{name}', 'rt') as f:
            for line in f:
                if line[0] != '#':
                    chrom, pos = line.split('\t', 2)[:2]
                    c[(chrom, int(pos) // block)][i] += 1
    return c


def f1(t):
    tpb, tpc, fp, fn = t[..., 0], t[..., 1], t[..., 2], t[..., 3]
    p = tpc / (tpc + fp)
    r = tpb / (tpb + fn)
    return 2 * p * r / (p + r)


def main():
    base, cand = sys.argv[1], sys.argv[2]
    block = int(sys.argv[3]) if len(sys.argv) > 3 else 10_000_000
    draws = int(sys.argv[4]) if len(sys.argv) > 4 else 20_000
    a, b = counts(base, block), counts(cand, block)
    keys = sorted(set(a) | set(b))
    A = np.array([a[k] for k in keys])
    B = np.array([b[k] for k in keys])
    rng = np.random.default_rng(20261108)
    idx = rng.integers(0, len(keys), size=(draws, len(keys)))
    d = f1(B[idx].sum(1)) - f1(A[idx].sum(1))
    observed = float(f1(B.sum(0)) - f1(A.sum(0)))
    print(json.dumps(dict(
        base=base, candidate=cand, blocks=len(keys), block_bp=block, draws=draws,
        delta_f1_pp=100 * observed,
        ci95_two_sided_pp=[100 * float(np.percentile(d, 2.5)), 100 * float(np.percentile(d, 97.5))],
        lower95_one_sided_pp=100 * float(np.percentile(d, 5)))))


if __name__ == '__main__':
    main()
