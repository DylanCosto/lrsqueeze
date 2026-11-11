"""Compare two sorted bedMethyl files (modkit columns) site by site, streaming.

Sites are matched by (chrom, start, mod code, strand) with a merge join over files in the same
sort order. For each mod code, at sites with valid coverage >= MIN_COV in both, reports the mean
and 95th-percentile absolute change in percent modified, the share of sites changing by more
than 10 points, agreement of >= 50% calls, Pearson correlation, and total valid coverage.

usage: python -I compare_bedmethyl_stream.py A.bed.gz B.bed.gz [MIN_COV]
"""
import collections
import gzip
import json
import math
import sys


def rows(path):
    with gzip.open(path, 'rt') as f:
        for line in f:
            c = line.split()
            yield (c[0], int(c[1]), c[3], c[5]), int(c[9]), float(c[10])


def main():
    a_path, b_path = sys.argv[1], sys.argv[2]
    min_cov = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    chrom_order = {}

    def order(key):
        rank = chrom_order.setdefault(key[0], len(chrom_order))
        return (rank, key[1], key[2], key[3])

    stats = collections.defaultdict(lambda: dict(n=0, sad=0.0, gt10=0, agree=0, sx=0.0, sy=0.0,
                                                 sxx=0.0, syy=0.0, sxy=0.0, hist=[0] * 1001,
                                                 cov_a=0, cov_b=0, only_a=0, only_b=0))
    ia, ib = rows(a_path), rows(b_path)
    ra, rb = next(ia, None), next(ib, None)
    while ra or rb:
        if rb is None or (ra and order(ra[0]) < order(rb[0])):
            stats[ra[0][2]]['only_a'] += 1; stats[ra[0][2]]['cov_a'] += ra[1]; ra = next(ia, None); continue
        if ra is None or order(rb[0]) < order(ra[0]):
            stats[rb[0][2]]['only_b'] += 1; stats[rb[0][2]]['cov_b'] += rb[1]; rb = next(ib, None); continue
        s = stats[ra[0][2]]
        s['cov_a'] += ra[1]; s['cov_b'] += rb[1]
        if ra[1] >= min_cov and rb[1] >= min_cov:
            x, y = ra[2], rb[2]; d = abs(x - y)
            s['n'] += 1; s['sad'] += d; s['gt10'] += d > 10; s['agree'] += (x >= 50) == (y >= 50)
            s['sx'] += x; s['sy'] += y; s['sxx'] += x * x; s['syy'] += y * y; s['sxy'] += x * y
            s['hist'][min(int(d * 10), 1000)] += 1
        ra, rb = next(ia, None), next(ib, None)
    out = {}
    for code, s in stats.items():
        n = s['n']
        if not n:
            continue
        cov = s['sxy'] / n - (s['sx'] / n) * (s['sy'] / n)
        var = math.sqrt(max(s['sxx'] / n - (s['sx'] / n) ** 2, 0) * max(s['syy'] / n - (s['sy'] / n) ** 2, 0))
        acc, p95 = 0, None
        for i, h in enumerate(s['hist']):
            acc += h
            if acc >= 0.95 * n:
                p95 = i / 10
                break
        out[code] = dict(sites=n, pearson=cov / var if var else None, mean_abs_diff_pp=s['sad'] / n,
                         p95_abs_diff_pp=p95, frac_diff_gt10pp=s['gt10'] / n, call50_agreement=s['agree'] / n,
                         valid_cov_ratio=s['cov_b'] / s['cov_a'] if s['cov_a'] else None,
                         sites_only_a=s['only_a'], sites_only_b=s['only_b'])
    print(json.dumps(out))


if __name__ == '__main__':
    main()
