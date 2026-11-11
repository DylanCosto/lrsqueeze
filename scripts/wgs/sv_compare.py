"""Concordance of two SV VCFs (e.g. original vs archive): exact and near matches, genotypes.

PASS records only. Exact: same CHROM, POS, SVTYPE, SVLEN, END. Near (for the rest): same CHROM
and SVTYPE, POS within 50 bp, |SVLEN| within 10% (greedy, one-to-one). Genotypes are compared
on matched pairs.

usage: python -I sv_compare.py BASE.vcf.gz CANDIDATE.vcf.gz
"""
import collections
import gzip
import json
import sys


def load(path):
    out = []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            c = line.rstrip('\n').split('\t')
            if c[6] != 'PASS':
                continue
            info = dict(kv.split('=', 1) if '=' in kv else (kv, True) for kv in c[7].split(';'))
            gt = c[9].split(':')[0] if len(c) > 9 else '.'
            out.append(dict(chrom=c[0], pos=int(c[1]), type=info.get('SVTYPE', '?'),
                            len=abs(int(info.get('SVLEN', 0))), end=info.get('END'), gt=gt))
    return out


def main():
    a, b = load(sys.argv[1]), load(sys.argv[2])
    key = lambda r: (r['chrom'], r['pos'], r['type'], r['len'], r['end'])
    exact_pairs, rest_a = [], []
    b_by_key = collections.defaultdict(list)
    for r in b:
        b_by_key[key(r)].append(r)
    for r in a:
        if b_by_key[key(r)]:
            exact_pairs.append((r, b_by_key[key(r)].pop()))
        else:
            rest_a.append(r)
    rest_b = [r for lst in b_by_key.values() for r in lst]
    idx = collections.defaultdict(list)
    for r in rest_b:
        idx[(r['chrom'], r['type'])].append(r)
    near_pairs, used = [], set()
    for r in rest_a:
        best = None
        for s in idx[(r['chrom'], r['type'])]:
            if id(s) in used or abs(s['pos'] - r['pos']) > 50:
                continue
            if max(r['len'], s['len']) and abs(r['len'] - s['len']) > 0.1 * max(r['len'], s['len']):
                continue
            if best is None or abs(s['pos'] - r['pos']) < abs(best['pos'] - r['pos']):
                best = s
        if best is not None:
            used.add(id(best))
            near_pairs.append((r, best))
    pairs = exact_pairs + near_pairs
    # Unordered and unphased: 1|0, 0|1 and 0/1 are the same genotype.
    unordered = lambda gt: sorted(gt.replace('|', '/').split('/'))
    gt_same = sum(unordered(x['gt']) == unordered(y['gt']) for x, y in pairs)
    by_type = collections.Counter(r['type'] for r in a)
    print(json.dumps(dict(
        base=sys.argv[1], candidate=sys.argv[2], base_pass=len(a), candidate_pass=len(b),
        exact=len(exact_pairs), near=len(near_pairs),
        only_base=len(a) - len(pairs), only_candidate=len(b) - len(pairs),
        matched_fraction_of_base=len(pairs) / len(a) if a else None,
        genotype_agreement=gt_same / len(pairs) if pairs else None,
        base_types=dict(by_type))))


if __name__ == '__main__':
    main()
