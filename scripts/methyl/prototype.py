"""Prototype methylation-storage variants for measurement (not used by lrsqueeze itself).

Rewrites MM/ML on every record and writes BAM:
  no6ma   drop the A+a (6mA) part
  cpg     keep C-modification calls only at CpG sites (read CpG, or aligned to a reference CpG),
          as an explicit ('?') part: implicit ('.') unmodified CpGs become explicit 0 calls, so
          CpG pileups keep their coverage and non-CpG Cs read as "no call", not "unmodified"
  cpg_no6ma  both
  none    drop MM, ML (and MN)

usage: python -I prototype.py VARIANT IN.bam OUT.bam REF.fa
"""
import array
import sys

import pysam

COMP = str.maketrans('ACGTN', 'TGCAN')


def parts(mm, ml):
    """Yield (head, deltas, values-per-call list) for each MM part."""
    k = 0
    for part in mm.rstrip(';').split(';'):
        if not part:
            continue
        head, *deltas = part.split(',')
        code = head[2:].rstrip('.?')
        n = 1 if code.isdigit() else len(code)
        vals = [tuple(ml[k + i * n: k + (i + 1) * n]) for i in range(len(deltas))]
        k += n * len(deltas)
        yield head, [int(d) for d in deltas], n, vals


def cpg_part(head, deltas, n, vals, orig, is_cpg):
    """Rebuild a C part keeping CpG positions only, explicit mode."""
    implicit = not head.endswith('?')
    c_idx = [i for i, b in enumerate(orig) if b == 'C']
    called = {}
    idx = -1
    for d, v in zip(deltas, vals):
        idx += d + 1
        called[idx] = v
    out_d, out_v, gap = [], [], 0
    for j, p in enumerate(c_idx):
        if not is_cpg(p):
            gap += 1
            continue
        if j in called:
            v = called[j]
        elif implicit:
            v = (0,) * n
        else:
            gap += 1
            continue
        out_d.append(gap)
        out_v.extend(v)
        gap = 0
    return head.rstrip('.?') + '?', out_d, out_v


def transform(r, variant, ref):
    if variant == 'none':
        for t in ('MM', 'ML', 'MN'):
            if r.has_tag(t):
                r.set_tag(t, None)
        return
    if not (r.has_tag('MM') and r.has_tag('ML')) or r.query_sequence is None:
        return
    mm, ml = r.get_tag('MM'), list(r.get_tag('ML'))
    seq = r.query_sequence
    orig = seq[::-1].translate(COMP) if r.is_reverse else seq
    L = len(seq)
    ref_cpg = set()
    if variant in ('cpg', 'cpg_no6ma') and not r.is_unmapped:
        # One base of context on each side, so CpGs at either alignment end are seen.
        lo = max(r.reference_start - 1, 0)
        rseq = ref.fetch(r.reference_name, lo, r.reference_end + 1).upper()
        for q, rp in r.get_aligned_pairs(matches_only=True):
            o = rp - lo
            # Forward: read C on reference C of CG. Reverse: the original read's C is a G in
            # SEQ, on the reference G of CG.
            if (not r.is_reverse and rseq[o:o + 2] == 'CG') or (r.is_reverse and o > 0 and rseq[o - 1:o + 1] == 'CG'):
                ref_cpg.add(L - 1 - q if r.is_reverse else q)

    def is_cpg(p):
        return orig[p + 1:p + 2] == 'G' or p in ref_cpg

    new_mm, new_ml = [], []
    for head, deltas, n, vals in parts(mm, ml):
        base = head[0]
        if base == 'A' and variant in ('no6ma', 'cpg_no6ma'):
            continue
        if base == 'C' and variant in ('cpg', 'cpg_no6ma'):
            head, deltas, flat = cpg_part(head, deltas, n, vals, orig, is_cpg)
        else:
            flat = [x for v in vals for x in v]
        new_mm.append(','.join([head] + [str(d) for d in deltas]))
        new_ml.extend(flat)
    r.set_tag('MM', ';'.join(new_mm) + ';' if new_mm else None, 'Z' if new_mm else None)
    r.set_tag('ML', array.array('B', new_ml) if new_mm else None)


def main():
    variant, src, dst, ref_path = sys.argv[1:]
    ref = pysam.FastaFile(ref_path)
    with pysam.AlignmentFile(src) as a, pysam.AlignmentFile(dst, 'wb', template=a, threads=4) as b:
        for r in a:
            transform(r, variant, ref)
            b.write(r)
    pysam.index(dst)


if __name__ == '__main__':
    main()
