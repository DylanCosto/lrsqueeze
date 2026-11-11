"""Independent check of lrsqueeze's ONT policies (bins + homopolymer protection) on real data.

Applies the policy to every source record in Python, with the homopolymer mask and projection
of the study the ONT policies were validated in (`protection()` below), and compares each decoded record
(CRAM or BAM) field by field. Also prints the SHA-256 of all expected quality bytes in source
order (for --bins q4 this is the study's hybrid4 checksum, 036d2739...).

  --bins q8 --keep-secondary   the default mode (8 bins, protection, secondaries kept; MD/NM
                               dropped from records without SEQ)
  --bins q4                    --max (4 bins, protection, secondaries dropped)

usage: python -I ont_hybrid.py SOURCE.bam OUTPUT REF.fa [--bins q4|q8] [--keep-secondary]
       [--unmapped read|binned]
"""
import argparse
import hashlib
import json
import re

import numpy as np
import pysam

p = argparse.ArgumentParser()
p.add_argument('source')
p.add_argument('output')
p.add_argument('ref')
p.add_argument('--bins', choices=['q4', 'q8'], default='q4')
p.add_argument('--keep-secondary', action='store_true')
p.add_argument('--unmapped', choices=['read', 'binned'], default='read',
               help="unmapped reads: native qualities in runs of their own bases (lrsqueeze 0.3.0's "
                    "homopolymers_in = reference-or-read), or binned throughout (0.2.1 and earlier)")
args = p.parse_args()

BINS = {
    'q4': [(0, 9, 7), (10, 19, 15), (20, 29, 27), (30, 93, 38)],
    'q8': [(0, 4, 4), (5, 9, 7), (10, 14, 12), (15, 19, 17), (20, 24, 23), (25, 29, 28), (30, 34, 32), (35, 93, 40)],
}[args.bins]
LUT = np.zeros(256, dtype=np.uint8)
for lo, hi, v in BINS:
    LUT[lo:hi + 1] = v
LUT[94:] = BINS[-1][2]
ML16 = np.array([(v // 16) * 16 + 8 for v in range(256)], dtype=np.uint8)

masks = {}


def protection(r, hp, starts, ends):
    """Query bases kept native: aligned (M/=/X) to a reference run in [starts, ends), or inserted
    next to one (`hp` is the run mask). Soft clips and unmapped reads are never protected."""
    keep = np.zeros(r.query_length, dtype=bool)
    if r.is_unmapped:
        return keep
    rp, qp = r.reference_start, 0
    j = int(np.searchsorted(ends, rp, side='right'))
    for op, n in r.cigartuples:
        if op in (0, 7, 8):
            stop = rp + n
            while j < len(ends) and ends[j] <= rp:
                j += 1
            while j < len(starts) and starts[j] < stop:
                lo, hi = max(rp, int(starts[j])), min(stop, int(ends[j]))
                if hi > lo:
                    keep[qp + lo - rp:qp + hi - rp] = True
                if ends[j] > stop:
                    break
                j += 1
            rp, qp = stop, qp + n
        elif op == 1:
            if (rp > 0 and hp[rp - 1]) or (rp < len(hp) and hp[rp]):
                keep[qp:qp + n] = True
            qp += n
        elif op in (2, 3):
            rp += n
        elif op == 4:
            qp += n
        elif op not in (5, 6):
            raise ValueError(f'unexpected CIGAR operation {op}')
    assert qp == r.query_length
    return keep


def hp_mask(name):
    if name not in masks:
        seq = pysam.FastaFile(args.ref).fetch(name).upper()
        m = np.zeros(len(seq), dtype=bool)
        for match in re.finditer(r'A{5,}|C{5,}|G{5,}|T{5,}', seq):
            m[match.start():match.end()] = True
        edges = np.diff(np.r_[False, m, False].astype(np.int8))
        masks[name] = (m, np.flatnonzero(edges == 1), np.flatnonzero(edges == -1))
    return masks[name]


def runs_mask(sequence):
    """Bases in runs of 5 or more of one of A, C, G, T (any case) in a read's own sequence."""
    m = np.zeros(len(sequence), dtype=bool)
    for match in re.finditer(r'A{5,}|C{5,}|G{5,}|T{5,}', sequence.upper()):
        m[match.start():match.end()] = True
    return m


def tags(r):
    return {k: (list(v) if not isinstance(v, (str, int, float)) else v) for k, v in r.get_tags()}


qual_hash = hashlib.sha256()
checked = skipped = 0
with pysam.AlignmentFile(args.source, check_sq=False) as src, \
        pysam.AlignmentFile(args.output, reference_filename=args.ref) as out:
    decoded = iter(out)
    for a in src:
        q = a.query_qualities
        if q is not None:
            q = np.asarray(q, dtype=np.uint8)
            mapped = LUT[q]
            if not a.is_unmapped:
                m, starts, ends = hp_mask(a.reference_name)
                keep = protection(a, m, starts, ends)
                mapped[keep] = q[keep]
            elif args.unmapped == 'read' and a.query_sequence:
                keep = runs_mask(a.query_sequence)
                mapped[keep] = q[keep]
            qual_hash.update(mapped.tobytes())
        if a.is_secondary and not args.keep_secondary:
            skipped += 1
            continue
        b = next(decoded)
        for field in ('query_name', 'flag', 'reference_id', 'reference_start', 'mapping_quality',
                      'cigarstring', 'next_reference_id', 'next_reference_start', 'template_length',
                      'query_sequence'):
            assert getattr(a, field) == getattr(b, field), (a.query_name, field)
        got = b.query_qualities
        assert (q is None and got is None) or np.array_equal(mapped, np.asarray(got, dtype=np.uint8)), \
            (a.query_name, 'QUAL')
        ta, tb = tags(a), tags(b)
        if 'ML' in ta:
            ta['ML'] = ML16[np.asarray(ta['ML'], dtype=np.uint8)].tolist()
        if a.query_sequence is None:
            ta.pop('MD', None)
            ta.pop('NM', None)
        assert ta.keys() == tb.keys(), (a.query_name, sorted(ta.keys() ^ tb.keys()))
        for k in ta:
            assert ta[k] == tb[k], (a.query_name, k, str(ta[k])[:60], str(tb[k])[:60])
        checked += 1
    assert next(decoded, None) is None, 'output has extra records'

print(json.dumps(dict(bins=args.bins, records_checked=checked, secondary_skipped=skipped,
                      quality_sha256=qual_hash.hexdigest())))
