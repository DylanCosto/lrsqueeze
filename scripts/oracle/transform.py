"""Lossy transforms for long-read BAM records, streamed as uncompressed BAM to stdout.

Options compose: quality removal/LUT, dropping low-probability calls from
implicit ('.') MM parts, and uniform ML quantization. Other fields are untouched.
"""
import argparse
import array
import functools
import operator

import numpy as np
import pysam


def ml_lut(levels):
    """Map 0..255 to the midpoint of one of `levels` equal-width bins."""
    if not levels:
        return None
    width = 256 // levels
    v = np.arange(256)
    return np.minimum((v // width) * width + width // 2, 255).astype(np.uint8)


def rewrite_mods(mm, ml, drop_below, lut):
    """Return (mm, ml) after dropping low calls in single-mod '.' parts and quantizing."""
    ml = np.asarray(ml, dtype=np.uint8)
    out_parts, out_ml, k = [], [], 0
    for part in mm.rstrip(';').split(';'):
        head, *deltas = part.split(',')
        code = head[2:].rstrip('.?')
        nmods = 1 if code.isdigit() else len(code)
        vals = ml[k:k + nmods * len(deltas)].reshape(-1, nmods)
        k += nmods * len(deltas)
        if drop_below and nmods == 1 and head[-1] != '?':
            # Implicit mode ('.' or unmarked): an omitted base means "not modified", so a skipped
            # low-probability call is folded into the next delta.
            keep_d, keep_v, carry = [], [], 0
            for d, v in zip(deltas, vals[:, 0]):
                if v < drop_below:
                    carry += int(d) + 1
                else:
                    keep_d.append(str(int(d) + carry)); keep_v.append(v); carry = 0
            deltas, vals = keep_d, np.asarray(keep_v, dtype=np.uint8).reshape(-1, 1)
        out_parts.append(','.join([head] + deltas))
        out_ml.append(vals.reshape(-1))
    assert k == len(ml), 'ML length does not match MM'
    new_ml = np.concatenate(out_ml) if out_ml else ml[:0]
    if lut is not None:
        new_ml = lut[new_ml]
    return ';'.join(out_parts) + ';', new_ml


def rq_phred(rq, top):
    """lrsqueeze's Phred score of a read accuracy: error clamped to 1e-4, rounded half to even."""
    return max(0, min(top, round(-10 * np.log10(max(1 - rq, 1e-4)))))


def mean_error_phred(quals, top):
    """The Phred score of the mean error probability, summed in order as lrsqueeze does."""
    if not quals:
        return top
    mean = functools.reduce(operator.add, (10.0 ** (-q / 10.0) for q in quals), 0.0) / len(quals)
    return max(0, min(top, round(-10 * np.log10(mean))))


def mean_phred(quals, top):
    """The mean of the Phred scores, rounded half to even."""
    return min(top, round(sum(quals) / len(quals))) if quals else top


def make_transform(args):
    qlut = None
    if args.qual_lut:
        qlut = np.arange(256, dtype=np.uint8)
        for i, q in enumerate(map(int, args.qual_lut.split(','))):
            qlut[i] = q
        qlut[len(args.qual_lut.split(',')):94] = qlut[len(args.qual_lut.split(',')) - 1]
    lut = ml_lut(args.ml_levels)
    drop = int(round(args.ml_drop * 256)) if args.ml_drop else 0

    def transform(r):
        if args.qual == 'none':
            r.query_qualities = None
        elif args.qual in ('rq', 'mean-error', 'mean-phred') and r.query_qualities is not None:
            # One score per read: from PacBio's read-accuracy tag (reads without it, and
            # mean-error, take their mean error probability's), or the mean Phred score.
            quals = list(r.query_qualities)
            if args.qual == 'rq' and r.has_tag('rq'):
                q = rq_phred(r.get_tag('rq'), args.qual_max)
            elif args.qual in ('rq', 'mean-error'):
                q = mean_error_phred(quals, args.qual_max)
            else:
                q = mean_phred(quals, args.qual_max)
            r.query_qualities = [q] * r.query_length
        elif qlut is not None and r.query_qualities is not None:
            r.query_qualities = qlut[np.asarray(r.query_qualities, dtype=np.uint8)].tolist()
        for t in args.drop_tags:
            if r.has_tag(t):
                r.set_tag(t, None)
        if args.strip_md_noseq and r.query_sequence is None:
            # Without SEQ, CRAM cannot regenerate MD/NM and stores them verbatim.
            for t in ('MD', 'NM'):
                if r.has_tag(t):
                    r.set_tag(t, None)
        if (drop or lut is not None) and r.has_tag('MM') and r.has_tag('ML'):
            mm, ml = rewrite_mods(r.get_tag('MM'), r.get_tag('ML'), drop, lut)
            r.set_tag('MM', mm, 'Z')
            r.set_tag('ML', array.array('B', ml.tobytes()))
        return r
    return transform


def parser():
    p = argparse.ArgumentParser()
    p.add_argument('src')
    p.add_argument('--qual', choices=['keep', 'none', 'rq', 'mean-error', 'mean-phred'],
                   default='keep', help='per-read sources as lrsqueeze names them (PacBio --max: mean-phred)')
    p.add_argument('--qual-max', type=int, default=40, help='highest per-read score')
    p.add_argument('--qual-lut', help='comma-separated output Q for input Q=0..n (last value extends to Q93)')
    p.add_argument('--ml-drop', type=float, default=0.0, help="drop '.'-mode calls with p below this")
    p.add_argument('--ml-levels', type=int, default=0, help='uniform ML quantization levels (0=off)')
    p.add_argument('--drop-tags', nargs='*', default=[])
    p.add_argument('--strip-md-noseq', action='store_true', help='drop MD/NM on records without SEQ')
    p.add_argument('--drop-secondary', action='store_true', help='omit secondary alignments')
    return p


if __name__ == '__main__':
    args = parser().parse_args()
    transform = make_transform(args)
    with pysam.AlignmentFile(args.src, check_sq=False) as src, \
            pysam.AlignmentFile('-', 'wbu', template=src) as out:
        for r in src:
            if args.drop_secondary and r.is_secondary:
                continue
            out.write(transform(r))
