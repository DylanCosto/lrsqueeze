"""Read mean Q (mean error probability, as chopper/NanoFilt compute it) before and after lrsqueeze.

Reads the original BAM and an lrsqueeze archive in lockstep (the archive must keep secondaries, as
the default mode does), compares each primary read's mean Q, and counts reads that cross common
filter cutoffs.

usage: python -I meanq_pair.py ORIGINAL.bam ARCHIVE REF.fa [MAX_READS]
"""
import json
import sys

import numpy as np
import pysam

src, arc, ref = sys.argv[1:4]
limit = int(sys.argv[4]) if len(sys.argv) > 4 else 50_000
err = 10 ** (-np.arange(256) / 10)
o, b = [], []
with pysam.AlignmentFile(src, check_sq=False) as fa, pysam.AlignmentFile(arc, reference_filename=ref) as fb:
    for ra, rb in zip(fa, fb):
        assert ra.query_name == rb.query_name, 'archives must keep record order (default mode)'
        if ra.is_secondary or ra.is_supplementary or ra.query_qualities is None:
            continue
        o.append(-10 * np.log10(err[np.asarray(ra.query_qualities)].mean()))
        b.append(-10 * np.log10(err[np.asarray(rb.query_qualities)].mean()))
        if len(o) >= limit:
            break
o, b = np.array(o), np.array(b)
d = b - o
print(json.dumps(dict(
    reads=len(o), median_shift=float(np.median(d)),
    p5_p95=[float(np.percentile(d, 5)), float(np.percentile(d, 95))],
    crossing_pct={f'Q{c}': float(100 * ((o >= c) != (b >= c)).mean()) for c in (7, 10, 15, 20, 25)})))
