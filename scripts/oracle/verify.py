"""Record-by-record check: decoded CRAM == transform(source BAM), with regenerated MD/NM.

The expected records come from the Python reference transform (`transform.py`), independent of
lrsqueeze's Rust code. CRAM stores alignments as read features, so `=`/`X` CIGAR operations decode
as `M` (as with samtools); CIGARs are compared after the same normalization. A decoder may add
MD/NM to records whose source had none; when the source has them, they must match.

usage: python -I verify.py SOURCE.bam OUTPUT.cram REF.fa [transform.py options...]
"""
import hashlib
import json
import pathlib
import sys

import numpy as np
import pysam

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from transform import make_transform, parser  # noqa: E402

REGENERATED = {'MD', 'NM'}


def cigar_m(record):
    """CIGAR with =/X folded into M and adjacent equal operations merged."""
    out = []
    for op, n in record.cigartuples or []:
        op = 0 if op in (7, 8) else op
        if out and out[-1][0] == op:
            out[-1][1] += n
        else:
            out.append([op, n])
    return out


src_path, cram_path, ref, *targs = sys.argv[1:]
args = parser().parse_args([src_path, *targs])
transform = make_transform(args)
n = skipped = 0
digest = hashlib.sha256()
with pysam.AlignmentFile(src_path, check_sq=False) as a, \
        pysam.AlignmentFile(cram_path, reference_filename=ref) as b:
    source = (r for r in a if not (args.drop_secondary and r.is_secondary))
    for ra, rb in zip(source, b, strict=True):
        ra = transform(ra)
        for field in ('query_name', 'flag', 'reference_id', 'reference_start', 'mapping_quality',
                      'next_reference_id', 'next_reference_start', 'template_length', 'query_sequence'):
            assert getattr(ra, field) == getattr(rb, field), (n, ra.query_name, field)
        assert cigar_m(ra) == cigar_m(rb), (n, ra.query_name, 'CIGAR')
        qa, qb = ra.query_qualities, rb.query_qualities
        assert (qa is None and qb is None) or list(qa) == list(qb), (n, ra.query_name, 'QUAL')
        ta = {k: (list(v) if not isinstance(v, (str, int, float)) else v) for k, v in ra.get_tags()}
        tb = {k: (list(v) if not isinstance(v, (str, int, float)) else v) for k, v in rb.get_tags()}
        for k in ta.keys() | tb.keys():
            if k in REGENERATED and k not in ta:
                continue
            va, vb = ta.get(k), tb.get(k)
            ok = np.allclose(va, vb, rtol=0, atol=0) if isinstance(va, float) else va == vb
            assert ok, (n, ra.query_name, k, str(va)[:80], str(vb)[:80])
        digest.update(rb.to_string().encode())
        n += 1
print(json.dumps(dict(cram=cram_path, records=n, decoded_sha256=digest.hexdigest(),
                      status='identical-to-expected')))
