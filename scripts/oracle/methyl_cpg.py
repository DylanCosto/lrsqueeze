"""Check `lrsqueeze pack --methylation cpg|cpg-only` against the Python prototype.

The expected records come from scripts/methyl/prototype.py (variant `cpg` or `cpg_no6ma`) on the
same input; their ML values are quantized to 16 levels as lrsqueeze's default does. MM and ML must
be identical record by record.

usage: python -I methyl_cpg.py PROTOTYPE_OUTPUT.bam LRSQUEEZE_OUTPUT.cram REF.fa
"""
import sys
import pysam
proto, out, ref = sys.argv[1:4]
q16 = [(v // 16) * 16 + 8 for v in range(256)]
n = same = 0; diffs = []
with pysam.AlignmentFile(proto) as a, pysam.AlignmentFile(out, reference_filename=ref) as b:
    for ra, rb in zip(a, b, strict=True):
        assert ra.query_name == rb.query_name
        n += 1
        ea = (ra.get_tag('MM') if ra.has_tag('MM') else None, [q16[v] for v in ra.get_tag('ML')] if ra.has_tag('ML') else None)
        eb = (rb.get_tag('MM') if rb.has_tag('MM') else None, list(rb.get_tag('ML')) if rb.has_tag('ML') else None)
        if ea == eb: same += 1
        elif len(diffs) < 3: diffs.append((ra.query_name, str(ea[0])[:80], str(eb[0])[:80]))
print(dict(records=n, identical_MM_ML=same, differing=n - same, examples=diffs))
