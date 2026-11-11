"""How quality binning moves each read's mean Q (mean error probability, as chopper/NanoFilt)."""
import sys
import numpy as np, pysam
path, n, *luts = sys.argv[1], int(sys.argv[2]), *sys.argv[3:]
def lut(spec):
    t = np.zeros(256, np.uint8)
    for item in spec.split(','):
        r, v = item.split(':'); lo, hi = map(int, r.split('-')); t[lo:hi + 1] = int(v)
    t[94:] = t[93]; return t
tables = {s: lut(s) for s in luts}
err = 10 ** (-np.arange(256) / 10)
rows = {k: [] for k in ['orig', *luts]}
with pysam.AlignmentFile(path, check_sq=False) as f:
    for r in f:
        if r.is_secondary or r.is_supplementary or r.query_qualities is None: continue
        q = np.asarray(r.query_qualities, np.uint8)
        rows['orig'].append(-10 * np.log10(err[q].mean()))
        for s, t in tables.items(): rows[s].append(-10 * np.log10(err[t[q]].mean()))
        if len(rows['orig']) >= n: break
o = np.array(rows['orig'])
print(f"{path.split('/')[-1]}: {len(o)} primary reads, original mean-Q median {np.median(o):.2f}")
for s in luts:
    b = np.array(rows[s]); d = b - o
    cross = {c: int(((o >= c) != (b >= c)).sum()) for c in (7, 10, 15, 20, 25, 30)}
    print(f"  bins {s}\n    mean-Q shift: median {np.median(d):+.2f}, 5th..95th pct {np.percentile(d,5):+.2f}..{np.percentile(d,95):+.2f}")
    print("    reads crossing a cutoff: " + ", ".join(f"Q{c}: {v} ({100*v/len(o):.2f}%)" for c, v in cross.items()))
