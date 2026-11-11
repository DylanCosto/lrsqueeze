#!/usr/bin/env python3
"""snv_concordance.py NAME ORIGINAL.vcf.gz ARCHIVE.vcf.gz — PASS genotype concordance, SNPs and indels.

Streams both VCFs (bcftools query, PASS non-reference genotypes) side by side in their shared contig
order, so memory stays constant. Prints one TSV line: name, then per class (snp, indel):
both_same_gt, both_diff_gt, only_original, only_archive, concordance (same / union).
"""
import subprocess
import sys

FMT = "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n"
INC = 'FILTER="PASS" && GT!="RR" && GT!="mis"'


def contigs(*paths):
    """Contig order: the first VCF's, then any contig only a later one declares."""
    order = {}
    for path in paths:
        out = subprocess.run(["bcftools", "view", "-h", path], capture_output=True, text=True, check=True).stdout
        for l in out.splitlines():
            if l.startswith("##contig="):
                order.setdefault(l.split("ID=", 1)[1].split(",", 1)[0].rstrip(">"), len(order))
    return order


def records(path, order):
    p = subprocess.Popen(["bcftools", "query", "-i", INC, "-f", FMT, path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        c, pos, ref, alt, gt = line.rstrip("\n").split("\t")
        # Unordered: which haplotype is 1 or 2 in a phase block is arbitrary, so 0|1 equals 1|0.
        yield (order[c], int(pos), ref, alt), "/".join(sorted(gt.replace("|", "/").split("/")))
    if p.wait():
        sys.exit(f"bcftools query failed on {path}")


def kind(key):
    return "snp" if len(key[2]) == 1 and all(len(a) == 1 for a in key[3].split(",")) else "indel"


def by_position(it):
    """Group consecutive records at the same (contig, pos) into a dict key -> gt."""
    group, at = {}, None
    for key, gt in it:
        if key[:2] != at and group:
            yield at, group
            group = {}
        at = key[:2]
        group[key] = gt
    if group:
        yield at, group


name, a, b = sys.argv[1:4]
order = contigs(a, b)
n = {k: [0, 0, 0, 0] for k in ("snp", "indel")}  # same, diff, only original, only archive
ga, gb = by_position(records(a, order)), by_position(records(b, order))
pa, pb = next(ga, None), next(gb, None)
while pa or pb:
    if pb is None or (pa and pa[0] < pb[0]):
        for k in pa[1]:
            n[kind(k)][2] += 1
        pa = next(ga, None)
    elif pa is None or pb[0] < pa[0]:
        for k in pb[1]:
            n[kind(k)][3] += 1
        pb = next(gb, None)
    else:
        da, db = pa[1], pb[1]
        for k in da.keys() | db.keys():
            c = n[kind(k)]
            if k in da and k in db:
                c[0 if da[k] == db[k] else 1] += 1
            else:
                c[2 if k in da else 3] += 1
        pa, pb = next(ga, None), next(gb, None)
out = [name]
for k in ("snp", "indel"):
    same, diff, oo, oa = n[k]
    out += [same, diff, oo, oa, f"{same / max(1, same + diff + oo + oa):.6f}"]
print("\t".join(map(str, out)), flush=True)
