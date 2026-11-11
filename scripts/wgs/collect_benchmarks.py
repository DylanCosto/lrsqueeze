"""Collect the whole-genome validation outputs into the tables under benchmarks/.

Reads the work directory that validate_dataset.sh, sv.sh and the methylation runs wrote, renames
run labels to the published dataset and condition names, and writes:

  small_variants.tsv       vcfeval counts, precision, recall and F1 per run, scope and type
  bootstrap.tsv            paired block-bootstrap F1 changes with 95% intervals
  structural_variants.tsv  SV concordance against the original
  methylation.tsv          CpG pileup comparisons (modkit-style and pb-CpG-tools count mode)
  sizes.tsv                archive sizes, encode wall time and peak memory

usage: python -I collect_benchmarks.py WORK_DIR OUT_DIR
"""
import csv
import gzip
import json
import pathlib
import re
import sys

# Run-label prefixes, longest first, and the published dataset names.
DATASETS = [('ont_hac', 'ont_hac'), ('pb_revio', 'pb_revio'), ('ont', 'ont_sup'), ('pb', 'pb_hifi')]
# Conditions whose labels differ from the published names.
CONDITIONS = {
    'ont_sup': {
        # The shipped default: 8 bins with homopolymer protection.
        'q8hp': 'default',
        # 8 bins without homopolymer protection, the default before protection was added.
        'default': 'q8_no_hp',
        # The shipped --max; run 'max' is the same policy with lzma and bzip2 enabled (see sizes).
        'max.shipped': 'max',
        'cpgonly': 'default_cpg_only',
        'cpgonly.parallel': 'default_cpg_only',
    },
    # The shipped PacBio --max: one quality per read, the mean of its Phred scores.
    'pb_revio': {'max': 'max_rq', 'max_meanphred': 'max', 'max_meanerr': 'max_mean_error'},
    'pb_hifi': {'max_meanphred': 'max', 'max_meanerr': 'max_mean_error'},
}
# Each dataset's original BAM, as named in the work directory's eval/sizes.tsv.
ORIGINAL_NAMES = {'ont_sup': 'original.bam', 'ont_hac': 'ont_hac.original',
                  'pb_revio': 'pb_revio.original', 'pb_hifi': 'HG002.hifi.bonefishY.bam'}
# The encode runs that used the parallel record transform (as 0.1.0 and later do); the others
# ran it on one thread.
PARALLEL_TRANSFORM_RUNS = {'ont.cpgonly.parallel', 'ont.read5', 'ont_hac.read5', 'ont_hac.read4',
                           'ont_hac.nohp', 'pb_revio.max_meanerr', 'pb_revio.max_meanphred',
                           'pb.max_meanerr', 'pb.max_meanphred'}
SCOPES = {'primary', 'autosomes', 'regions', 'longhp'}


def label(run):
    """(dataset, condition) for a run label such as 'ont.q8hp' or 'ont_hac.l5_f2'."""
    for prefix, dataset in DATASETS:
        if run.startswith(prefix + '.'):
            cond = run[len(prefix) + 1:]
            cond = CONDITIONS.get(dataset, {}).get(cond, cond)
            m = re.fullmatch(r'l(\d+)_f(\d+)', cond)
            if m:
                cond = f'hp{m[1]}_flank{m[2]}'
            m = re.fullmatch(r'read(\d+)', cond)
            if m:
                cond = f'read_hp{m[1]}'
            if cond == 'nohp':
                cond = 'q8_no_hp'
            return dataset, cond
    raise ValueError(f'unknown run label {run}')


def split_eval(name):
    """(dataset, condition, scope, type) for an eval name such as 'ont.q8hp.longhp.indels'."""
    parts = name.split('.')
    kind, scope = parts[-1], parts[-2]
    assert scope in SCOPES, name
    return (*label('.'.join(parts[:-2])), scope, kind)


def f1(tpb, tpc, fp, fn):
    p = tpc / (tpc + fp) if tpc + fp else 0.0
    r = tpb / (tpb + fn) if tpb + fn else 0.0
    return p, r, 2 * p * r / (p + r) if p + r else 0.0


def small_variants(work):
    rows = []
    for summary in sorted(work.glob('eval/*/summary.txt')):
        name = summary.parent.name
        if name.split('.')[-2] not in SCOPES or '.screen.' in name:
            continue
        for line in summary.read_text().splitlines():
            f = line.split()
            if f and f[0] == 'None':
                tpb, tpc, fp, fn = map(int, f[1:5])
                p, r, f1_ = f1(tpb, tpc, fp, fn)
                rows.append([*split_eval(name), tpb, tpc, fp, fn, f'{p:.6f}', f'{r:.6f}', f'{f1_:.6f}'])
    return sorted(rows), ['dataset', 'condition', 'scope', 'type', 'tp_baseline', 'tp_call', 'fp',
                          'fn', 'precision', 'recall', 'f1']


def bootstrap(work):
    rows = []
    for path in sorted(work.glob('eval/bootstrap.*.jsonl')):
        for line in path.read_text().splitlines():
            d = json.loads(line)
            ds, base, scope, kind = split_eval(d['base'])
            ds2, cand, scope2, kind2 = split_eval(d['candidate'])
            assert (scope, kind) == (scope2, kind2), d
            lo, hi = d['ci95_two_sided_pp']
            rows.append([ds, cand, base, scope, kind, f"{d['delta_f1_pp']:+.4f}", f'{lo:+.4f}',
                         f'{hi:+.4f}', d['blocks'], d['block_bp'], d['draws']])
    return sorted(set(map(tuple, rows))), ['dataset', 'condition', 'against', 'scope', 'type',
                                           'delta_f1_pp', 'ci95_low_pp', 'ci95_high_pp', 'blocks',
                                           'block_bp', 'draws']


def vcf_records(path):
    with gzip.open(path, 'rt') as f:
        return [line for line in f if line[0] != '#']


def structural_variants(work):
    rows = []
    for line in (work / 'eval/sv_concordance.jsonl').read_text().splitlines():
        d = json.loads(line)
        ds, base = label(d['base'].removesuffix('.vcf.gz'))
        _, cand = label(d['candidate'].removesuffix('.vcf.gz'))
        a, b = vcf_records(work / 'sv' / d['base']), vcf_records(work / 'sv' / d['candidate'])
        rows.append([ds, cand, base, d['base_pass'], d['candidate_pass'], d['exact'],
                     d['only_base'], d['only_candidate'], len(a), len(b),
                     'yes' if a == b else 'no'])
    return rows, ['dataset', 'condition', 'against', 'pass_against', 'pass_condition',
                  'pass_matched', 'pass_only_against', 'pass_only_condition', 'records_against',
                  'records_condition', 'vcf_records_identical']


def pbcpg_count_mode(work):
    """Compare pb-CpG-tools count-mode pileups of the Revio original and default archive."""
    def rows(path):
        with gzip.open(path, 'rt') as f:
            for line in f:
                if line[0] != '#':
                    yield line.rstrip('\n').split('\t')
    sites = identical = 0
    for a, b in zip(rows(work / 'methyl/pbcpg.orig.combined.bed.gz'),
                    rows(work / 'methyl/pbcpg.default.combined.bed.gz'), strict=True):
        sites += 1
        # chrom, begin, end, mod_score (methylation %), type, cov, mod_count, unmod_count
        identical += a[:8] == b[:8]
    return sites, identical


def methylation(work):
    header = ['dataset', 'condition', 'against', 'method', 'modification', 'sites', 'pearson',
              'mean_abs_diff_pp', 'p95_abs_diff_pp', 'call50_agreement', 'valid_coverage_ratio',
              'identical_sites']
    rows = []
    pileup = 'pileup --cpg --combine-strands --filter-threshold C:0.8 (modkit 0.6.4)'
    for file, ds, cond, against in [('methyl.hg003.json', 'ont_sup', 'default', 'original'),
                                    ('methyl.pb_revio.json', 'pb_revio', 'default', 'original'),
                                    ('methyl.hg003.cpgonly_vs_default.json', 'ont_sup',
                                     'default_cpg_only', 'default')]:
        for mod, d in json.loads((work / 'eval' / file).read_text()).items():
            rows.append([ds, cond, against, pileup, {'m': '5mC', 'h': '5hmC'}[mod], d['sites'],
                         f"{d['pearson']:.6f}", f"{d['mean_abs_diff_pp']:.4f}",
                         d['p95_abs_diff_pp'], f"{d['call50_agreement']:.6f}",
                         f"{d['valid_cov_ratio']:.6f}", ''])
    # Identical sites: position, methylation %, coverage and counts all equal. Only the
    # average-probability columns can differ.
    sites, identical = pbcpg_count_mode(work)
    rows.append(['pb_revio', 'default', 'original', 'cpg --pileup-mode count (pb-CpG-tools 3.0.0)',
                 '5mC', sites, '', '', '', '', '', identical])
    return rows, header


def elapsed_seconds(text):
    parts = [float(x) for x in text.split(':')]
    return sum(v * 60 ** i for i, v in enumerate(reversed(parts)))


def logged_sizes(work):
    """Bytes per name from eval/sizes.tsv, where inputs and deleted archives were recorded (the
    last entry for a name wins)."""
    sizes = {}
    for line in (work / 'eval/sizes.tsv').read_text().splitlines():
        name, value = line.split('\t')
        if name != 'file':
            sizes[name.removesuffix('.cram')] = int(value)
    return sizes


def sizes(work):
    logged = logged_sizes(work)
    original = {ds: logged[name] for ds, name in ORIGINAL_NAMES.items()}
    rows = []
    for time in sorted(work.glob('logs/encode_*.time')):
        run = time.name.removeprefix('encode_').removesuffix('.time')
        archive = work / 'archives' / f'{run}.cram'
        t = time.read_text()
        wall = re.search(r'Elapsed \(wall clock\).*: (\S+)', t)[1]
        rss = int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', t)[1])
        command = re.search(r'Command being timed: "(.*)"', t)[1]
        # An archive deleted after use had its size logged.
        size = archive.stat().st_size if archive.exists() else logged[run]
        ds, cond = label(run)
        notes = ['lzma and bzip2 enabled'] if run == 'ont.max' else []
        if run not in PARALLEL_TRANSFORM_RUNS:
            notes.append('single-threaded transform')
        rows.append([ds, cond, '; '.join(notes), size, f'{100 * size / original[ds]:.2f}',
                     round(elapsed_seconds(wall)), round(rss / 1024 ** 2, 2),
                     '-@ 16' if '-@ 16' in command else ''])
    for ds, n in original.items():
        rows.append([ds, 'original', '', n, '100.00', '', '', ''])
    return sorted(rows, key=lambda r: (r[0], r[1] != 'original', r[1], r[2])), [
        'dataset', 'condition', 'note', 'bytes', 'percent_of_original', 'encode_seconds',
        'peak_rss_gib', 'threads']


def main():
    work, out = map(pathlib.Path, sys.argv[1:])
    out.mkdir(parents=True, exist_ok=True)
    for name, make in [('small_variants', small_variants), ('bootstrap', bootstrap),
                       ('structural_variants', structural_variants), ('methylation', methylation),
                       ('sizes', sizes)]:
        rows, header = make(work)
        with open(out / f'{name}.tsv', 'w', newline='') as f:
            w = csv.writer(f, delimiter='\t', lineterminator='\n')
            w.writerow(header)
            w.writerows(rows)
        print(f'{name}.tsv: {len(rows)} rows')


if __name__ == '__main__':
    main()
