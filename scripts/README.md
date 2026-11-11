# Scripts

Tooling behind the published numbers. None of it is built or shipped with lrsqueeze, and none of
it is needed to use lrsqueeze.

## `oracle/`: independent checks

Python implementations of lrsqueeze's policies that read archives through htslib (pysam), so an
archive is checked against code that shares nothing with lrsqueeze or pasta. Needs Python 3 with
pysam and numpy.

- `verify.py`: quality bins, one quality per read (`--qual mean-phred|mean-error|rq`), `ML`
  levels, `MD`/`NM`, secondary alignments.
- `ont_hybrid.py`: the ONT homopolymer rule, for mapped and unmapped reads.
- `methyl_cpg.py`, `transform.py`: the `--methylation cpg` rewrite, and the record transforms the
  others share.

```sh
python -I scripts/oracle/verify.py in.bam out.cram ref.fa --qual mean-phred --ml-levels 16 --strip-md-noseq --drop-secondary
python -I scripts/oracle/ont_hybrid.py in.bam out.cram ref.fa --bins q8 --keep-secondary
```

## `wgs/`: whole-genome validation

Encodes the validation datasets, calls variants on the originals and the archives, scores them,
and collects the tables in [benchmarks/](../benchmarks/README.md). The method and results are in
[docs/validation-wgs.md](../docs/validation-wgs.md).

Paths and tools come from the environment or from `scripts/wgs/local.env` (untracked): copy
[`local.env.example`](wgs/local.env.example) and fill it in. The scripts need `rtg` (RTG Tools),
`bcftools`, `tabix`, `samtools`, lrpipe for calling (a Rust port of Clair3 and Sniffles2),
Clair3 models, and Python 3 with numpy.

- `validate.sh`, `validate_dataset.sh`: encode, call (lrpipe snv) and score (vcfeval).
- `screen_read_homopolymers.sh`: the homopolymer and per-read quality screens.
- `sv.sh`, `sv_compare.py`: structural-variant calls and concordance.
- `compare_bedmethyl_stream.py`: CpG pileups compared site by site.
- `bootstrap.py`: paired block bootstrap of F1 changes.
- `snv_concordance.py`: genotype concordance between two VCFs without a truth set.
- `collect_benchmarks.py`: the tables in `benchmarks/`.

## Others

- `methyl/`: the methylation-option measurements in
  [docs/methylation-options.md](../docs/methylation-options.md).
- `usability/`: read mean-Q shifts in [docs/usability.md](../docs/usability.md).
- `package-release.sh`: packages a release build (used by the release workflow).
