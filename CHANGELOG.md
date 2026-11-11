# Changelog

Notable changes to lrsqueeze (called lrslim up to 0.2.1). The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## 0.9.0 - 2026-11-11

The first public release, as lrsqueeze. Versions jump from 0.2.1 to 0.9.0.

### Changed

- Renamed lrslim to lrsqueeze: the binary, the crate, the repository, and the `@PG` and `@CO`
  lines archives get (`@CO lrsqueeze policy toml: ...`). `verify` still reads the policy that
  lrslim archives record.
- `verify` decodes CRAM on `-@` threads (pasta's multithreaded CRAM reader, after 0.9.1): on
  HG003, 8 minutes instead of 16, with peak memory 6.6 GB instead of 4 GB. When the original
  is CRAM too, the two files share those threads.
- PacBio `--max`: each read gets the mean of its Phred scores, instead of a quality from `rq` (or
  from its mean error probability where there is no `rq`, as in FASTQ and older HiFi BAMs). Its
  calls now match the original's: indel F1 −0.005 pp on Revio HG002 and −0.003 pp on HiFi HG002,
  within noise, where `rq` cost 0.17 pp and the mean error 0.82 pp. Archives are slightly
  smaller. Policy files can name `source = "mean-phred"`; archives from earlier versions record
  `"rq"` and still verify.
- ONT: unmapped reads, such as those of Dorado's unaligned BAM, keep native qualities in
  homopolymers of 5 bp or more found in their own bases, as FASTQ does, instead of being binned
  throughout. The built-in policies record `homopolymers_in = "reference-or-read"`; archives
  made by 0.2.1 and earlier record no source, mean `"reference"`, and still verify.

### Added

- CRAM input: `pack` and `verify` read an aligned CRAM (3.0 or 3.1, such as one written by
  samtools) wherever they read a BAM, decoded against `-T`. `MD` and `NM` are kept as the CRAM
  stored them. Both warn when a CRAM lacks its end-of-file container, as they do for BGZF. CRAM
  is decoded on half the `-@` threads, alongside the other pools. On HG003 at `-@ 16`, packing a
  CRAM took 7 minutes, against 5.5 from the BAM.
- `--keep-ml`: keep base-modification probabilities (`ML`) exactly instead of in 16 levels, for
  methylation called at a fixed threshold that is not a multiple of 1/16. On eight PacBio Revio
  genomes, 16 levels moved 716–18,676 CpG sites by more than 10 pp at such a threshold;
  `--keep-ml` moved none, for an archive about 7% larger. It applies to FASTQ too, and
  `verify --platform` takes it for archives that do not record their policy.

### Fixed

- `pack` refuses an output that is its input (under any name, links included), which it would
  have overwritten while reading it. With CRAM input, `-o sample.cram sample.cram` is an easy
  slip.
- `pack` removes an index of an earlier file at the output's name (`.crai`, `.bai`, `.csi`) when
  it writes none, so region queries never read a stale index.
- `--methylation cpg` and `cpg-only` leave `MM`/`ML` as they are on records whose `MN` is not
  their read length, such as hard-clipped supplementary alignments that kept the primary's
  tags. Their calls were matched to the wrong bases, or stopped the run.
- Malformed records (shorter than BAM's fixed fields, or with a CIGAR that does not cover the
  read's bases where homopolymers are projected) are reported as errors instead of stopping a
  thread.
- CRAM and gzip are recognized on standard input even when the pipe delivers the first bytes in
  pieces.

## 0.2.1 - 2026-11-10

### Fixed

- Reads with more than 65,535 CIGAR operations, whose real CIGAR BAM keeps in the `CG` tag, were
  treated as fully soft-clipped: ONT homopolymer qualities were binned, and `--methylation cpg`
  and `cpg-only` saw no reference CpGs. `lrslim verify` now reports such reads in archives made
  by 0.1.0 and 0.2.0.
- pasta 0.9.1: 0.2.0 wrote such reads to CRAM with their CIGAR and also the `CG` tag that held
  it, so `lrslim verify` rejected the archive and htslib showed an extra tag.

## 0.2.0 - 2026-11-10

### Added

- FASTQ: `lrslim pack --platform ont|pacbio` packs plain or gzip FASTQ into a lossy FASTQ.gz with
  the same quality bins and `ML` levels. Native qualities are kept in homopolymers of 5 bp or
  more found in each read's own bases, which measured the same as the reference-based rule on
  ONT HAC (indel F1 −0.084 pp either way). `--max` writes BGZF level 10, which on chr20 reads was
  2.6 times as fast as level 12 for files 2.4% larger.
- `lrslim verify --platform` (with `--max` and `--methylation`) for archives that do not record
  their policy, such as FASTQ.gz.
- Policy files: `homopolymers_in = "read"` applies the homopolymer rule to the read's bases
  instead of the reference's.
- `pack` and `verify` warn when a BAM or bgzipped input lacks BGZF's end-of-file marker, as a copy
  cut short between two blocks does: such a file reads as complete.

### Changed

- pasta 0.9.0, which writes CRAM as HTSlib does: on the HG003 ONT and HG002 Revio chr20 reads at
  `-@ 16`, packing to CRAM is 2.8–3.0 times as fast (1.8 times with `--max`), and archives are
  0.3–2.4% smaller. Their records are unchanged. The whole HG003 genome, where the rest of the
  pipeline sets the pace, packs in 5.5 minutes instead of 6.3, to 25.0 GB instead of 25.1.
- `lrslim verify` refuses a policy given with `-c` or `--platform` that differs from the one the
  archive records, instead of checking against it.

## 0.1.0 - 2026-11-09

The first release.

### Added

- `lrslim pack`: lossy CRAM 3.0 (or BAM) archives of aligned ONT and PacBio BAMs, 20–26% of the
  BAM's size, with changes measured as negligible for variant callers, structural-variant
  callers and CpG methylation pileups on four whole genomes (see `benchmarks/`).
  - ONT: 8 quality bins, keeping native qualities in reference homopolymers of 5 bp or more.
  - PacBio: 3 quality bins.
  - Base-modification probabilities (`ML`) quantized to 16 levels.
  - `MD`/`NM` removed from records without a sequence, which CRAM cannot regenerate.
- `--max`: smaller files in CRAM 3.1 (ONT 4 quality bins, PacBio one quality per read from `rq`,
  secondary alignments dropped), at a small measured accuracy cost.
- `--methylation cpg | cpg-only | none`: opt in to keeping 5mC/5hmC calls only at CpG sites, or
  dropping modification calls.
- `lrslim verify`: checks every record of an archive against what its policy makes of the
  original, before the original is deleted. Archives record their policy in the header for it.
- `lrslim policy` and `pack --config`: the built-in policies as TOML, to edit.
- A warning when most records are unmapped, since unaligned reads get no homopolymer protection
  and no reference-based compression.
- Indexes (`.crai`, `.bai`/`.csi`) for coordinate-sorted output; records are transformed on `-@`
  threads.
