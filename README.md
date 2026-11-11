# lrsqueeze

Lossy archives for Oxford Nanopore and PacBio long reads

lrsqueeze takes an aligned BAM or CRAM (or a FASTQ), makes a few measured lossy changes to each
record, and writes standard CRAM, BAM or FASTQ.gz. Variant callers, SV callers, methylation tools
and read filters read the archive as they read the original, with no conversion step. It is one
Rust binary with no C dependencies, built on [pasta](https://github.com/DylanCosto/pasta).

## At a glance

Measured on four whole genomes against GIAB v4.2.1
([validation](docs/validation-wgs.md), [raw tables](benchmarks/README.md)):

| | Default | `--max` |
|---|---|---|
| Archive size | 20–26% of the BAM | ONT 20–22%, PacBio 3–5% |
| Small-variant F1 change | SNP within 0.003 pp; indel within 0.011 pp (ONT HAC −0.084 pp) | SNP within 0.014 pp; indel −0.072 to +0.10 pp |
| Structural variants | Identical VCFs (ONT, HiFi); 10 of 24,507 differ on Revio, from CRAM itself | Same |
| CpG methylation | Identical calls in pb-CpG-tools count mode; `--keep-ml` for fixed thresholds | Same |
| Format | CRAM 3.0: every htslib-based tool | CRAM 3.1: htslib 1.12 or newer |

## Quick start

```sh
lrsqueeze pack -T GRCh38.fa -o sample.cram sample.bam      # make the archive
lrsqueeze verify -T GRCh38.fa sample.bam sample.cram       # check every record before deleting the BAM
samtools view -T GRCh38.fa sample.cram | head               # read it like any CRAM
```

## Install

Tagged releases have prebuilt binaries for Linux (x86-64, glibc 2.28 or newer) and macOS (Apple
silicon) on the [releases page](https://github.com/DylanCosto/lrsqueeze/releases). To build from
source, with Rust 1.91 or newer and no C compiler:

```sh
cargo install --locked --git https://github.com/DylanCosto/lrsqueeze
```

or `cargo install --locked --path .` in a checkout.

## Usage

```sh
lrsqueeze pack -T GRCh38.fa -o sample.cram sample.bam          # default
lrsqueeze pack -T GRCh38.fa --max -o sample.cram sample.bam    # smaller files
lrsqueeze pack -T GRCh38.fa -o sample.lossy.bam sample.bam     # BAM output
lrsqueeze pack -T GRCh38.fa --methylation cpg-only -o sample.cram sample.bam   # opt-in, see below
lrsqueeze verify -T GRCh38.fa sample.bam sample.cram           # before deleting sample.bam
lrsqueeze pack -T GRCh38.fa -o sample.lossy.cram sample.cram   # CRAM input
```

**Run `lrsqueeze verify` before deleting an original.** It reads both files and checks every
record: names, flags, alignments, sequence, qualities and tags must be exactly what the policy
makes of the original. The only differences it allows are CRAM's own: `=`/`X` CIGAR operations
read back as `M`, and `MD`/`NM` are added on decoding. It exits with an error and lists the first
differing records otherwise. It also checks that the original's header lines (`@RG`, `@PG`,
`@CO`) are all still there. Both `pack` and `verify` warn when the original lacks its end-of-file
marker (BGZF's, or CRAM's end-of-file container): a copy cut short between two blocks reads as
complete, and its archive would verify.

On the 110 GB HG003 genome at `-@ 16`, verifying the CRAM took 8 minutes and 6.6 GB of memory,
against 5.5 minutes to pack it. BAM archives verify about as fast
as they are read.

An aligned CRAM (3.0 or 3.1, from samtools or any other writer) is read wherever a BAM is, by
`pack` and `verify`, decoded against `-T`. Its `MD` and `NM` tags stay as the CRAM stored them.
On the HG003 chr20 slice, an htslib CRAM packed to the same archive as its BAM, in 4.9 s against
4.4 s. On the whole HG003 genome at `-@ 16`, packing a CRAM took 7 minutes and 3.8 GB of memory,
and verifying against a CRAM original 8 minutes and 6.2 GB.

The platform comes from the header's `@RG PL` (`ONT` or `PACBIO`). Pass
`--platform ont|pacbio` if the header lacks it.

`-T` is the reference the reads are aligned to. CRAM needs it, and so do ONT's homopolymer
rules. PacBio BAM output is the only case that can be written without it.

Every output is indexed: `.crai`, or `.bai`/`.csi` for BAM. Every output records how it was made
in its header: an `@PG` line, a `@CO lrsqueeze lossy policy: ...` line in words, and a
`@CO lrsqueeze policy toml: ...` line that `verify` reads. `-@` sets the thread count, and `-` reads
BAM or CRAM from standard input.

lrsqueeze is for **aligned** reads. Unaligned reads (such as Dorado's default unaligned BAM) are
packed, with a warning: CRAM stores their bases in full, so on the HG003 chr20 slice the archive
was 59% larger. Their ONT qualities are protected as in [FASTQ](#fastq), in homopolymers found in
each read's own bases, which measured as accurate as the reference rule once the reads are
aligned. Align first, or keep unaligned reads as FASTQ.

To change the quality bins or anything else, start from the built-in policy as TOML:

```sh
lrsqueeze policy --platform ont > my-policy.toml   # add --max for the --max policy
$EDITOR my-policy.toml
lrsqueeze pack -T GRCh38.fa -c my-policy.toml -o sample.cram sample.bam
```

## FASTQ

```sh
lrsqueeze pack --platform ont -o reads.lossy.fastq.gz reads.fastq.gz
lrsqueeze verify --platform ont reads.fastq.gz reads.lossy.fastq.gz
```

Reads kept as FASTQ (plain or gzip, four-line records) are packed to FASTQ.gz with the same
quality bins and `ML` levels. With no alignment, the homopolymer rule looks at the read itself:
bases in a run of 5 or more identical bases keep their native quality. Names, comments and
sequences are kept byte for byte, except that the values of an `ML:B:C` field in the comment
(as `samtools fastq -T MM,ML` writes it) are quantized. The output is BGZF, which every gzip
reader opens.

On the HG003 chr20 reads (1.2 Gbases, 1.3 GB gzip) at `-@ 8`, packing took 7 s and verifying 6 s.
Plain gzip input is decompressed on one thread, which limits packing; BGZF input (from `bgzip`,
or an lrsqueeze archive) is decompressed on `-@` threads.

- A FASTQ records neither its platform nor its policy, so pass `--platform` to both `pack` and
  `verify`, and give `verify` the same `--max`, `--methylation`, `--keep-ml` or `-c` as `pack`.
- `--methylation` takes `all` or `none`; `cpg` and `cpg-only` need aligned reads.

Sizes on chr20 reads with `MM`/`ML`, against the original gzip:

| | Default | `--max` |
|---|---|---|
| ONT HG003 SUP (1.2 Gbases) | 702 MB of 1,293 MB (54%) | 593 MB (46%; 87 s to pack) |
| PacBio HG002 Revio (2.2 Gbases) | 844 MB of 968 MB (87%; Revio already uses 7 quality values) | 562 MB (58%; 53 s to pack) |

Aligners such as minimap2 ignore base qualities, so a packed FASTQ aligns exactly as the original
does: with `minimap2 -ax map-ont`, all 141,939 alignment records of the ONT chr20 reads were
identical before and after packing, apart from the qualities. Its BAM is therefore what packing
the aligned BAM with the read-based rule gives. Measured that way on ONT HAC chr1–6 with Clair3, the read-based rule changed F1 as the reference-based default
does: SNP −0.001 pp, indel −0.084 pp (−0.454 pp without protection). On ONT SUP HG003 it matched
the default within ±0.003 pp ([benchmarks](benchmarks/README.md#summary-fastq)).

## Changes

| | Default | `--max` |
|---|---|---|
| ONT qualities | 8 bins `0-4:4,5-9:7,10-14:12,15-19:17,20-24:23,25-29:28,30-34:32,35-93:40`; native qualities kept in reference homopolymers ≥ 5 bp | 4 bins `0-9:7,10-19:15,20-29:27,30-93:38`; native qualities kept in reference homopolymers ≥ 5 bp |
| PacBio qualities | 3 bins `0-17:10,18-29:25,30-93:40` (Revio already uses 7) | One quality per read, the mean of its Phred scores |
| Methylation (`ML`) | 16 levels (`--keep-ml` to keep exactly) | 16 levels (`--keep-ml` to keep exactly) |
| Methylation calls | All kept (`--methylation` to opt in to less) | All kept |
| Secondary alignments | Kept | Dropped |
| Container | CRAM 3.0 | CRAM 3.1 without the slow lzma and bzip2 codecs |

The homopolymer rule:

- Bases aligned to an exact A/C/G/T run of the reference keep their native quality, with no flanks.
- Inserted bases keep it too when the reference base on either side is in such a run.
- Soft clips are always binned.
- Unmapped reads keep native qualities in runs of 5 or more identical bases in the read itself, as
  FASTQ does. Archives made by 0.2.1 and earlier binned them, and their recorded policy says so.

**Opt-in: `--methylation cpg | cpg-only | none`.** On data called with all-context models,
`MM`/`ML` can be most of the archive.

- `cpg` keeps 5mC/5hmC calls only at CpG sites, where the read has CpG or the base is aligned to
  a reference CpG. Calls the basecaller left implicit become explicit 0-probability calls, so
  CpG pileups keep their coverage. Other modifications, such as 6mA, are kept.
- `cpg-only` also drops the other modifications.
- `none` removes the tags.

Non-CpG methylation (and with `cpg-only`, 6mA) is lost, so do not use these for Fiber-seq or
non-CpG methylation work. Measurements are in
[docs/methylation-options.md](docs/methylation-options.md). In a policy file, set
`[modifications] calls = "cpg-only"`.

**Opt-in: `--keep-ml`.** Keeps every `ML` probability exactly, instead of in 16 levels. Use it
when methylation will be called at a fixed threshold that is not a multiple of 1/16, such as
modkit's `--filter-threshold` or a threshold estimated on the original BAM. At such thresholds
the 16-level archives shift some sites by 20 pp or more (see
[docs/validation-wgs.md](docs/validation-wgs.md#methylation-cpg-genome-wide)). It costs about 7%
on PacBio Revio. In a policy file, remove `ml_levels`.

**Never changed (with the default methylation setting):** read names, flags, alignments, sequence, the `MM` tag, the read-level
quality tags `qs` and `rq`, `mv` move tables, and every other tag.

**One removal in both modes:** `MD`/`NM` on records without a sequence (`SEQ` = `*`, i.e.
secondaries). CRAM cannot regenerate these tags, so it would store them verbatim; every other
`MD`/`NM` is regenerated on decoding.

## Pipelines

Measured on whole genomes against GIAB v4.2.1:

- ONT HG003, Dorado SUP v5.0 (110 GB BAM)
- ONT HG002, HAC v5.0, chr1–6 (82 GB)
- PacBio HG002 Revio (45 GB, with `rq` and methylation tags)
- PacBio HG002 HiFi (39 GB)

Calls were made with Clair3, Sniffles2, and modkit/pb-CpG-tools equivalents (via lrpipe). Changes
are in F1, in percentage points. Details are in [docs/validation-wgs.md](docs/validation-wgs.md),
and the raw results in [benchmarks/](benchmarks/).

| | Default | `--max` |
|---|---|---|
| ONT small variants, SUP | SNP −0.003; indel −0.0004 (95% CI ±0.010), also in long homopolymers | SNP +0.014; indel −0.072 |
| ONT small variants, HAC | SNP −0.001; indel −0.084 (in GIAB's imperfect homopolymers) | SNP −0.002; indel +0.10 |
| PacBio small variants | SNP −0.001; indel −0.011 | SNP +0.0005; indel −0.005 (Revio; HiFi the same) |
| Structural variants | Identical VCF records (ONT); PacBio with `=`/`X` CIGARs: 10 of 24,507 differ, from CRAM itself | Same |
| Read mean Q recomputed from qualities (chopper, NanoFilt, NanoPlot) | ONT +0.20 Q (1.8% of reads cross a Q20 cutoff); PacBio +0.28 Q (1.75%) | ONT +0.86 Q (6.3% cross Q20); PacBio: every base carries the read's mean Phred score |
| Filters on `qs` or `rq` | Unchanged | Unchanged |
| CpG methylation | pb-CpG-tools count mode: identical calls at all 28.9 M sites (PacBio). modkit at a fixed 0.8: mean change 0.4 pp (ONT) and 1.0 pp (PacBio); at estimated thresholds, up to 29 pp at single sites. With `--keep-ml`: identical | Same |
| `flagstat`, read counts | Unchanged | Secondary alignments gone |
| Readers | Any htslib; Clair3 reads it directly | htslib ≥ 1.12 or pasta; Clair3 2.0.3 cannot, so restore to BAM first |

Both modes, when written as CRAM:

- Reading needs the same reference (`-T`/`--reference` for every tool).
- `=`/`X` CIGAR operations come back as `M`, as with any CRAM; the mismatches are recoverable from
  the sequence and the reference.
- Decoding is about 2× slower than the original BAM (19 s against 9 s single-threaded per 2 GB).

BAM output avoids these CRAM costs, at much larger sizes.

## Sizes

| Data | Output | Size | % of BAM | Time |
|---|---|---:|---:|---:|
| ONT HG003 whole genome (110.1 GB BAM) | default, CRAM 3.0 | 25.0 GB | 22.7% | 5.5 min |
| ONT HG003 whole genome | `--max`, CRAM 3.1 | 22.6 GB | 20.5% | 7 min |
| ONT HG003 whole genome | default + `--methylation cpg-only` | 10.9 GB | 9.9% | 5.5 min |
| PacBio HG002 whole genome (39.5 GB BAM) | default, CRAM 3.0 | 8.0 GB | 20.2% | 3 min |
| ONT HG002 HAC chr1–6 (82.4 GB BAM) | default, CRAM 3.0 | 21.8 GB | 26.4% | 5 min |
| ONT HG002 HAC chr1–6 | `--max`, CRAM 3.1 | 18.2 GB | 22.1% | 5.5 min |
| PacBio HG002 Revio whole genome (45.1 GB BAM) | default, CRAM 3.0 | 9.25 GB | 20.5% | 4.5 min |
| PacBio HG002 Revio whole genome | `--max`, CRAM 3.1 | 2.44 GB | 5.4% | 5 min |
| ONT HG002 chr20 (2125.1 MB BAM, `-@ 8`) | default, CRAM 3.0 | 465.3 MB | 21.9% | 17 s |
| ONT HG002 chr20 | `--max`, CRAM 3.1 | 313.4 MB | 14.7% | 46 s |
| ONT HG002 chr20 | default, BAM | 1133.7 MB | 53.3% | 9 s |
| PacBio Revio regions (2156.0 MB BAM, `-@ 8`) | default, CRAM 3.0 | 424.8 MB | 19.7% | 13 s |
| PacBio Revio regions | `--max`, CRAM 3.1 | 101.6 MB | 4.7% | 43 s |
| PacBio Revio regions | default, BAM | 1781.3 MB | 82.6% | 12 s |

Whole-genome encodes ran at `-@ 16`. Peak memory is 1–3 GB, except PacBio `--max` (10–13 GB at
`-@ 16`, scaling with threads). The HG003 default and Revio `--max` rows are from lrsqueeze 0.2.0
or later (pasta 0.9.0); the others are from 0.1.0, which packed CRAM 1.8–3 times slower on chr20 and 1.15 times slower on
the whole HG003 genome, with files up to 2.4% larger.

Use CRAM unless a tool cannot take it. With lossy qualities, BAM's size is dominated by the bases
themselves, which only CRAM stores as differences from the reference.

**Methylation tags can dominate.** On the HG003 genome (all-context 5mC/5hmC and 6mA calls),
`MM`/`ML` are about two thirds of the archive, so `--max` saves little there. On CpG-only or
6mA-free data, qualities dominate and `--max` saves more. Opting in to
`--methylation cpg-only` halves that genome's archive, with genome-wide CpG pileups
byte-identical.

The accuracy and methylation studies are in [docs/research-report.md](docs/research-report.md),
the read-filter and reader checks in [docs/usability.md](docs/usability.md), and the
whole-genome validation in [docs/validation-wgs.md](docs/validation-wgs.md).

## Verification

`scripts/oracle/` checks archives against independent Python implementations of the policies.
Both scripts decode through htslib (pysam):

- `verify.py` covers bins, one quality per read (`--qual mean-phred` for PacBio `--max`, or
  `rq` and `mean-error` for policies that name them), `ML`, and MD/NM.
- `ont_hybrid.py` covers ONT's homopolymer rule, with the projection of the study the ONT
  policies were validated in for mapped reads, and runs in each unmapped read's own bases
  (`--unmapped binned` for archives from 0.2.1 and earlier).

Every record matched for each mode and platform above, including PacBio `--max` (the mean Phred
score per read) on Revio and HiFi chr20 slices (6,125 and 6,010 records) and the ONT rule for
50,000 unmapped HG003 reads. The ONT `--max` quality bytes have the study's SHA-256
(`036d273941d61591d4396fc6c46e07e3664f5e6cdeb5e09ec89fc66cd0be7c42`), with `--unmapped binned`.

```sh
python -I scripts/oracle/ont_hybrid.py in.bam out.cram ref.fa --bins q8 --keep-secondary
python -I scripts/oracle/verify.py in.bam out.cram ref.fa --qual mean-phred --ml-levels 16 --strip-md-noseq --drop-secondary
```

## Notes

- **Index.** An index is written only when the records are in coordinate order: reference
  sequences in header order, positions ascending. Otherwise `lrsqueeze` warns and skips it, because
  region queries on such files return wrong results even with samtools' own index. An index left
  at the output's name by an earlier file is removed.
- **Calls are never dropped from `MM`/`ML`.** Once `ML` is quantized, dropping low-probability
  calls makes archives larger.
- **Methylation thresholds.** 16-level `ML` keeps calls identical at thresholds on multiples
  of 1/16 (0.5, 0.75, 0.8125). A fixed threshold inside a level, such as modkit's 0.8, moves
  the calls just above it, and one above the top level (248/256, 0.969) drops every call. Use
  `--keep-ml` when such thresholds are possible.
- **CRAM and `=`/`X` CIGARs.** pbmm2 writes `=`/`X`, which CRAM stores as `M`. Sniffles2 then
  calls about 0.04% of SVs differently; any CRAM, samtools' included, does the same.
- **6mA has not been evaluated** against a 6mA-positive (Fiber-seq) sample.
- **Earlier names.** lrsqueeze was called lrpack until 2026-11-08 and lrslim until 2026-11-11,
  and archives name the version that made them in their `@PG` and `@CO` lines. `verify` reads
  the policy that lrslim archives record; lrpack archives record none, so pass `--platform` (and
  `--max`) to check them.

## Documentation

- [docs/validation-wgs.md](docs/validation-wgs.md): the whole-genome validation (small variants,
  structural variants, methylation, sizes, a 19-genome clinical cohort).
- [benchmarks/](benchmarks/README.md): the raw tables behind it.
- [docs/usability.md](docs/usability.md): read filters, readers and tools.
- [docs/methylation-options.md](docs/methylation-options.md): the `--methylation` options,
  measured.
- [docs/research-report.md](docs/research-report.md): the study the presets come from.
- [scripts/](scripts/README.md): the validation scripts and the independent oracle.
- [CHANGELOG.md](CHANGELOG.md).

## License

MIT; see [LICENSE](LICENSE).
