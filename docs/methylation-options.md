# Methylation storage (`--methylation`)

On data called with all-context modification models, `MM`/`ML` dominate lrsqueeze's archives. On
the HG003 ONT genome (Dorado SUP v5.0, all-context 5mC/5hmC and 6mA), they make up about two
thirds of the default archive.

`--methylation` (or `[modifications] calls = ...` in a policy file) opts in to keeping fewer
calls. The default keeps them all.

| Mode | Keeps |
|---|---|
| `all` (default) | every call |
| `cpg` | 5mC/5hmC calls only at CpG sites; other modifications (e.g. 6mA) unchanged |
| `cpg-only` | 5mC/5hmC calls only at CpG sites, and nothing else |
| `none` | no `MM`/`ML`/`MN` |

A C call is kept at a CpG site: where the read has CpG, or where the base is aligned to a
reference CpG (the C of CG; for reverse-strand reads, the G). That is why the CpG modes need
`-T`.

Kept C parts become explicit (`?`). Calls the basecaller left out of an implicit (`.`) part mean
"unmodified", so at CpGs they become explicit 0-probability calls: CpG pileups keep their
coverage, and non-CpG Cs read as "no call", not "unmodified".

**What is lost:** non-CpG 5mC/5hmC and, with `cpg-only`, 6mA and any other modifications. Do not
use these modes for Fiber-seq or non-CpG methylation studies.

## Whole-genome result (ONT HG003)

The input is Dorado SUP v5.0, 110.1 GB BAM, all-context 5mC/5hmC and 6mA. Both archives use
lrsqueeze's default mode.

| | Archive | % of BAM | Encode (`-@ 16`) |
|---|---:|---:|---:|
| `--methylation all` (default) | 25.1 GB | 22.8% | 7 min |
| `--methylation cpg-only` | **10.9 GB** | **9.9%** | 5.5 min |

Genome-wide CpG pileups from the two archives were **byte-identical**: 28,675,315 sites,
5mC and 5hmC, `lrpipe methyl pileup --cpg --combine-strands --filter-threshold C:0.8`, the same
as modkit 0.6.4.

Records are rewritten on `-@` threads; when the rewrite ran on one thread, the `cpg-only`
encode took 19 min. The archive's records are identical either way.

## Slice measurements (HG003 chr20:10–40 Mb)

From the prototype, each variant then packed with 8 bins and 16-level `ML`:

| Variant | Archive | Change |
|---|---:|---:|
| all calls | 272.4 MB | — |
| drop 6mA only | 191.0 MB | −30% |
| `cpg` | 194.0 MB | −29% |
| `cpg-only` | 101.9 MB | −63% |
| `none` | 88.0 MB | −68% |

## Verification

The Rust implementation matches the Python prototype (`scripts/methyl/prototype.py`): `MM` and
`ML` are identical on all 4,301 records of HG003 chr20:20–22 Mb for both `cpg` and `cpg-only`
(`scripts/oracle/methyl_cpg.py`).

A 0.5 Mb check of the prototype against pysam's `modified_bases` found:

- 328,486 CpG calls kept with identical values
- 56,614 implicit CpG calls made explicit
- 6,218 calls kept only through the reference CpG rule
- 0 mismatches
