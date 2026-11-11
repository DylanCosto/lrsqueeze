# Whole-genome validation

The validation checks whether variant calling on lrsqueeze archives gives the same answers as
calling on the original BAMs. It was measured on 2026-11-07 and 2026-11-08 with lrsqueeze at
pasta `3cea7a3`. The scripts are in `scripts/wgs/`:

- `validate.sh` and `validate_dataset.sh`: encode, call, score.
- `bootstrap.py`: paired block bootstrap.
- `sv.sh` and `sv_compare.py`: SV concordance.
- `compare_bedmethyl_stream.py`: methylation pileups.
- `snv_concordance.py`: genotype concordance between two lrpipe VCFs (clinical cohort).
- `collect_benchmarks.py`: the tables in [benchmarks/](../benchmarks/).

## Data

| Dataset | Basecalling and modifications | Original BAM | Scored regions |
|---|---|---:|---|
| ONT HG003 SUP | Dorado SUP v5.0; all-context 5mC/5hmC and 6mA ([ONT GIAB 2025.01](https://epi2me.nanoporetech.com/giab-2025.01/), PAY87954) | 110.1 GB | confident autosomes without chr20 (2.47 Gb; chr20 was used in the pilots) |
| ONT HG002 HAC | Dorado HAC v5.0; all-context 5mC/5hmC (same release, PAW71238); chr1–6 only | 82.4 GB | confident chr1–6 |
| PacBio HG002 Revio | Revio rep1, aligned by PacBio (pbmm2); `rq`, 5mC MM/ML, 7-level qualities | 45.1 GB | confident autosomes |
| PacBio HG002 HiFi | aligned with bonefish; Revio qualities, no `rq`/MM/ML | 39.5 GB | confident autosomes |

The truth is GIAB v4.2.1 for HG002 and HG003. Each dataset uses the reference it was aligned
to: ONT's release GRCh38 for ONT, and the GRCh38 no-alt analysis set for PacBio.

## Method

Small variants were called with lrpipe 0.1.8 `snv`, a Clair3 2.0.3 port, on the CPU in exact
mode, with the same settings for every input. The models were `r1041_e82_400bps_sup_v500`,
`r1041_e82_400bps_hac_v500`, and `hifi_revio`. lrpipe read the CRAM archives directly.

Scoring used RTG vcfeval 3.13 on PASS calls, at the threshold-free summary row. Confidence
intervals come from a paired bootstrap over 10 Mb genome blocks (20,000 draws): TP, FP, and FN
counts are pooled per draw, and F1 is recomputed from them. The long-homopolymer stratum is GIAB
v3.3's `AllHomopolymers_ge7bp_imperfectge11bp_slop5`, intersected with the confident regions.

## Small variants: change in F1 against the original

The tables give percentage points with 95% intervals; an interval containing 0 is no
detectable change.

### Default

| Dataset | SNP | Indel | Indel, long homopolymers |
|---|---:|---:|---:|
| ONT HG003 SUP | −0.003 [−0.004, −0.002] | **−0.0004 [−0.010, +0.010]** | −0.006 [−0.027, +0.016] |
| ONT HG002 HAC (chr1–6) | −0.001 [−0.003, +0.000] | −0.084 [−0.106, −0.061] | −0.189 [−0.237, −0.139] |
| PacBio HG002 Revio | −0.001 [−0.002, −0.001] | −0.011 [−0.017, −0.005] | |
| PacBio HG002 HiFi (bonefish) | −0.0006 [−0.0012, +0.0001] | −0.011 [−0.016, −0.005] | |

### `--max`

| Dataset | SNP | Indel |
|---|---:|---:|
| ONT HG003 SUP | +0.014 [+0.012, +0.016] | −0.072 [−0.094, −0.051] |
| ONT HG002 HAC (chr1–6) | −0.002 [−0.003, +0.000] | **+0.102 [+0.071, +0.133]** |
| PacBio HG002 Revio | +0.0005 [−0.0006, +0.0016] | −0.005 [−0.014, +0.005] |
| PacBio HG002 HiFi | −0.0002 [−0.0014, +0.0011] | −0.003 [−0.013, +0.007] |

### For reference: ONT HG003 with 8 bins and no homopolymer protection

| SNP | Indel | Indel, long homopolymers | Extra indel FP |
|---:|---:|---:|---:|
| −0.003 | −0.156 [−0.177, −0.137] | −0.320 [−0.365, −0.276] | +1,522 |

The default with protection gave −78 extra indel FP instead.

### Original F1

| Dataset | SNP | Indel |
|---|---:|---:|
| ONT HG003 SUP | 0.996910 | 0.891892 |
| ONT HG002 HAC | 0.997991 | 0.851182 |
| PacBio Revio | 0.998463 | 0.991564 |
| PacBio HiFi | 0.998288 | 0.992202 |

### Findings

- **Binning costs ONT indels almost only in long homopolymers.** Without protection, that
  stratum gained 1,529 false positives on HG003, while the rest of the genome lost 7.
- **Protection fixes SUP fully.** Keeping native qualities on exact reference A/C/G/T runs of
  5 bp or more makes the SUP default indistinguishable from the original.
- **HAC keeps a small loss.** It sits in GIAB's imperfect homopolymers and their ±5 bp flanks,
  which the exact-run rule does not cover: 249 of the 257 lost TP indels fall there. Wider
  protection does not fix it (next section).
- **`--max` moves indels in either direction.** 4 bins with protection improved HAC indels (392
  fewer FP) and cost SUP indels slightly. PacBio `--max` (one quality per read, its mean
  Phred score) matches the original; see "One quality per read" below.

### Wider homopolymer protection on HAC (not adopted)

Variants of the default's rule, each encoded and called genome-wide on HAC chr1–6. Native
qualities were kept on exact reference runs of at least *length* bp, plus *flank* bp on each side.
F1 changes are against the original:

| Length, flank | Indel | Indel, long homopolymers | Extra indel FP | Archive |
|---|---:|---:|---:|---:|
| 5, 0 (default) | −0.084 | −0.189 | +113 | 21.75 GB |
| 4, 0 | −0.074 | −0.166 | +90 | 23.76 GB (+9%) |
| 5, 2 | −0.161 | −0.340 | +471 | 22.45 GB (+3%) |
| 5, 5 | −0.131 | −0.272 | +410 | 23.41 GB (+8%) |
| 4, 2 | −0.142 | −0.300 | +414 | 25.57 GB (+18%) |

SNPs changed by −0.001 pp in every variant.

- **Flanks hurt.** Every flanked variant added about 300 indel false positives over the default,
  nearly all in long homopolymers. Paired against the default, flank 5 is −0.048 pp
  [−0.067, −0.028]. A plausible cause, not tested, is that the caller reads the step from native to
  binned qualities next to a run as signal.
- **4 bp runs gain little.** +0.009 pp [+0.002, +0.016] over the default, about a tenth of its
  loss, for 9% more space.

The default stays at 5 bp runs with no flank. Calls restricted to the stratum (Clair3
`--bed_fn`, even padded by 100 bp) found about 3,000 fewer true indels there than genome-wide
calls, most likely because phasing then sees fewer SNPs, so these screens must call genome-wide.

### Homopolymers found in the read, for FASTQ (adopted for FASTQ)

A FASTQ has no alignment, so `lrsqueeze` keeps native qualities in homopolymers found in each
read's own bases (`homopolymers_in = "read"`). Aligners such as minimap2 ignore base qualities,
so a packed FASTQ aligns exactly as the original does, and its calls are those of the aligned
BAM packed with the read-based rule. That is how the rule was measured: each variant was encoded
from the aligned BAM and called genome-wide. F1 changes against the original:

| Rule | SNP | Indel | Indel, long homopolymers | Extra indel FP | Archive |
|---|---:|---:|---:|---:|---:|
| HAC: reference runs ≥ 5 bp (default) | −0.001 | −0.084 | −0.189 | +113 | 21.75 GB |
| HAC: read runs ≥ 5 bp | −0.001 | −0.084 [−0.106, −0.063] | −0.191 [−0.240, −0.142] | +118 | 21.79 GB |
| HAC: read runs ≥ 4 bp | −0.001 | −0.075 [−0.097, −0.053] | −0.167 [−0.216, −0.117] | +93 | 23.85 GB |
| HAC: no protection | −0.001 | −0.454 [−0.491, −0.417] | −0.913 [−0.989, −0.836] | +1,720 | 20.15 GB |
| SUP: read runs ≥ 5 bp | −0.003 [−0.004, −0.002] | −0.002 [−0.012, +0.009] | −0.006 [−0.027, +0.016] | −74 | 25.27 GB |

SUP is scored on all confident autosomes here (the default's headline uses `primary`).

- **Read runs match reference runs.** Paired against the default, read runs of 5 bp change indel
  F1 by −0.0007 pp [−0.008, +0.007] on HAC and −0.0003 pp [−0.003, +0.003] on SUP; the archives
  are within 0.5% of each other in size.
- **4 bp runs gain little**, as with reference runs: +0.010 pp over 5 bp, for 9.5% more space.
- **Protection is what keeps binning safe.** Without it, HAC indel F1 drops 5 times as much,
  from about 1,600 extra false positives.

FASTQ uses read runs of 5 bp, the same length as the reference rule.

The premise was checked directly: the HG003 SUP chr20 reads (72,135 reads, 1.2 Gbases) were
aligned with minimap2 2.31 (`-ax map-ont --eqx`) before and after packing, and all 141,939
alignment records (primary, secondary and supplementary) were identical in every field but
`QUAL`. PacBio FASTQ (HG002 Revio chr20, 141,297 reads) packs and verifies in both modes; its
default is the BAM policy.

### One quality per read

PacBio `--max` gives each read one quality. It first came from `rq`, which a FASTQ, and some BAMs,
do not have. Three sources were measured, encoded from the aligned BAM and called genome-wide
(F1 changes against the original):

| One quality per read from | SNP | Indel | Archive |
|---|---:|---:|---:|
| Revio: `rq` (0.1.0–0.2.1) | −0.001 | −0.166 | 2.51 GB |
| Revio: mean error probability | −0.020 [−0.024, −0.016] | −0.827 [−0.877, −0.778] | 2.45 GB |
| **Revio: mean Phred score (adopted)** | +0.0005 [−0.0006, +0.0016] | −0.005 [−0.014, +0.005] | 2.44 GB |
| HiFi, no `rq`: mean error probability (0.1.0–0.2.1) | −0.021 [−0.028, −0.015] | −0.821 [−0.870, −0.773] | 1.13 GB |
| **HiFi, no `rq`: mean Phred score (adopted)** | −0.0002 [−0.0014, +0.0011] | −0.003 [−0.013, +0.007] | 1.13 GB |

The mean error probability is dominated by a read's few low-quality bases: on 15,436 Revio reads
its median was Q30, against Q34 from `rq` and Q39 for the mean Phred score. The mean Phred score
is close to what the default's 3 bins give most bases (Q30 and above become Q40), and matched the
original's calls on both genomes. Reads keep their `rq` tags either way.

## Structural variants

Calls were made with lrpipe 0.1.8 `sv`, a Sniffles2 2.8.1 port.

| Dataset | Default | `--max` |
|---|---|---|
| ONT HG003 SUP | identical VCF records (27,046) | identical (27,046) |
| ONT HG002 HAC | identical (9,430) | identical (9,430) |
| PacBio HG002 HiFi (bonefish) | identical (26,287) | |
| PacBio HG002 Revio | 24,497 of 24,507 PASS SVs identical; 10 lost | the same 10 lost |

The 10 Revio SVs are lost by CRAM itself, not by lrsqueeze. A lossless samtools CRAM of the same
region loses the same calls.

Revio's pbmm2 alignments use `=`/`X` CIGAR operations, e.g. `7I20X2=1X…` beside indels. CRAM
stores those as `M`, and Sniffles2's indel clustering uses the `X` runs. All 10 calls had
4–8 supporting reads.

## Methylation (CpG, genome-wide)

The default archive (16-level `ML`) was compared with the original BAM.

**modkit-style pileup.** `lrpipe methyl pileup --cpg --combine-strands --filter-threshold C:0.8`,
which is byte-identical to modkit 0.6.4, at sites with coverage ≥ 5:

| | Sites | Pearson | Mean change | p95 change | ≥ 50% calls unchanged | Valid coverage |
|---|---:|---:|---:|---:|---:|---:|
| ONT HG003, 5mC | 28.56 M | 0.99936 | 0.39 pp | 2.8 pp | 99.57% | +0.7% |
| ONT HG003, 5hmC | 28.56 M | 0.977 | 0.05 pp | 0.0 pp | 99.9995% | |
| PacBio Revio, 5mC | 28.03 M | 0.99758 | 0.98 pp | 5.2 pp | 99.20% | −2.8% |

**pb-CpG-tools count mode.** `lrpipe methyl cpg --pileup-mode count`, which matches
pb-CpG-tools 3.0.0, on PacBio Revio. The methylation %, coverage, and methylated and unmethylated
counts were **byte-identical at all 28,934,010 CpG sites**. Only the average-probability columns
changed, by 0.005 on average.

Why the two methods differ: 16 levels split probabilities into bins 1/16 wide, each stored as
its midpoint. A threshold on a bin boundary gives identical calls. pb-CpG-tools counts a call as
methylated at ≥ 0.5, which is a boundary. A threshold inside a bin moves the calls just above it
below it. modkit's fixed 0.8 lies inside the bin from 0.75 to 0.81, which is stored as 0.78, and
that cost PacBio 2.8% of valid coverage.

Prefer thresholds at multiples of 1/16 when analysing these archives with a fixed cutoff
(0.5, 0.75, 0.8125); 0.8 is not one. modkit's automatic per-sample threshold re-estimates from
the archive itself.

**The tail, at thresholds estimated on the original.** Eight further Revio genomes (clinical,
not public) were packed with the default policy. Each was compared with its original using one
threshold estimated on the original's chr1 (0.631–0.643, the 10th percentile of call
confidences, as modkit does), with `--cpg --combine-strands`, at sites with coverage ≥ 10:

| | Range over the 8 genomes |
|---|---:|
| Sites with any change in valid coverage | 2.3–6.7 M of 28.8 M |
| p99 change | 3.5–6.1 pp |
| Largest change | 19.5–28.6 pp |
| Sites changed by more than 10 pp | 716–18,676 |
| Calls filtered as low-confidence | 3–12% fewer |

The mean change is small, but site-level work at such a threshold (differential methylation,
imprinting, single-locus checks) sees thousands of sites move by more than 10 pp.

**Thresholds above the top level.** The largest value 16 levels can store is 248/256 (0.97). On
two lower-coverage Revio genomes from the same cohort, the threshold estimated on the original was
0.998, which no archive probability reaches: about 10 M of 11.5 M sites lost all their calls, and
the p99 change at the remaining sites was 96 pp. Any fixed threshold above 0.97 does this. Use
`--keep-ml` when such thresholds are possible.

**`--keep-ml`** keeps `ML` exactly. One of the two genomes with the largest shift above was
packed again with it:

| | Default | `--keep-ml` |
|---|---:|---:|
| Archive (of a 37.6 GB BAM) | 7.28 GB (19.4%) | 7.82 GB (20.8%) |
| CpG sites with any change | 6.67 M of 28.8 M | **0** |
| Sites changed by more than 10 pp | 18,278 | 0 |

Base qualities are binned the same way, so the variant results above apply to both.

## Clinical cohort: the whole pipeline on original and archive

Measured on 2026-11-09 and 2026-11-10 with lrsqueeze 0.1.0 (then called lrslim; default policy) and
lrpipe 0.1.3. Nineteen clinical PacBio Revio genomes (not public) were packed: 723 GB of BAM
became 115 GB of archive (16%). Each archive passed `lrsqueeze verify`, `samtools quickcheck`, and a
separate record-by-record comparison with the original.

On ten of them, lrpipe ran the same analyses (`snv`, `phase`, `sv`, `cnv`, `tr`, `depth`,
`methyl`) on the original and on the archive, and `lrpipe compare` compared the outputs. There is
no truth set for these genomes, so this measures agreement, not accuracy.

**Small variants.** PASS, non-reference calls; genotypes are compared unordered, because which
haplotype a phase block labels 1 is arbitrary (`0|1` and `1|0` are the same genotype).
`snv_concordance.py` streams both VCFs and uses about 23 MB. About 4.25 M SNPs and 0.95 M indels
per genome:

| | Range over the 10 genomes |
|---|---:|
| SNP concordance (same genotype / union) | 0.99740–0.99803 |
| SNPs at a shared site with a different genotype | 669–1,045 |
| SNPs called only on the original / only on the archive | 4,133–5,605 / 3,483–4,539 |
| Indel concordance | 0.97947–0.99393 |
| Indels called only on the original / only on the archive | 2,805–7,579 / 3,045–12,577 |

The calls that differ are borderline. In one genome, 99.6% of the SNPs called only on the
original have QUAL < 10, against 7% of the calls both share.

**Everything else.** In all ten genomes, the SV VCF, the tandem-repeat genotypes and the depth
windows and summary were identical. CNV calls differed only in the random IDs of LOH records and
in one or two LOH boundaries per genome, which follow low-quality heterozygous calls inside long
homozygous stretches.

**Operational notes.**

- **Keep the reference with the archives.** A CRAM decodes only against the reference it was
  written with. The FASTA the BAM was aligned to (and its `.fai`) should sit next to the archives,
  with its checksum, before the originals are deleted.
- **A wrong reference can go unnoticed.** When the input BAM's `@SQ` lines have no `M5`, the
  writer takes the reference it is given, records that reference's checksums and encodes against
  it. Nothing fails, but the policy's reference-based decisions (homopolymer protection, CpG
  context) are made at the wrong bases, and `lrsqueeze verify` given the same reference agrees.
  Readers reject the file against the right reference (checksum mismatch).
- **Sex detection from an index.** lrpipe's automatic sex call reads per-contig counts from the
  index. A `.crai` is coarser than a `.bai`, and on two lower-coverage archives it fit neither
  XX nor XY. Pass `--sex` when running lrpipe on such archives.

## Sizes and times

The whole-genome encodes ran at `-@ 16`; calls at 14 threads.

| Dataset | Default | `--max` |
|---|---:|---:|
| ONT HG003 SUP (110.1 GB) | 25.1 GB (22.8%), 7 min, 1.6 GB peak | 22.6 GB (20.5%), 7 min, 3.2 GB peak |
| ONT HG002 HAC chr1–6 (82.4 GB) | 21.8 GB (26.4%), 5 min, 1.6 GB peak | 18.2 GB (22.1%), 5.5 min, 2.8 GB peak |
| PacBio HG002 Revio (45.1 GB) | 9.25 GB (20.5%), 4.5 min, 1.7 GB peak | **2.44 GB (5.4%)**, 4.7 min, 12.6 GB peak |
| PacBio HG002 HiFi (39.5 GB) | 8.0 GB (20.2%), 3.3 min, 1.3 GB peak | 1.13 GB (2.9%), 3.3 min, 9.9 GB peak |

PacBio `--max` keeps the archive profile's 50 Mbase slices, because its qualities are per read.
Its peak memory grows with the thread count (10–13 GB at `-@ 16`).

The earlier `--max` with lzma and bzip2 took 62 min on HG003, for a 21.5 GB archive.

On the all-context ONT data, `MM`/`ML` make up about two thirds of the archive. That is why
`--max` is only 10–17% smaller than the default there. See
[methylation-options.md](methylation-options.md) for CpG-only storage.
