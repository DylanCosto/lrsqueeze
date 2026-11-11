# Lossy long-read archive format: ONT and PacBio

Completed 7 November 2026, before lrsqueeze existed: this is the study its presets come from. It builds on an earlier ONT quality-binning study (October 2026), which supplies the ONT quality-binning accuracy results. No new container format is needed. Every recipe here is a set of record transforms followed by standard **CRAM 3.1 (archive profile)**, so archives remain readable by samtools/htslib 1.24. (lrsqueeze's default since 0.2.0 writes CRAM 3.0 instead; see [validation-wgs.md](validation-wgs.md) for the whole-genome results of what ships.)

## Recommendation

| Platform | Recipe | Archive vs source BAM | vs lossless CRAM 3.1 | Measured accuracy cost |
|---|---|---:|---:|---|
| **ONT** (R10.4.1, Dorado HAC v5.2) | hybrid-q4 quals + ML 16 levels + no secondary alignments | **296.9 MB / 2125.1 MB = 14.0%** (7.2× smaller) | −75.9% vs earlier lossless CRAM (1234.5 MB) | SNP F1 −0.002 pp; indel F1 −0.22 pp (95% CI −0.49…+0.04); 5mC CpG mean abs diff 0.24 pp |
| **PacBio** (Revio HiFi) | 3-level quals (Q10/Q25/Q40) + ML 16 levels | **357.0 MB / 2156.0 MB = 16.6%** | −44.2% vs lossless CRAM (639.6 MB) | none measurable (indel +2 TP, SNP −1 FP) |
| **PacBio**, aggressive | one constant Q per read from `rq` + ML 16 levels | **96.8 MB / 2156.0 MB = 4.5%** (22× smaller) | −84.9% | SNP F1 −0.02 pp; indel F1 −0.19 pp (13 more FN of 4054) |

Two ONT figures are exact rather than estimated. Restoring the recommended ONT archive to BAM and re-running Clair3 gave **byte-identical VCF records** (87,025) to the earlier hybrid4 study. Every record in `ont.hybrid4.ml16.cram` also decodes to exactly the expected transform of the source (216,428 records, `results/verify.hybrid4.ml16.json`).

## Where the bytes go

CRAM 3.1 archive, compressed MB per stream group (`results/streams.json`):

| Recipe | quality | ML | MM | MD/NM | seq-diff | names | other | **total** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ONT lossless (MD not forced) | 929.6 | 135.4 | 13.4 | 23.2 | 74.1 | 3.1 | 12.7 | **1191.5** |
| ONT hybrid-q4 | 172.3 | 135.4 | 13.4 | 23.2 | 74.1 | 3.1 | 12.7 | **434.2** |
| ONT hybrid-q4 + ML16 | 172.3 | 36.9 | 13.3 | 23.2 | 75.0 | 3.1 | 12.7 | **336.5** |
| ONT + strip MD on SEQ-less records | 172.3 | 36.9 | 13.3 | 0.0 | 75.2 | 3.1 | 12.7 | **313.5** |
| **ONT + drop secondary alignments** | 167.2 | 36.9 | 13.3 | 0.0 | 67.4 | 2.1 | 9.9 | **296.7** |
| ONT no quals at all (floor) | 0 | 26.8 | 29.0 | 23.2 | 76.4 | 3.0 | 12.6 | 171.0 |
| PacBio lossless | 520.8 | 45.0 | 10.4 | 0 | 43.2 | 1.0 | 19.0 | **639.4** |
| PacBio 3-level + ML16 | 261.7 | 22.1 | 10.4 | 0 | 42.6 | 1.0 | 19.0 | **356.8** |
| PacBio per-read Q + ML16 | 0.8 | 22.1 | 10.4 | 0 | 43.3 | 1.0 | 19.0 | **96.6** |

Findings:

1. **Quality scores dominate both platforms.** They are 78% of a lossless ONT CRAM and 81% of a PacBio one. Revio already bins to 7 levels, and those quals are still 521 of 639 MB.
2. **After binning, ONT's methylation ML tag becomes the next target.** It was the largest stream (135 MB) because this run used `--modified-bases-threshold 0.0` and stored a 6mA probability for every A. Quantizing ML to 16 levels cuts it 73%.
3. **Dropping low-probability calls is counter-productive once ML is quantized.** Removing implicit `.`-mode calls below p=0.05 (Dorado-style) shrinks ML but doubles MM (13→29 MB), because runs of `0,` deltas become varied gaps. Quantize-only (1092.9 MB) beats drop+quantize (1101.1 MB).
4. **The earlier archives wasted space on MD/NM in two ways.** First, `store_md=1,store_nm=1` stored tags CRAM can regenerate (−43 MB). Second, even with defaults, 78,565 secondary alignments (36% of records) carry `SEQ=*`, so CRAM must store their MD verbatim (23 MB). These secondaries are ignored by Clair3, Sniffles and modkit.
5. **Codec tuning has little left.** `level=9` saves another 1.4% (296.9→292.9 MB). Very large slices made encoding over 20× slower (aborted after more than 25 minutes vs 80 s) and are not worth it.

## Accuracy evidence

**ONT small variants** come from the earlier HG002 study on the same BAM: Clair3 2.0.3 HAC v5.2, 21.8 Mb GIAB v4.2.1 confident, 1-Mb block bootstrap. See `research/ont-quality-binning-2026-10-22/native-giab/followup/results/REPORT.md`.

- Hybrid-q4: SNP F1 0.997815 vs 0.997835; indel F1 0.880052 vs 0.882236 (−0.218 pp). That clears a −0.5 pp non-inferiority margin overall (one-sided lower bound −0.447 pp).
- **Inside long homopolymers it does not clear the margin** (−0.41 pp, lower bound −0.91 pp).
- Median-q8 is the alternative if homopolymer indels matter most: −0.19 pp overall, with a stronger bound. It costs about 64 MB more on this sample.
- ML quantization, MD handling and dropping secondaries do not change Clair3 input. This was verified as identical calls, so these accuracy numbers transfer exactly.

**PacBio small variants** (new): Clair3 `hifi_revio`, HG002 chr20:12.1–40 Mb, same truth, confident regions and RTG settings.

| Quals | SNP TP/FP/FN | SNP F1 | Indel TP/FP/FN | Indel F1 |
|---|---|---:|---|---:|
| Revio native 7-level | 25601/42/32 | 0.998557 | 3998/40/56 | 0.988136 |
| 3-level 10/25/40 | 25601/41/32 | 0.998576 | 4000/41/54 | 0.988264 |
| per-read constant from `rq` | 25598/49/35 | 0.998362 | 3985/42/69 | 0.986264 |

**Methylation (ONT 5mC/5hmC at CpG, combined strands)**: modkit 0.6.4, fixed `--filter-threshold C:0.8`, 291,581 sites with coverage ≥5 in chr20:12.2–39.9 Mb, compared against original ML.

| ML levels | ML stream | 5mC Pearson | mean abs diff | p95 abs diff | ≥50% call agreement |
|---|---:|---:|---:|---:|---:|
| 8 | 25.0 MB | 0.99931 | 0.61 pp | 2.79 pp | 99.50% |
| **16** | **36.9 MB** | **0.99978** | **0.24 pp** | **1.85 pp** | **99.77%** |
| 32 | 52.2 MB | 0.99986 | 0.17 pp | 1.35 pp | 99.86% |
| original | 135.4 MB | 1 | 0 | 0 | 100% |

For scale, binomial sampling noise at a 64%-methylated CpG at 30× is about ±9 pp (1 SD). With modkit's automatic threshold, quantization also nudges the estimated threshold (0.799→0.811). Pin `--filter-threshold` when comparing archived and original data.

## Caveats

- **Scope.** One ONT HG002 sample (chr20:10–40 Mb) and one PacBio HG002 sample (chr20). This is not genome-wide or clinical validation; [validation-wgs.md](validation-wgs.md) has the whole-genome results. The HG003 whole-genome run in the earlier study did not complete (next point).
- **CRAM 3.1 archive codecs are not universally readable.** The Clair3 2.0.3 install here fails on these archives with "Failure to decode slice". That caused the earlier WGS failure, and samtools 1.24 decodes the same files fine. Either restore to BAM before calling, as done here, or make sure every consumer links htslib with CRAM 3.1 codec support.
- **6mA could not be evaluated.** This ONT run looks like a no-methyltransferase control (mean 6mA about 0.02%). Fiber-seq or m6A footprinting data needs its own check before quantizing `A+a`.
- **Lossy items you should choose deliberately.** Secondary alignments are dropped; keep `--strip-md-noseq` (313.5 MB) instead if any tool needs them. Move tables (`mv`) are now used by Clair3 v2 `*_with_mv` models; this sample had none, but do not strip `mv` without checking. PacBio kinetics were already absent from this BAM. Kinetics can make HiFi BAMs up to about 5× larger per PacBio, so strip them unless you will re-call methylation or run Fiber-seq.
- The PacBio per-read constant was tested with Clair3 only. DeepVariant, assembly (hifiasm ignores QV anyway) and SV callers were not tested.
- POD5 or raw signal retention is out of scope.

## Literature context

Sources were checked in November 2026.

- [HTSlib's CRAM benchmark](https://www.htslib.org/benchmarks/CRAM.html) shows the same breakdown (ONT: qualities 50%, MM+ML 44%; Revio: qualities 80%). It notes that ML "could benefit from a similar quantisation" but none exists. We found no published study of lossy ML.
- [CoLoRd](https://github.com/refresh-bio/colord) is the closest lossy long-read compressor (ONT about 6.9× with 4 quality bins), but it is FASTQ-only: no alignments, tags or methylation.
- [Rivara-Espasandín et al. 2022](https://doi.org/10.1093/bioadv/vbac054) reported ONT 4-bin quals as near-neutral for SNPs, with small indel losses (R9, PEPPER-DeepVariant). That is consistent with the HG002 results here.

## Reproduce

The study ran its record transforms in Python, then wrote CRAM 3.1 (archive profile) with samtools and measured streams with `samtools cram-size`. Those transforms carried over into this repository as the independent oracle, `scripts/oracle/transform.py` (`--qual-lut`, `--qual`, `--ml-levels`, `--ml-drop`, `--strip-md-noseq`, `--drop-secondary`), with `scripts/oracle/verify.py` for the record-level decode check. Hybrid-q4 keeps native qualities on reference homopolymers of 5 bp or more and uses median-q4 (7/15/27/38) elsewhere; its projection is `protection()` in `scripts/oracle/ont_hybrid.py`. lrsqueeze now applies the same transforms itself, in Rust, and its archives match the oracle record for record.

Calls used Clair3 2.0.3 with RTG vcfeval against GIAB v4.2.1, and methylation used modkit 0.6.4 pileups compared site by site.
