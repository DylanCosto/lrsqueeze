# Usability checks

These checks measure what changes for downstream tools, beyond caller accuracy. Caller accuracy
and methylation are in [research-report.md](research-report.md). All data are HG002 GRCh38:

- ONT: chr20:10–40 Mb, Dorado HAC v5.2, 216,428 records.
- PacBio: Revio HiFi, 334,161 records.

Measured 2026-11-07 with samtools/htslib 1.24 and Clair3 2.0.3.

## Read mean Q after binning

Read-level filters and QC tools such as chopper, NanoFilt, and NanoPlot recompute each read's
mean Q from its per-base qualities. They use the mean error probability:
`-10 log10(mean(10^(-Q/10)))`.

Binning stores each bin as its median. A bin's mean error is dominated by its lowest scores, so
binned reads look slightly better. Method: 30,000 primary reads per platform,
`scripts/usability/meanq.py` (bins alone), and for the default as stored,
`scripts/usability/meanq_pair.py` (original against the archive).

| Bins | Shift in read mean Q: median (5th–95th percentile) | Q10 | Q15 | Q20 | Q25 |
|---|---|---:|---:|---:|---:|
| ONT default, 8 bins + homopolymer protection (measured on the archive) | +0.20 (−0.06 to +0.60) | 0.35% | 0.87% | 1.78% | 1.44% |
| ONT `--max`, 4 bins | +0.86 (+0.40 to +1.77) | 1.56% | 2.64% | 6.31% | 8.34% |
| PacBio default, 3 bins | +0.28 (−0.72 to +1.46) | 0% | 0.70% | 1.75% | 1.44% |

The cutoff columns give the share of reads that move across that cutoff.

The `--max` 4-bin figures ignore its homopolymer protection, which keeps 3.8% of ONT bases
native; protection would move them slightly toward the original.

Filters on Dorado's `qs` or PacBio's `rq` tag are unaffected, since `lrsqueeze` never changes those
tags. Storing each bin as its mean error instead of its median would remove most of the shift,
but those representatives have not been tested with variant callers.

## Readers

Clair3 2.0.3 was run on chr20:30–31 Mb with the same model and settings for each input:

| Input | Result |
|---|---|
| CRAM 3.0, written by samtools | 4,177 variants |
| CRAM 3.1 (default `normal` profile), written by samtools | "Failure to decode slice" |
| CRAM 3.1 normal and archive profiles, written by lrsqueeze | "Failure to decode slice" |
| CRAM 3.0, written by lrsqueeze (default mode) | 4,203 variants, no errors |

Clair3 2.0.3's decoder lacks the CRAM 3.1 codecs. This is why the default mode writes CRAM 3.0.
The variant counts differ because the inputs differ (original against binned qualities); this
check is about readability, not accuracy.

## Decoding speed

Full decode with `samtools view -u -o /dev/null`, single-threaded, on the same ONT records:

| Input | Size | Decode |
|---|---:|---:|
| Original BAM | 2125.1 MB | 8.8 s |
| lrsqueeze default, CRAM 3.0 (8 bins, before homopolymer protection) | 416.8 MB | 18.7 s |
| lrsqueeze `--max`, CRAM 3.1 archive | 288.5 MB | 20.8 s |
| lrsqueeze default, BAM | 1059.8 MB | ≈ 5 s (a 4-bin BAM measured 5.2 s) |

For PacBio, the original BAM decodes in 8.9 s and the default CRAM in 16.6 s. With 8 threads,
every format decodes in a few seconds.

Region queries cost about the same either way: a 1 Mb query takes 0.29 s on the BAM and 0.25 s
on the CRAM.
