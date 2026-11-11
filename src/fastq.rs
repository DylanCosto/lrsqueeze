//! Lossy FASTQ, for reads without an alignment.
//!
//! The policy's quality bins and `ML` levels apply as they do to BAM records, but with nothing to
//! align against, protected homopolymers are found in each read's own bases
//! ([`Policy::for_fastq`]). Records must have four lines. Everything except the qualities, and
//! the values of an `ML:B:C,...` field in the comment (as `samtools fastq -T MM,ML` writes it), is
//! kept byte for byte. Output is BGZF, which every gzip reader opens.

use std::{
    fs::File,
    io::{self, BufRead, BufReader, Read, Write},
    mem,
    num::NonZero,
    path::Path,
    sync::mpsc::{self, Receiver, SyncSender},
    thread,
};

use anyhow::{Context, bail};
use flate2::read::MultiGzDecoder;
use pasta::bgzf::{self, io::writer::CompressionLevel};

use crate::{
    alignments::peek,
    homopolymer::Runs,
    pack::{BATCH_BYTES, IndexOutcome, Input, PackOptions, Report},
    pipeline,
    policy::{Mode, ModificationCalls, Platform, Policy},
    quality::MAX_PHRED,
    transform::{QualityOp, Stats, bin_qualities, fill_quality, ml_lut, read_quality},
    verify::{VerifyOptions, VerifyReport},
};

/// FASTQ qualities are Phred scores plus this.
const OFFSET: u8 = 33;

/// The BGZF level of the output: pasta's fast default, or with `--max` 10. On chr20 reads at
/// `-@ 8`, 10 was 2.6 times (ONT) and 1.8 times (Revio) as fast as the highest level, 12, for
/// files 2.4% and 1% larger; 9 and below were 6–7% larger again.
fn level(mode: Mode) -> CompressionLevel {
    match mode {
        Mode::Safe => CompressionLevel::default(),
        Mode::Max => CompressionLevel::new(10).expect("10 is a DEFLATE level"),
    }
}

/// Whether `path` names a FASTQ file: `.fastq`, `.fq`, or either with `.gz` (any case).
pub fn is_fastq_path(path: &Path) -> bool {
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    name.ends_with(".fastq") || name.ends_with(".fq")
}

/// One four-line record, without the line ends.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct Record {
    /// The header line after `@`: the name, then any comment.
    header: Vec<u8>,
    sequence: Vec<u8>,
    /// The separator line after `+` (usually empty).
    separator: Vec<u8>,
    /// Phred+33 qualities.
    quality: Vec<u8>,
}

impl Record {
    fn name(&self) -> String {
        let end = self
            .header
            .iter()
            .position(|b| matches!(b, b' ' | b'\t'))
            .unwrap_or(self.header.len());
        String::from_utf8_lossy(&self.header[..end]).into_owned()
    }

    fn write(&self, out: &mut impl Write) -> io::Result<()> {
        out.write_all(b"@")?;
        out.write_all(&self.header)?;
        out.write_all(b"\n")?;
        out.write_all(&self.sequence)?;
        out.write_all(b"\n+")?;
        out.write_all(&self.separator)?;
        out.write_all(b"\n")?;
        out.write_all(&self.quality)?;
        out.write_all(b"\n")
    }
}

/// Reads four-line FASTQ records, plain or gzip-compressed (detected from the content).
pub(crate) struct Reader {
    inner: Box<dyn BufRead + Send>,
    records: u64,
}

impl Reader {
    /// A reader of `source`, the file at `path` if given. BGZF (as `bgzip` and lrsqueeze write it)
    /// is decompressed on `threads` threads with pasta; other gzip on one thread of its own.
    pub(crate) fn new(
        source: Box<dyn Read + Send>,
        path: Option<&Path>,
        threads: NonZero<usize>,
    ) -> io::Result<Self> {
        // Enough for a BGZF header (is_bgzf).
        let (head, source) = peek(source, 16)?;
        let buffered = BufReader::with_capacity(1 << 20, source);
        let inner: Box<dyn BufRead + Send> = if is_bgzf(&head) {
            Box::new(BufReader::with_capacity(
                1 << 20,
                Bgzf {
                    blocks: Some(bgzf::io::MultithreadedReader::with_worker_count(
                        threads, buffered,
                    )),
                    path: path.map(Path::to_path_buf),
                    delivered: 0,
                    fallback: None,
                },
            ))
        } else if head.starts_with(&GZIP_MAGIC) {
            Box::new(BufReader::with_capacity(
                1 << 20,
                Decompressor::new(buffered),
            ))
        } else {
            Box::new(buffered)
        };
        Ok(Self { inner, records: 0 })
    }

    fn open(input: &Input, threads: NonZero<usize>) -> anyhow::Result<Self> {
        let (source, path): (Box<dyn Read + Send>, _) = match input {
            Input::Path(path) => (
                Box::new(File::open(path).with_context(|| format!("opening {}", path.display()))?),
                Some(path.as_path()),
            ),
            Input::Stdin => (Box::new(io::stdin()), None),
        };
        Self::new(source, path, threads).context("reading the FASTQ input")
    }

    /// Reads the next record into `record`; `false` at the end of the input.
    pub(crate) fn read_record(&mut self, record: &mut Record) -> anyhow::Result<bool> {
        let n = self.records + 1;
        let mut line = |buf: &mut Vec<u8>, what: &str| -> anyhow::Result<bool> {
            buf.clear();
            let read = self
                .inner
                .read_until(b'\n', buf)
                .with_context(|| format!("reading FASTQ record {n}"))?;
            if buf.ends_with(b"\n") {
                buf.pop();
                if buf.ends_with(b"\r") {
                    buf.pop();
                }
            }
            if read == 0 && what != "header" {
                bail!("FASTQ record {n} ends before its {what} line");
            }
            Ok(read > 0)
        };
        let mut header = std::mem::take(&mut record.header);
        if !line(&mut header, "header")? {
            record.header = header;
            return Ok(false);
        }
        let Some(rest) = header.strip_prefix(b"@") else {
            bail!("FASTQ record {n} does not start with '@' (only four-line FASTQ is read)");
        };
        record.header = rest.to_vec();
        line(&mut record.sequence, "sequence")?;
        let mut separator = std::mem::take(&mut record.separator);
        line(&mut separator, "separator")?;
        let Some(rest) = separator.strip_prefix(b"+") else {
            bail!("FASTQ record {n} has no '+' line (only four-line FASTQ is read)");
        };
        record.separator = rest.to_vec();
        line(&mut record.quality, "quality")?;
        if record.quality.len() != record.sequence.len() {
            bail!(
                "FASTQ record {n} ({}) has {} bases but {} qualities",
                record.name(),
                record.sequence.len(),
                record.quality.len()
            );
        }
        self.records = n;
        Ok(true)
    }
}

const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// Whether `head`, the start of a file, is a BGZF block: gzip with the extra field's first
/// subfield `BC` holding the block size (SAM specification, section 4.1). A start too short to
/// tell is read as plain gzip, which BGZF also is.
pub(crate) fn is_bgzf(head: &[u8]) -> bool {
    const FEXTRA: u8 = 1 << 2;
    head.len() >= 16
        && head.starts_with(&GZIP_MAGIC)
        && head[3] & FEXTRA != 0
        && u16::from_le_bytes([head[10], head[11]]) >= 6
        && head[12..16] == [b'B', b'C', 2, 0]
}

/// BGZF, decompressed in parallel. A gzip member that is not BGZF (a BGZF file with gzip
/// appended, say) stops pasta's reader; a file is then read again from the start as plain gzip,
/// skipping what was already read.
struct Bgzf {
    blocks: Option<bgzf::io::MultithreadedReader<BufReader<Box<dyn Read + Send>>>>,
    path: Option<std::path::PathBuf>,
    /// Decompressed bytes read so far.
    delivered: u64,
    fallback: Option<Decompressor>,
}

impl Read for Bgzf {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if let Some(fallback) = &mut self.fallback {
            return fallback.read(buf);
        }
        let blocks = self
            .blocks
            .as_mut()
            .expect("BGZF is read until the fallback");
        match blocks.read(buf) {
            Ok(n) => {
                self.delivered += n as u64;
                Ok(n)
            }
            Err(e) => {
                let Some(path) = &self.path else {
                    return Err(io::Error::new(
                        e.kind(),
                        format!(
                            "{e} (a stream that starts as BGZF must be BGZF throughout; read it \
                             from a file instead)"
                        ),
                    ));
                };
                self.blocks = None;
                let mut fallback = Decompressor::new(File::open(path)?);
                let skipped = io::copy(&mut (&mut fallback).take(self.delivered), &mut io::sink())?;
                if skipped < self.delivered {
                    return Err(e);
                }
                self.fallback.insert(fallback).read(buf)
            }
        }
    }
}

/// Decompresses gzip on its own thread, so that decompression and parsing overlap.
struct Decompressor {
    /// Decompressed chunks; `None` marks the end of the input.
    chunks: Receiver<io::Result<Option<Vec<u8>>>>,
    /// Chunks that have been read, for the thread to fill again.
    spent: SyncSender<Vec<u8>>,
    chunk: Vec<u8>,
    position: usize,
    finished: bool,
}

impl Decompressor {
    const CHUNK: usize = 1 << 20;
    /// Chunks filled ahead of the reader.
    const AHEAD: usize = 4;

    fn new(source: impl Read + Send + 'static) -> Self {
        let (sender, chunks) = mpsc::sync_channel(Self::AHEAD);
        let (spent, recycled) = mpsc::sync_channel::<Vec<u8>>(Self::AHEAD + 1);
        thread::spawn(move || {
            let mut decoder = MultiGzDecoder::new(source);
            loop {
                let chunk = Self::fill(&mut decoder, recycled.try_recv().unwrap_or_default());
                let last = !matches!(chunk, Ok(Some(_)));
                // A failed send means the reader was dropped.
                if sender.send(chunk).is_err() || last {
                    return;
                }
            }
        });
        Self {
            chunks,
            spent,
            chunk: Vec::new(),
            position: 0,
            finished: false,
        }
    }

    /// Reads the next chunk into `chunk`, or `None` at the end of the input.
    fn fill(decoder: &mut impl Read, mut chunk: Vec<u8>) -> io::Result<Option<Vec<u8>>> {
        chunk.resize(Self::CHUNK, 0);
        let mut filled = 0;
        while filled < chunk.len() {
            match decoder.read(&mut chunk[filled..]) {
                Ok(0) => break,
                Ok(n) => filled += n,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
        chunk.truncate(filled);
        Ok((filled > 0).then_some(chunk))
    }
}

impl Read for Decompressor {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.finished {
            return Ok(0);
        }
        if self.position == self.chunk.len() {
            match self.chunks.recv() {
                Ok(Ok(Some(chunk))) => {
                    let spent = mem::replace(&mut self.chunk, chunk);
                    self.position = 0;
                    // The thread allocates a new chunk if this one is not taken back.
                    let _ = self.spent.try_send(spent);
                }
                Ok(Ok(None)) => {
                    self.finished = true;
                    return Ok(0);
                }
                Ok(Err(e)) => return Err(e),
                // Without the end marker, the thread did not finish the input.
                Err(_) => return Err(io::Error::other("the gzip decompression thread stopped")),
            }
        }
        let n = buf.len().min(self.chunk.len() - self.position);
        buf[..n].copy_from_slice(&self.chunk[self.position..self.position + n]);
        self.position += n;
        Ok(n)
    }
}

/// Applies a policy to FASTQ records.
#[derive(Clone, Debug)]
pub(crate) struct Transformer {
    quality: QualityOp,
    ml_lut: Option<Box<[u8; 256]>>,
    drop_modifications: bool,
    keep: Vec<bool>,
    stats: Stats,
}

impl Transformer {
    /// A transformer for `policy` as it applies to FASTQ ([`Policy::for_fastq`]).
    ///
    /// # Errors
    ///
    /// Fails for what FASTQ cannot do: dropping qualities, or keeping methylation calls at CpG
    /// sites, which needs the reference.
    pub(crate) fn new(policy: &Policy) -> anyhow::Result<Self> {
        let quality = QualityOp::new(&policy.quality);
        if let QualityOp::Drop = quality {
            bail!("FASTQ records need qualities; the policy drops them");
        }
        let drop_modifications = match policy.modifications.calls {
            ModificationCalls::All => false,
            ModificationCalls::None => true,
            ModificationCalls::Cpg | ModificationCalls::CpgOnly => bail!(
                "keeping methylation calls at CpG sites needs aligned reads and the reference; \
                 FASTQ supports --methylation all or none"
            ),
        };
        Ok(Self {
            quality,
            ml_lut: policy.modifications.ml_levels.map(|n| Box::new(ml_lut(n))),
            drop_modifications,
            keep: Vec::new(),
            stats: Stats::default(),
        })
    }

    pub(crate) fn stats(&self) -> Stats {
        self.stats
    }

    /// Rewrites the record's qualities and the `ML`/`MM`/`MN` fields of its comment.
    pub(crate) fn apply(&mut self, record: &mut Record) -> anyhow::Result<()> {
        self.stats.records_in += 1;
        self.stats.records_out += 1;
        self.rewrite_comment(record)?;
        self.stats.quality_bases += record.quality.len() as u64;
        if let QualityOp::Keep = self.quality {
            return Ok(());
        }
        if let Some(&c) = record
            .quality
            .iter()
            .find(|&&c| !(OFFSET..=OFFSET + MAX_PHRED).contains(&c))
        {
            bail!(
                "read {} has the quality character {:?}, outside Phred+33 ('!' to '~')",
                record.name(),
                char::from(c)
            );
        }
        // Phred+33 characters to Phred scores, and back once rewritten.
        let quals = &mut record.quality;
        quals.iter_mut().for_each(|q| *q -= OFFSET);
        match &self.quality {
            QualityOp::Bins { lut, preserve } => {
                self.keep.clear();
                self.keep.resize(quals.len(), false);
                // Every protected homopolymer is in the read: there is no reference.
                if let Some((n, _)) = preserve {
                    Runs::find(&record.sequence, *n).mark(&mut self.keep);
                }
                bin_qualities(quals, &self.keep, lut, &mut self.stats);
            }
            QualityOp::PerRead { source, max } => {
                let q = read_quality(*source, quals, None, *max);
                fill_quality(quals, q, &mut self.stats);
            }
            QualityOp::Keep | QualityOp::Drop => unreachable!("handled above"),
        }
        quals.iter_mut().for_each(|q| *q += OFFSET);
        Ok(())
    }

    /// Quantizes `ML:B:C` values, or removes `MM`, `ML` and `MN` fields, in the comment. Fields
    /// are separated by tabs or spaces; the first is the read name.
    fn rewrite_comment(&mut self, record: &mut Record) -> anyhow::Result<()> {
        let header = &record.header;
        let needs_work = |(_, text): (&[u8], &[u8])| {
            (self.drop_modifications && is_modification_field(text))
                || (self.ml_lut.is_some() && text.starts_with(b"ML:B:C,"))
        };
        if !split_with_separators(header).skip(1).any(needs_work) {
            return Ok(());
        }
        let mut out = Vec::with_capacity(header.len());
        let mut dropped = false;
        for (i, (separator, text)) in split_with_separators(header).enumerate() {
            if i > 0 && self.drop_modifications && is_modification_field(text) {
                dropped = true;
                continue;
            }
            out.extend_from_slice(separator);
            match (&self.ml_lut, text.strip_prefix(b"ML:B:C,")) {
                (Some(lut), Some(values)) if i > 0 => {
                    out.extend_from_slice(b"ML:B:C,");
                    for (j, value) in values.split(|&b| b == b',').enumerate() {
                        let v: u8 = std::str::from_utf8(value)
                            .ok()
                            .and_then(|v| v.parse().ok())
                            .with_context(|| {
                                format!("an ML value of read {} is not 0-255", record.name())
                            })?;
                        if j > 0 {
                            out.push(b',');
                        }
                        push_decimal(&mut out, lut[usize::from(v)]);
                        self.stats.ml_values += 1;
                    }
                }
                _ => out.extend_from_slice(text),
            }
        }
        self.stats.modifications_rewritten += u64::from(dropped);
        record.header = out;
        Ok(())
    }
}

/// Appends `v` in decimal.
fn push_decimal(out: &mut Vec<u8>, v: u8) {
    if v >= 100 {
        out.push(b'0' + v / 100);
    }
    if v >= 10 {
        out.push(b'0' + v / 10 % 10);
    }
    out.push(b'0' + v % 10);
}

fn is_modification_field(field: &[u8]) -> bool {
    field.starts_with(b"MM:") || field.starts_with(b"ML:") || field.starts_with(b"MN:")
}

/// The header's fields, each with the separator before it (empty for the first).
fn split_with_separators(header: &[u8]) -> impl Iterator<Item = (&[u8], &[u8])> {
    let mut rest = header;
    let mut first = true;
    std::iter::from_fn(move || {
        if rest.is_empty() && !first {
            return None;
        }
        let separator_len = if first {
            0
        } else {
            rest.iter()
                .position(|b| !matches!(b, b' ' | b'\t'))
                .unwrap_or(rest.len())
        };
        first = false;
        let (separator, after) = rest.split_at(separator_len);
        let text_len = after
            .iter()
            .position(|b| matches!(b, b' ' | b'\t'))
            .unwrap_or(after.len());
        let (text, remaining) = after.split_at(text_len);
        rest = remaining;
        Some((separator, text))
    })
}

/// The policy a FASTQ is packed with: the custom one, or the platform's built-in one, with the
/// command line's overrides, as it applies to reads without an alignment.
fn resolve_policy(options: &PackOptions) -> anyhow::Result<(Policy, Option<Platform>)> {
    let (policy, platform) = match (&options.custom_policy, options.platform) {
        (Some(policy), _) => (policy.clone(), None),
        (None, Some(platform)) => (Policy::preset(platform, options.mode), Some(platform)),
        (None, None) => bail!(
            "a FASTQ does not say which platform made it; pass --platform ont or --platform pacbio"
        ),
    };
    let policy = policy.with_overrides(options.overrides()).for_fastq();
    policy.validate().map_err(anyhow::Error::msg)?;
    Ok((policy, platform))
}

/// Packs a FASTQ (plain or gzip) into a lossy FASTQ.gz.
pub(crate) fn pack(options: &PackOptions) -> anyhow::Result<Report> {
    if let Input::Path(path) = &options.input
        && !is_fastq_path(path)
    {
        bail!(
            "FASTQ output needs FASTQ input (.fastq, .fq, or either with .gz), not {}",
            path.display()
        );
    }
    if !options
        .output
        .to_string_lossy()
        .to_ascii_lowercase()
        .ends_with(".gz")
    {
        bail!("FASTQ output is gzip-compressed; name it .fastq.gz or .fq.gz");
    }
    let (policy, platform) = resolve_policy(options)?;
    let transformer = Transformer::new(&policy)?;
    let mut reader = Reader::open(&options.input, options.threads)?;
    let file = File::create(&options.output)
        .with_context(|| format!("creating {}", options.output.display()))?;
    let mut writer = bgzf::io::multithreaded_writer::Builder::default()
        .set_compression_level(level(options.mode))
        .set_worker_count(options.threads)
        .build_from_writer(file);
    let workers = pipeline::ordered(
        |send| read_batches(&mut reader, send),
        vec![transformer; options.threads.get()],
        transform_batch,
        |text: Vec<u8>| writer.write_all(&text).context("writing FASTQ records"),
    )?;
    writer.finish().context("finishing the output")?;
    Ok(Report {
        index: IndexOutcome::NotRequested,
        policy,
        platform,
        stats: workers.iter().fold(Stats::default(), |mut stats, worker| {
            stats += worker.stats();
            stats
        }),
    })
}

/// Reads records into batches of about [`BATCH_BYTES`] and sends them until the input ends or the
/// workers stop.
fn read_batches(
    reader: &mut Reader,
    send: &mut dyn FnMut(Vec<Record>) -> bool,
) -> anyhow::Result<()> {
    let mut batch = Vec::new();
    let mut bytes = 0;
    let mut record = Record::default();
    while reader.read_record(&mut record)? {
        bytes += record.header.len() + 2 * record.sequence.len();
        batch.push(mem::take(&mut record));
        if bytes >= BATCH_BYTES {
            if !send(mem::take(&mut batch)) {
                // The workers stopped; the pipeline reports why.
                return Ok(());
            }
            bytes = 0;
        }
    }
    if !batch.is_empty() {
        send(batch);
    }
    Ok(())
}

/// Transforms a batch of records into FASTQ text.
fn transform_batch(transformer: &mut Transformer, batch: Vec<Record>) -> anyhow::Result<Vec<u8>> {
    let mut text = Vec::with_capacity(
        batch
            .iter()
            .map(|r| r.header.len() + r.sequence.len() + r.separator.len() + r.quality.len() + 6)
            .sum(),
    );
    for mut record in batch {
        transformer.apply(&mut record)?;
        record.write(&mut text).context("writing a FASTQ record")?;
    }
    Ok(text)
}

/// Compares a FASTQ.gz archive with what its policy makes of the original FASTQ. A FASTQ cannot
/// record its policy, so `options.policy` must hold it.
pub(crate) fn verify(options: &VerifyOptions) -> anyhow::Result<VerifyReport> {
    let policy = options.policy.as_ref().context(
        "a FASTQ does not record its policy; pass the --platform (and --max or --methylation) it \
         was packed with, or the policy file",
    )?;
    let policy = policy.clone().for_fastq();
    policy.validate().map_err(anyhow::Error::msg)?;
    let transformer = Transformer::new(&policy)?;
    let mut original = Reader::open(&Input::Path(options.original.clone()), options.threads)?;
    let mut archive = Reader::open(&Input::Path(options.archive.clone()), options.threads)?;
    let mut report = VerifyReport::new(policy, Vec::new());
    pipeline::ordered(
        |send| read_pairs(&mut original, &mut archive, send),
        vec![transformer; options.threads.get()],
        compare_batch,
        |checked: Vec<Checked>| {
            for (in_original, difference) in checked {
                let index = report.records;
                if in_original {
                    report.records += 1;
                }
                if let Some((name, what)) = difference {
                    report.record(index, name, what);
                }
            }
            Ok(())
        },
    )?;
    Ok(report)
}

/// A record of the original with the archive's record in the same place, or a record that only
/// one of them has.
enum Pair {
    Both(Record, Record),
    Missing(Record),
    Extra(Record),
}

/// Whether a compared record is in the original, and the name and difference if it differs.
type Checked = (bool, Option<(String, &'static str)>);

/// Reads the two files' records in step into batches of about [`BATCH_BYTES`] and sends them
/// until both end or the workers stop.
fn read_pairs(
    original: &mut Reader,
    archive: &mut Reader,
    send: &mut dyn FnMut(Vec<Pair>) -> bool,
) -> anyhow::Result<()> {
    let mut batch = Vec::new();
    let mut bytes = 0;
    loop {
        let (mut expected, mut actual) = (Record::default(), Record::default());
        let pair = match (
            original.read_record(&mut expected)?,
            archive.read_record(&mut actual)?,
        ) {
            (true, true) => Pair::Both(expected, actual),
            (true, false) => Pair::Missing(expected),
            (false, true) => Pair::Extra(actual),
            (false, false) => break,
        };
        bytes += match &pair {
            Pair::Both(r, _) | Pair::Missing(r) | Pair::Extra(r) => {
                r.header.len() + 2 * r.sequence.len()
            }
        };
        batch.push(pair);
        if bytes >= BATCH_BYTES {
            if !send(mem::take(&mut batch)) {
                // The workers stopped; the pipeline reports why.
                return Ok(());
            }
            bytes = 0;
        }
    }
    if !batch.is_empty() {
        send(batch);
    }
    Ok(())
}

fn compare_batch(transformer: &mut Transformer, batch: Vec<Pair>) -> anyhow::Result<Vec<Checked>> {
    batch
        .into_iter()
        .map(|pair| {
            Ok(match pair {
                Pair::Both(mut expected, actual) => {
                    transformer.apply(&mut expected)?;
                    let what = [
                        (expected.header != actual.header, "header line"),
                        (expected.sequence != actual.sequence, "sequence"),
                        (expected.separator != actual.separator, "separator line"),
                        (expected.quality != actual.quality, "quality scores"),
                    ]
                    .into_iter()
                    .find_map(|(differs, what)| differs.then_some(what));
                    (true, what.map(|what| (expected.name(), what)))
                }
                Pair::Missing(expected) => {
                    (true, Some((expected.name(), "missing from the archive")))
                }
                Pair::Extra(actual) => (false, Some((actual.name(), "extra in the archive"))),
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn open_bytes(bytes: Vec<u8>) -> io::Result<Reader> {
        Reader::new(
            Box::new(io::Cursor::new(bytes)),
            None,
            NonZero::new(2).expect("2 is nonzero"),
        )
    }

    fn records(text: &str) -> Vec<Record> {
        let mut reader = open_bytes(text.as_bytes().to_vec()).expect("readable");
        let mut out = Vec::new();
        let mut record = Record::default();
        while reader.read_record(&mut record).expect("valid") {
            out.push(record.clone());
        }
        out
    }

    #[test]
    fn reads_four_line_records_and_gzip() {
        let text = "@r1 c1\nACGT\n+\n!!II\n@r2\nAC\n+r2\nII\n";
        let plain = records(text);
        assert_eq!(plain.len(), 2);
        assert_eq!(
            (plain[0].name(), &plain[1].separator[..]),
            ("r1".into(), &b"r2"[..])
        );
        let gzip = |text: &str| {
            let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            gz.write_all(text.as_bytes()).expect("compress");
            gz.finish().expect("compress")
        };
        // Two gzip members, each over a decompression chunk.
        let many = text.repeat(150_000);
        let compressed = [gzip(&many), gzip(&many)].concat();
        let mut reader = open_bytes(compressed.clone()).expect("readable");
        let mut record = Record::default();
        for i in 0..600_000 {
            assert!(reader.read_record(&mut record).expect("valid"));
            assert_eq!(record, plain[i % 2]);
        }
        assert!(!reader.read_record(&mut record).expect("the end"));
        assert!(!reader.read_record(&mut record).expect("still the end"));

        let truncated = compressed[..compressed.len() - 100].to_vec();
        let mut reader = open_bytes(truncated).expect("readable");
        let error = loop {
            match reader.read_record(&mut record) {
                Ok(true) => {}
                Ok(false) => panic!("a truncated gzip reads to the end"),
                Err(error) => break error,
            }
        };
        assert!(
            format!("{error:#}").contains("incomplete deflate stream"),
            "{error:#}"
        );
    }

    #[test]
    fn reads_bgzf_in_parallel() -> anyhow::Result<()> {
        let text = "@r1 c1\nACGT\n+\n!!II\n@r2\nAC\n+r2\nII\n".repeat(100_000);
        let mut writer = bgzf::io::Writer::new(Vec::new());
        writer.write_all(text.as_bytes())?;
        let compressed = writer.finish()?;
        assert!(is_bgzf(&compressed));
        let mut reader = open_bytes(compressed.clone())?;
        let mut record = Record::default();
        for i in 0..200_000 {
            assert!(reader.read_record(&mut record)?);
            assert_eq!(record.name(), if i % 2 == 0 { "r1" } else { "r2" });
        }
        assert!(!reader.read_record(&mut record)?);

        // Cut inside a data block, past the end-of-file marker.
        let truncated = compressed[..compressed.len() - 100].to_vec();
        let mut reader = open_bytes(truncated)?;
        loop {
            match reader.read_record(&mut record) {
                Ok(true) => {}
                Ok(false) => panic!("a truncated BGZF reads to the end"),
                Err(_) => break,
            }
        }

        // Plain gzip, and a start too short to tell, are not BGZF.
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        gz.write_all(text.as_bytes())?;
        assert!(!is_bgzf(&gz.finish()?) && !is_bgzf(&compressed[..15]));
        Ok(())
    }

    #[test]
    fn reads_bgzf_followed_by_plain_gzip_from_a_file() -> anyhow::Result<()> {
        let (a, b) = ("@r1\nACGT\n+\nIIII\n", "@r2\nAC\n+\nII\n");
        let mut writer = bgzf::io::Writer::new(Vec::new());
        writer.write_all(a.as_bytes())?;
        let mut bytes = writer.finish()?;
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        gz.write_all(b.as_bytes())?;
        bytes.extend(gz.finish()?);
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("mixed.fastq.gz");
        std::fs::write(&path, &bytes)?;

        let mut reader = Reader::open(&Input::Path(path), NonZero::new(2).expect("nonzero"))?;
        let mut record = Record::default();
        let mut names = Vec::new();
        while reader.read_record(&mut record)? {
            names.push(record.name());
        }
        assert_eq!(names, ["r1", "r2"]);

        // A stream cannot be read again, so it fails, saying why.
        let mut reader = open_bytes(bytes)?;
        assert!(reader.read_record(&mut record)?);
        let error = reader
            .read_record(&mut record)
            .expect_err("not BGZF throughout");
        assert!(
            format!("{error:#}").contains("BGZF throughout"),
            "{error:#}"
        );
        Ok(())
    }

    #[test]
    fn rejects_quality_characters_outside_phred_33() {
        let mut t = Transformer::new(&Policy::preset(Platform::Ont, Mode::Safe).for_fastq())
            .expect("valid");
        for bad in ["@r\nACGT\n+\nII I\n", "@r\nACGT\n+\nII\x7fI\n"] {
            let mut record = records(bad).remove(0);
            let error = t.apply(&mut record).expect_err("outside Phred+33");
            assert!(error.to_string().contains("outside Phred+33"), "{error}");
        }
    }

    #[test]
    fn writes_decimals() {
        let mut out = Vec::new();
        for v in [0, 7, 10, 99, 100, 255] {
            push_decimal(&mut out, v);
            out.push(b',');
        }
        assert_eq!(out, b"0,7,10,99,100,255,");
    }

    #[test]
    fn rejects_mismatched_lengths_and_multiline_records() {
        let mut reader = open_bytes(b"@r\nACGT\n+\nII\n".to_vec()).expect("ok");
        assert!(reader.read_record(&mut Record::default()).is_err());
        let mut reader = open_bytes(b"r\nACGT\n+\nIIII\n".to_vec()).expect("ok");
        assert!(reader.read_record(&mut Record::default()).is_err());
    }

    #[test]
    fn bins_qualities_outside_read_homopolymers_and_quantizes_ml() {
        let policy = Policy::preset(Platform::Ont, Mode::Safe).for_fastq();
        let mut t = Transformer::new(&policy).expect("FASTQ policy");
        // Q13 ('.') everywhere; AAAAA is a 5 bp run in the read and keeps Q13, the rest goes to
        // the 10-14 bin (12, '-').
        let mut record =
            records("@r\tMM:Z:C+m?,0;\tML:B:C,0,100,255\nCAAAAAG\n+\n.......\n").remove(0);
        t.apply(&mut record).expect("valid");
        assert_eq!(record.quality, b"-.....-");
        assert_eq!(record.header, b"r\tMM:Z:C+m?,0;\tML:B:C,8,104,248");
        assert_eq!(t.stats().ml_values, 3);
    }

    #[test]
    fn drops_modification_fields_and_keeps_the_rest() {
        let mut policy = Policy::preset(Platform::Pacbio, Mode::Safe);
        policy.modifications.calls = ModificationCalls::None;
        let mut t = Transformer::new(&policy.for_fastq()).expect("FASTQ policy");
        let mut record = records("@r runid=x\tMM:Z:C+m?,0;\tML:B:C,7 ch=3\nAC\n+\nII\n").remove(0);
        t.apply(&mut record).expect("valid");
        assert_eq!(record.header, b"r runid=x ch=3");
    }

    #[test]
    fn cpg_methylation_needs_alignments() {
        let mut policy = Policy::preset(Platform::Ont, Mode::Safe);
        policy.modifications.calls = ModificationCalls::Cpg;
        assert!(Transformer::new(&policy.for_fastq()).is_err());
    }

    #[test]
    fn recognizes_fastq_names() {
        for name in ["a.fastq", "a.FQ", "a.fastq.gz", "dir/a.fq.gz"] {
            assert!(is_fastq_path(Path::new(name)), "{name}");
        }
        for name in ["a.bam", "a.cram", "a.gz", "fastq"] {
            assert!(!is_fastq_path(Path::new(name)), "{name}");
        }
    }
}
