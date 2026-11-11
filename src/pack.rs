//! Streams a BAM or CRAM file through a [`Transformer`] into a CRAM or BAM archive.

use std::{
    collections::HashMap,
    ffi::OsString,
    fs::File,
    io::{self, Read},
    mem,
    num::NonZero,
    path::{Path, PathBuf},
    sync::Arc,
};

use anyhow::{Context, bail};
use pasta::{
    bam,
    bgzf::{self, io::writer::CompressionLevel},
    cram::{
        self,
        file_definition::Version,
        io::writer::compression::{Options, Profile},
    },
    csi, fasta,
    sam::{
        self,
        header::record::value::{
            Map,
            map::{Program, header::tag as header_tag, program::tag, read_group},
        },
    },
};

use crate::{
    alignments::{Alignments, CRAM_MAGIC, Source},
    homopolymer::Runs,
    pipeline,
    policy::{Mode, ModificationCalls, Overrides, Platform, Policy, Quality},
    raw,
    transform::{Outcome, RecordContext, Stats, Transformer},
};

/// Records are transformed in parallel, in batches of about this many encoded bytes.
pub(crate) const BATCH_BYTES: usize = 4 << 20;

/// Where the input comes from.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Input {
    /// A BAM, CRAM, or FASTQ file.
    Path(PathBuf),
    /// BAM, CRAM, or FASTQ on standard input.
    Stdin,
}

/// The archive format.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Format {
    /// CRAM: smallest, and needs the reference to decode. CRAM 3.0 by default, CRAM 3.1
    /// (archive profile) with [`Mode::Max`].
    Cram,
    /// BAM: larger, readable by every tool without a reference.
    Bam,
    /// FASTQ (BGZF-compressed), from FASTQ input: for reads without an alignment.
    Fastq,
}

impl Format {
    /// FASTQ for a `.fastq`/`.fq` path (with or without `.gz`), BAM for `.bam` (any case), and
    /// CRAM otherwise.
    pub fn from_path(path: &Path) -> Self {
        if crate::fastq::is_fastq_path(path) {
            return Self::Fastq;
        }
        match path.extension().and_then(|e| e.to_str()) {
            Some(ext) if ext.eq_ignore_ascii_case("bam") => Self::Bam,
            _ => Self::Cram,
        }
    }
}

/// The BGZF level of BAM output. Level 12 ([`Mode::Max`]) is 8-10% smaller but 6-30x slower
/// on long reads (measured in the README).
fn bam_level(mode: Mode) -> CompressionLevel {
    match mode {
        Mode::Safe => CompressionLevel::default(),
        Mode::Max => CompressionLevel::BEST,
    }
}

/// Bases per slice for `--max` when per-base qualities are stored. The archive profile's
/// 50 Mbases make the quality block's codec trials long, and pasta waits for trials in
/// container order. On ONT chr20 (without lzma and bzip2), 10 Mbases was 1% smaller, faster
/// (49 s against 62 s), and used a third of the memory. With per-read qualities the quality
/// block is tiny and 10 Mbases only cost size (PacBio +3%), so the profile's slices stay.
const MAX_BASES_PER_SLICE_PER_BASE_QUALITIES: usize = 10_000_000;

/// The CRAM version and profile: 3.0 is read by every htslib-based tool, 3.1's codecs need
/// htslib 1.12 or newer (some tools bundle older decoders).
///
/// `--max` leaves out the archive profile's lzma and bzip2, the slowest codecs to try: on
/// chr20 that cut CPU by about 70% (ONT, 2.5x faster) and 65% (PacBio, 2.2x faster) for 8%
/// (ONT) and 19% (PacBio) larger output.
fn cram_options(mode: Mode, policy: &Policy) -> Options {
    match mode {
        Mode::Safe => Options::new(Version::new(3, 0), Profile::Normal),
        Mode::Max => {
            let mut options = Options::new(Version::new(3, 1), Profile::Archive);
            options.use_lzma = false;
            options.use_bzip2 = false;
            if matches!(policy.quality, Quality::Keep | Quality::Bins { .. }) {
                options.bases_per_slice = MAX_BASES_PER_SLICE_PER_BASE_QUALITIES;
            }
            options
        }
    }
}

/// Everything [`pack`] needs.
#[derive(Clone, Debug)]
pub struct PackOptions {
    /// The BAM input.
    pub input: Input,
    /// The output path.
    pub output: PathBuf,
    /// The output format.
    pub format: Format,
    /// The indexed FASTA the reads are aligned to. Required for CRAM output and for
    /// homopolymer protection; otherwise unused.
    pub reference: Option<PathBuf>,
    /// Safe (default) or maximum compression: picks the built-in policy and the container.
    pub mode: Mode,
    /// The platform, when the header's `@RG` `PL` fields do not name exactly one.
    pub platform: Option<Platform>,
    /// A policy to use instead of the built-in one (the container still follows `mode`).
    pub custom_policy: Option<Policy>,
    /// Which base-modification calls to keep, overriding the policy's.
    pub methylation: Option<ModificationCalls>,
    /// Keep `ML` probabilities exactly, overriding the policy's `ml_levels`.
    pub keep_ml: bool,
    /// Threads for BGZF decompression and output compression.
    pub threads: NonZero<usize>,
    /// Whether to index the output (CRAI for CRAM; BAI, or CSI for long references, for BAM).
    pub write_index: bool,
    /// The full command line, for the `@PG` record.
    pub command_line: String,
}

impl PackOptions {
    /// The command line's changes to the policy.
    pub(crate) fn overrides(&self) -> Overrides {
        Overrides {
            methylation: self.methylation,
            keep_ml: self.keep_ml,
        }
    }
}

/// What [`pack`] did: the policy it applied and its counts.
#[derive(Clone, Debug)]
pub struct Report {
    /// The policy applied.
    pub policy: Policy,
    /// The platform whose built-in policy was used; `None` for a custom policy.
    pub platform: Option<Platform>,
    /// What was changed.
    pub stats: Stats,
    /// Whether an index was written.
    pub index: IndexOutcome,
}

/// What happened to the index.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IndexOutcome {
    /// Written to this path.
    Written(PathBuf),
    /// Not asked for.
    NotRequested,
    /// Not written because the records are not in coordinate order (reference sequences in
    /// header order, positions ascending), which region queries rely on.
    SkippedUnsorted,
    /// Not written because the BAM header does not declare `SO:coordinate`, which a BAM index
    /// needs.
    SkippedUndeclaredOrder,
}

/// Tracks whether records arrive in coordinate order: reference IDs ascending, positions
/// ascending within each, unplaced reads last.
#[derive(Clone, Copy, Debug)]
struct CoordinateOrder {
    last: (u64, u64),
    sorted: bool,
}

impl Default for CoordinateOrder {
    fn default() -> Self {
        Self {
            last: (0, 0),
            sorted: true,
        }
    }
}

impl CoordinateOrder {
    fn observe(&mut self, buf: &[u8]) {
        let key = match raw::reference_sequence_id(buf) {
            Some(id) => (id as u64, raw::position(buf).map_or(0, |p| p as u64)),
            None => (u64::MAX, 0),
        };
        self.sorted &= key >= self.last;
        self.last = key;
    }

    fn is_sorted(&self) -> bool {
        self.sorted
    }
}

/// The platform named by the header's `@RG` `PL` fields (`ONT` or `PACBIO`, any case). `None`
/// when there is no such platform or more than one.
pub fn detect_platform(header: &sam::Header) -> Option<Platform> {
    let mut platforms = header.read_groups().values().filter_map(|rg| {
        let platform = rg.other_fields().get(&read_group::tag::PLATFORM)?;
        match platform.to_ascii_uppercase().as_slice() {
            b"ONT" => Some(Platform::Ont),
            b"PACBIO" => Some(Platform::Pacbio),
            _ => None,
        }
    });
    let first = platforms.next()?;
    platforms.all(|p| p == first).then_some(first)
}

/// The policy for `options` and the input's header: the custom policy if given, otherwise the
/// built-in one for the platform (given, or detected from the header) and mode.
///
/// # Errors
///
/// Fails when no platform is given or detected, or a custom policy is invalid.
pub fn resolve_policy(
    options: &PackOptions,
    header: &sam::Header,
) -> anyhow::Result<(Policy, Option<Platform>)> {
    if let Some(policy) = &options.custom_policy {
        policy.validate().map_err(anyhow::Error::msg)?;
        return Ok((policy.clone().with_overrides(options.overrides()), None));
    }
    let platform = options
        .platform
        .or_else(|| detect_platform(header))
        .context(
            "cannot tell the platform: the header's @RG PL fields name neither ONT nor PACBIO \
         (or name both); pass --platform ont or --platform pacbio",
        )?;
    Ok((
        Policy::preset(platform, options.mode).with_overrides(options.overrides()),
        Some(platform),
    ))
}

/// The archive being written.
enum Output {
    Cram(Box<cram::io::Writer<File>>),
    Bam(bam::io::Writer<bgzf::io::MultithreadedWriter<File>>),
}

impl Output {
    fn create(
        options: &PackOptions,
        policy: &Policy,
        repository: Option<&fasta::Repository>,
        header: &sam::Header,
    ) -> anyhow::Result<Self> {
        let path = &options.output;
        let create = || format!("creating {}", path.display());
        let mut output = match options.format {
            Format::Cram => {
                let repository = repository.context("CRAM output needs the reference (-T)")?;
                let writer = cram::io::writer::Builder::default()
                    .set_reference_sequence_repository(repository.clone())
                    .set_compression_options(cram_options(options.mode, policy))
                    .set_worker_count(options.threads)
                    .build_from_path(path)
                    .with_context(create)?;
                Self::Cram(Box::new(writer))
            }
            Format::Bam => {
                let bgzf = bgzf::io::multithreaded_writer::Builder::default()
                    .set_compression_level(bam_level(options.mode))
                    .set_worker_count(options.threads)
                    .build_from_writer(File::create(path).with_context(create)?);
                Self::Bam(bam::io::Writer::from(bgzf))
            }
            Format::Fastq => bail!("FASTQ is written from FASTQ input"),
        };
        match &mut output {
            Self::Cram(w) => w.write_header(header),
            Self::Bam(w) => w.write_header(header),
        }
        .context("writing the header")?;
        Ok(output)
    }

    fn write(&mut self, header: &sam::Header, record: &bam::Record) -> io::Result<()> {
        match self {
            Self::Cram(w) => w.write_bam_record(header, record),
            Self::Bam(w) => w.write_record(header, record),
        }
    }

    fn finish(self, header: &sam::Header) -> io::Result<()> {
        match self {
            Self::Cram(mut w) => w.try_finish(header),
            Self::Bam(w) => w.into_inner().finish().map(drop),
        }
    }
}

/// Whether the header's `@HD` line says `SO:coordinate`.
fn declares_coordinate_order(header: &sam::Header) -> bool {
    header
        .header()
        .and_then(|hd| hd.other_fields().get(&header_tag::SORT_ORDER))
        .is_some_and(|order| order.as_slice() == b"coordinate")
}

/// `path` with `suffix` appended, as indexes are named.
fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut p = OsString::from(path);
    p.push(suffix);
    PathBuf::from(p)
}

/// The names an index of a `format` archive at `path` may have.
fn index_paths(path: &Path, format: Format) -> Vec<PathBuf> {
    let suffixes: &[&str] = match format {
        Format::Cram => &[".crai"],
        Format::Bam => &[".bai", ".csi"],
        Format::Fastq => &[],
    };
    suffixes.iter().map(|s| with_suffix(path, s)).collect()
}

/// Removes the indexes of whatever file was at `path` before, which no longer describe it.
fn remove_stale_indexes(path: &Path, format: Format) -> anyhow::Result<()> {
    for index in index_paths(path, format) {
        match std::fs::remove_file(&index) {
            Err(e) if e.kind() != io::ErrorKind::NotFound => {
                return Err(e).with_context(|| format!("removing the stale {}", index.display()));
            }
            _ => {}
        }
    }
    Ok(())
}

/// Indexes a finished archive beside it.
fn write_index(path: &Path, format: Format) -> anyhow::Result<PathBuf> {
    let with_suffix = |suffix: &str| with_suffix(path, suffix);
    let indexing = || format!("indexing {}", path.display());
    match format {
        Format::Cram => {
            let index = cram::fs::index(path).with_context(indexing)?;
            let dst = with_suffix(".crai");
            cram::crai::fs::write(&dst, &index).context("writing the CRAI index")?;
            Ok(dst)
        }
        // pasta builds a BAI, or a CSI when a reference sequence is too long for BAI (> 2^29).
        Format::Bam => match bam::fs::index(path).with_context(indexing)? {
            bam::Index::Bai(index) => {
                let dst = with_suffix(".bai");
                bam::bai::fs::write(&dst, &index).context("writing the BAI index")?;
                Ok(dst)
            }
            bam::Index::Csi(index) => {
                let dst = with_suffix(".csi");
                csi::fs::write(&dst, &index).context("writing the CSI index")?;
                Ok(dst)
            }
        },
        Format::Fastq => bail!("FASTQ archives are not indexed"),
    }
}

/// Packs `options.input` into an archive.
///
/// # Errors
///
/// Fails on unreadable input or reference, a missing reference that the format or policy
/// needs, malformed records, an unresolvable policy, or write errors.
pub fn pack(options: &PackOptions) -> anyhow::Result<Report> {
    if let Input::Path(path) = &options.input
        && same_file(path, &options.output)
    {
        bail!(
            "the output {} is the input; writing it would destroy the input before it is read",
            options.output.display()
        );
    }
    if options.format == Format::Fastq {
        return crate::fastq::pack(options);
    }
    if let Input::Path(path) = &options.input
        && crate::fastq::is_fastq_path(path)
    {
        bail!("FASTQ input is packed to FASTQ; name the output .fastq.gz");
    }
    pack_in_batches(options, BATCH_BYTES)
}

/// Whether `a` and `b` are the same existing file, under any names (links included).
fn same_file(a: &Path, b: &Path) -> bool {
    let (Ok(a_meta), Ok(b_meta)) = (std::fs::metadata(a), std::fs::metadata(b)) else {
        return false;
    };
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        a_meta.dev() == b_meta.dev() && a_meta.ino() == b_meta.ino()
    }
    #[cfg(not(unix))]
    {
        let _ = (a_meta, b_meta);
        matches!(
            (std::fs::canonicalize(a), std::fs::canonicalize(b)),
            (Ok(a), Ok(b)) if a == b
        )
    }
}

/// [`pack`], transforming batches of about `batch_bytes` encoded bytes.
fn pack_in_batches(options: &PackOptions, batch_bytes: usize) -> anyhow::Result<Report> {
    let threads = options.threads;
    let source: Source = match &options.input {
        Input::Path(path) => {
            Box::new(File::open(path).with_context(|| format!("opening {}", path.display()))?)
        }
        Input::Stdin => Box::new(io::stdin()),
    };
    let repository = options
        .reference
        .as_deref()
        .map(open_reference)
        .transpose()?;
    let (mut input, mut header) = Alignments::open(source, threads, repository.as_ref())?;
    let (policy, platform) = resolve_policy(options, &header)?;
    annotate_header(&mut header, &policy, &options.command_line)?;

    let transformer = Transformer::new(&policy);
    if transformer.homopolymer_length().is_some() && options.reference.is_none() {
        bail!("homopolymer protection needs the reference (-T)");
    }
    if transformer.needs_reference() && options.reference.is_none() {
        bail!("keeping methylation calls at CpG sites needs the reference (-T)");
    }
    let mut writer = Output::create(options, &policy, repository.as_ref(), &header)?;
    remove_stale_indexes(&options.output, options.format)?;
    let mut references = ReferenceCache::new(
        transformer.homopolymer_length(),
        transformer.needs_reference(),
        repository,
    );
    let mut order = CoordinateOrder::default();

    // CRAM is decoded while records are rewritten and written, so its decoders take half the
    // threads, keeping CPU use near -@.
    let decoders = if input.is_cram() {
        NonZero::new(threads.get().div_ceil(2)).expect("at least one thread")
    } else {
        threads
    };
    let stats = transform_in_parallel(
        input.records(&header, decoders),
        &header,
        &mut references,
        &transformer,
        threads,
        batch_bytes,
        |record| {
            order.observe(record.as_ref());
            writer.write(&header, record).context("writing a record")
        },
    )?;
    writer.finish(&header).context("finishing the output")?;

    let index = if !options.write_index {
        IndexOutcome::NotRequested
    } else if !order.is_sorted() {
        IndexOutcome::SkippedUnsorted
    } else if options.format == Format::Bam && !declares_coordinate_order(&header) {
        IndexOutcome::SkippedUndeclaredOrder
    } else {
        IndexOutcome::Written(write_index(&options.output, options.format)?)
    };

    Ok(Report {
        index,
        policy,
        platform,
        stats,
    })
}

/// The end-of-file container that ends every CRAM 3 file (CRAM specification, section 9).
const CRAM_EOF: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

/// Whether the file at `path` is BGZF (a BAM, or bgzipped FASTQ) or CRAM 3 that does not end
/// with its end-of-file marker, as a copy cut short between two blocks or containers does not.
/// Such a file reads as complete, so this is the only sign of the records it lost. Only regular
/// files are checked: reading a pipe would take the bytes its reader needs.
///
/// # Errors
///
/// Fails if the file cannot be read.
pub fn missing_eof_marker(path: &Path) -> io::Result<bool> {
    use io::{Seek, SeekFrom};
    const BGZF_EOF_LEN: usize = 28;
    if !std::fs::metadata(path)?.is_file() {
        return Ok(false);
    }
    let mut file = File::open(path)?;
    let mut head = Vec::with_capacity(18);
    file.by_ref().take(18).read_to_end(&mut head)?;
    // CRAM 2 and earlier end differently.
    let is_cram = head.starts_with(CRAM_MAGIC) && head.get(4) == Some(&3);
    let eof_len = if is_cram {
        CRAM_EOF.len()
    } else if crate::fastq::is_bgzf(&head) {
        BGZF_EOF_LEN
    } else {
        return Ok(false);
    };
    let len = file.metadata()?.len();
    let Some(start) = len.checked_sub(eof_len as u64) else {
        return Ok(true);
    };
    file.seek(SeekFrom::Start(start))?;
    let mut tail = vec![0; eof_len];
    file.read_exact(&mut tail)?;
    Ok(if is_cram {
        tail != CRAM_EOF
    } else {
        !bgzf::io::raw::is_eof_marker(&tail)
    })
}

/// Reads records, transforms them on `threads` workers, and passes the written ones to `write`
/// in input order. A reader thread batches records of about `batch_bytes` encoded bytes; the
/// calling thread runs `write`.
///
/// # Errors
///
/// Fails on unreadable or malformed records, a reference sequence that cannot be loaded, an error
/// from `write`, or a panicking worker. Every thread has stopped when this returns.
pub(crate) fn transform_in_parallel(
    records: impl Iterator<Item = anyhow::Result<Vec<u8>>> + Send,
    header: &sam::Header,
    references: &mut ReferenceCache,
    transformer: &Transformer,
    threads: NonZero<usize>,
    batch_bytes: usize,
    mut write: impl FnMut(&bam::Record) -> anyhow::Result<()>,
) -> anyhow::Result<Stats> {
    let workers = vec![transformer.clone(); threads.get()];
    let workers = pipeline::ordered(
        |send| read_batches(records, header, references, send, batch_bytes),
        workers,
        transform_batch,
        |records: Vec<bam::Record>| records.iter().try_for_each(&mut write),
    )?;
    Ok(workers.iter().fold(Stats::default(), |mut stats, worker| {
        stats += worker.stats();
        stats
    }))
}

/// A run of encoded records with the reference data they need.
#[derive(Default)]
struct Batch {
    records: Vec<Vec<u8>>,
    references: Vec<(usize, LoadedReference)>,
}

/// Reads records into batches and sends them until the input ends or the workers stop.
fn read_batches(
    records: impl Iterator<Item = anyhow::Result<Vec<u8>>>,
    header: &sam::Header,
    references: &mut ReferenceCache,
    send: &mut dyn FnMut(Batch) -> bool,
    batch_bytes: usize,
) -> anyhow::Result<()> {
    let mut batch = Batch::default();
    let mut bytes = 0;
    for buf in records {
        let buf = buf?;
        if let Some(id) = raw::reference_sequence_id(&buf)
            && !batch.references.iter().any(|(i, _)| *i == id)
        {
            batch.references.push((id, references.load(header, id)?));
        }
        bytes += buf.len();
        batch.records.push(buf);
        if bytes >= batch_bytes {
            if !send(mem::take(&mut batch)) {
                // The workers stopped; the pipeline reports why.
                return Ok(());
            }
            bytes = 0;
        }
    }
    if !batch.records.is_empty() {
        send(batch);
    }
    Ok(())
}

fn transform_batch(
    transformer: &mut Transformer,
    batch: Batch,
) -> anyhow::Result<Vec<bam::Record>> {
    let mut records = Vec::with_capacity(batch.records.len());
    for mut buf in batch.records {
        let context = raw::reference_sequence_id(&buf)
            .and_then(|id| batch.references.iter().find(|(i, _)| *i == id))
            .map_or_else(RecordContext::default, |(_, loaded)| loaded.context());
        let outcome = transformer
            .apply(&mut buf, context)
            .context("transforming a record")?;
        if outcome == Outcome::Write {
            records.push(bam::Record::try_from(buf).context("re-encoding a record")?);
        }
    }
    Ok(records)
}

/// Adds lrsqueeze's `@PG` record and `@CO` lines describing the policy (in words, and as one line of
/// TOML), so the archive says how it was made.
fn annotate_header(
    header: &mut sam::Header,
    policy: &Policy,
    command_line: &str,
) -> anyhow::Result<()> {
    let program = Map::<Program>::builder()
        .insert(tag::NAME, env!("CARGO_PKG_NAME"))
        .insert(tag::VERSION, env!("CARGO_PKG_VERSION"))
        .insert(tag::COMMAND_LINE, command_line)
        .build()
        .context("building the @PG record")?;
    header
        .add_program_record(env!("CARGO_PKG_NAME"), program)
        .context("adding the @PG record")?;
    let comments = header.comments_mut();
    comments.push(format!("lrsqueeze lossy policy: {}", policy.summary()).into());
    // The same policy, machine-readable, for `lrsqueeze verify`.
    comments.push(policy.to_header_line().into());
    Ok(())
}

pub(crate) fn open_reference(path: &Path) -> anyhow::Result<fasta::Repository> {
    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .with_context(|| {
            format!(
                "opening the reference {} (it needs a .fai index)",
                path.display()
            )
        })?;
    Ok(fasta::Repository::new(
        fasta::repository::adapters::IndexedReader::new(reader),
    ))
}

/// What records need from their reference sequences: homopolymer runs (found once per
/// sequence), and for CpG context the sequence itself (only the current one is kept).
pub(crate) struct ReferenceCache {
    min_length: Option<usize>,
    keep_sequence: bool,
    repository: Option<fasta::Repository>,
    runs: HashMap<usize, Arc<Runs>>,
    current: Option<(usize, Arc<fasta::record::Sequence>)>,
}

/// What records on one reference sequence get from [`ReferenceCache`].
#[derive(Clone, Debug, Default)]
struct LoadedReference {
    runs: Option<Arc<Runs>>,
    sequence: Option<Arc<fasta::record::Sequence>>,
}

impl LoadedReference {
    fn context(&self) -> RecordContext<'_> {
        RecordContext {
            runs: self.runs.as_deref(),
            reference: self.sequence.as_deref().map(|s| s.as_ref()),
        }
    }
}

impl ReferenceCache {
    /// `repository` must be set when either `min_length` or `keep_sequence` is.
    pub(crate) fn new(
        min_length: Option<usize>,
        keep_sequence: bool,
        repository: Option<fasta::Repository>,
    ) -> Self {
        Self {
            min_length,
            keep_sequence,
            repository,
            runs: HashMap::new(),
            current: None,
        }
    }

    /// What records on reference sequence `id` need.
    fn load(&mut self, header: &sam::Header, id: usize) -> anyhow::Result<LoadedReference> {
        let need_runs = self.min_length.is_some() && !self.runs.contains_key(&id);
        let need_sequence =
            self.keep_sequence && self.current.as_ref().is_none_or(|(c, _)| *c != id);
        if need_runs || need_sequence {
            let sequence = self.sequence(header, id)?;
            if let Some(min_length) = self.min_length.filter(|_| need_runs) {
                self.runs.insert(
                    id,
                    Arc::new(Runs::find(sequence.as_ref().as_ref(), min_length)),
                );
            }
            if need_sequence {
                self.current = Some((id, sequence));
            }
        }
        Ok(LoadedReference {
            runs: self.runs.get(&id).cloned(),
            sequence: self
                .current
                .as_ref()
                .filter(|(c, _)| *c == id)
                .map(|(_, s)| Arc::clone(s)),
        })
    }

    /// Reference sequence `id`, read from the FASTA.
    fn sequence(
        &self,
        header: &sam::Header,
        id: usize,
    ) -> anyhow::Result<Arc<fasta::record::Sequence>> {
        let (name, _) = header
            .reference_sequences()
            .get_index(id)
            .with_context(|| {
                format!("record on reference sequence {id}, which the header lacks")
            })?;
        self.repository
            .as_ref()
            .expect("pack checks the reference is given when it is needed")
            .get(name)
            .with_context(|| format!("reference sequence {name} is not in the FASTA"))?
            .with_context(|| format!("reading reference sequence {name}"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::ensure;
    use std::{fs, sync::mpsc, thread};

    /// A BAM of 200 unmapped 1 kb PacBio reads, `r0` to `r199`, in several BGZF blocks; read
    /// `bad`, if given, has malformed `MM`/`ML`.
    fn bam_bytes(bad: Option<usize>) -> anyhow::Result<Vec<u8>> {
        let mut sam = String::from("@HD\tVN:1.6\n@SQ\tSN:chrT\tLN:8\n@RG\tID:rg\tPL:PACBIO\n");
        let (bases, quals) = ("ACGT".repeat(250), "+".repeat(1000));
        for i in 0..200 {
            let mods = if Some(i) == bad {
                // The read has 250 Cs.
                "\tMM:Z:C+m?,999;\tML:B:C,1"
            } else {
                ""
            };
            sam.push_str(&format!(
                "r{i}\t4\t*\t0\t0\t*\t*\t0\t0\t{bases}\t{quals}\tRG:Z:rg{mods}\n"
            ));
        }
        let mut reader = sam::io::Reader::new(sam.as_bytes());
        let header = reader.read_header()?;
        let mut writer = bam::io::Writer::new(Vec::new());
        writer.write_header(&header)?;
        for record in reader.record_bufs(&header) {
            use sam::alignment::io::Write;
            writer.write_alignment_record(&header, &record?)?;
        }
        Ok(writer.into_inner().finish()?)
    }

    /// [`bam_bytes`] as `in.bam` in `dir`, with a reference.
    fn write_bam(dir: &Path, bad: Option<usize>) -> anyhow::Result<()> {
        fs::write(dir.join("ref.fa"), ">chrT\nACGTACGT\n")?;
        fs::write(dir.join("ref.fa.fai"), "chrT\t8\t6\t8\t9\n")?;
        fs::write(dir.join("in.bam"), bam_bytes(bad)?)?;
        Ok(())
    }

    /// Runs the parallel transform of `bam` (PacBio default policy, 4 workers, one record per
    /// batch) on another thread, failing the test if it does not finish within a minute.
    fn transform_with_timeout(
        bam: Vec<u8>,
        write: impl FnMut(&bam::Record) -> anyhow::Result<()> + Send + 'static,
    ) -> anyhow::Result<Stats> {
        let (done, finished) = mpsc::channel();
        thread::spawn(move || {
            let run = || {
                let (mut input, header) =
                    Alignments::open(Box::new(io::Cursor::new(bam)), NonZero::<usize>::MIN, None)?;
                let transformer = Transformer::new(&Policy::preset(Platform::Pacbio, Mode::Safe));
                transform_in_parallel(
                    input.records(&header, NonZero::<usize>::MIN),
                    &header,
                    &mut ReferenceCache::new(None, false, None),
                    &transformer,
                    NonZero::new(4).expect("4 is nonzero"),
                    1,
                    write,
                )
            };
            // The receiver is gone only if the test already failed.
            let _ = done.send(run());
        });
        finished
            .recv_timeout(std::time::Duration::from_secs(60))
            .expect("the transform stops within a minute")
    }

    #[test]
    fn a_failing_write_stops_every_thread() -> anyhow::Result<()> {
        let mut written = 0;
        let error = transform_with_timeout(bam_bytes(None)?, move |_| {
            written += 1;
            ensure!(written < 50, "disk full");
            Ok(())
        })
        .expect_err("the 50th write fails");
        assert_eq!(format!("{error:#}"), "disk full");
        Ok(())
    }

    #[test]
    fn a_missing_eof_marker_is_found() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        let path = dir.path().join("in.bam");
        let bam = bam_bytes(None)?;
        fs::write(&path, &bam)?;
        assert!(!missing_eof_marker(&path)?, "a complete BAM");
        // Cut off the end-of-file marker: the BAM still reads to its last whole block.
        fs::write(&path, &bam[..bam.len() - 28])?;
        assert!(missing_eof_marker(&path)?, "a BAM cut short between blocks");
        fs::write(&path, &bam[..10])?;
        assert!(!missing_eof_marker(&path)?, "too short to be BGZF");
        fs::write(&path, "@r\nACGT\n+\nIIII\n")?;
        assert!(!missing_eof_marker(&path)?, "plain FASTQ");
        Ok(())
    }

    #[test]
    fn truncated_input_is_an_error() -> anyhow::Result<()> {
        let mut bam = bam_bytes(None)?;
        // Cut into the last data block, past the end-of-file marker.
        bam.truncate(bam.len() - 60);
        let error = transform_with_timeout(bam, |_| Ok(())).expect_err("the input is truncated");
        assert!(
            format!("{error:#}").contains("reading a BAM record"),
            "{error:#}"
        );
        Ok(())
    }

    fn options(dir: &Path, methylation: Option<ModificationCalls>) -> PackOptions {
        PackOptions {
            input: Input::Path(dir.join("in.bam")),
            output: dir.join("out.bam"),
            format: Format::Bam,
            reference: Some(dir.join("ref.fa")),
            mode: Mode::Safe,
            platform: None,
            custom_policy: None,
            methylation,
            keep_ml: false,
            threads: NonZero::new(4).expect("4 is nonzero"),
            write_index: false,
            command_line: "test".to_owned(),
        }
    }

    #[test]
    fn batches_are_written_in_input_order() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        write_bam(dir.path(), None)?;
        // One record per batch, spread over 4 workers.
        let report = pack_in_batches(&options(dir.path(), None), 1)?;
        assert_eq!(
            (report.stats.records_in, report.stats.records_out),
            (200, 200)
        );
        let mut reader = bam::io::reader::Builder.build_from_path(dir.path().join("out.bam"))?;
        let header = reader.read_header()?;
        let names: Vec<String> = reader
            .record_bufs(&header)
            .map(|r| Ok(r?.name().map(|n| n.to_string()).unwrap_or_default()))
            .collect::<io::Result<_>>()?;
        let want: Vec<String> = (0..200).map(|i| format!("r{i}")).collect();
        assert_eq!(names, want);
        Ok(())
    }

    #[test]
    fn a_failing_record_stops_every_thread() -> anyhow::Result<()> {
        let dir = tempfile::tempdir()?;
        write_bam(dir.path(), Some(150))?;
        let error = pack_in_batches(&options(dir.path(), Some(ModificationCalls::CpgOnly)), 1)
            .expect_err("MM points past the read's Cs");
        assert!(
            format!("{error:#}").contains("transforming a record"),
            "{error:#}"
        );
        Ok(())
    }
}
