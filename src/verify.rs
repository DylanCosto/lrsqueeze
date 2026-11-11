//! Checks that an archive holds exactly what its policy makes of the original BAM or CRAM.
//!
//! The original's records go through the same transform as [`pack`](crate::pack::pack), and each
//! result is compared with the archive's record in the same position. Two differences are
//! expected of CRAM and allowed: `=`/`X` CIGAR operations decode as `M`, and decoding adds the
//! `MD`/`NM` tags a record lacked. Everything else must match: names, flags, positions, mapping
//! quality, CIGAR, mate fields, sequence, qualities, and every tag (in any order).

use std::{
    collections::{BTreeSet, VecDeque},
    fs::File,
    mem,
    num::NonZero,
    path::PathBuf,
    sync::mpsc::{self, Receiver, SyncSender},
    thread,
};

use anyhow::{Context, bail};
use pasta::{
    bam, bgzf, cram,
    sam::{
        self,
        alignment::{
            RecordBuf,
            record::{cigar::op::Kind, data::field::Tag},
            record_buf::data::field::{Value, value::Array},
        },
    },
};

use crate::{
    alignments::Alignments,
    pack::{BATCH_BYTES, Format, ReferenceCache, open_reference, transform_in_parallel},
    policy::Policy,
    transform::Transformer,
};

/// Archive records are read on their own thread (CRAM slices decoded on a pool of workers) and
/// passed on in batches of this many.
const DECODED_BATCH: usize = 1024;

/// At most this many differences are described; the rest are only counted.
const MAX_DESCRIBED: usize = 10;

/// What to verify.
#[derive(Clone, Debug)]
pub struct VerifyOptions {
    /// The BAM or CRAM that was packed.
    pub original: PathBuf,
    /// The archive made from it.
    pub archive: PathBuf,
    /// The archive's format.
    pub format: Format,
    /// The reference the reads are aligned to (needed for CRAM and for reference-based rules).
    pub reference: Option<PathBuf>,
    /// The policy to check against, for archives that do not record theirs.
    pub policy: Option<Policy>,
    /// Threads for reading the original and transforming its records.
    pub threads: NonZero<usize>,
}

/// One record that is not what the policy makes of the original.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Difference {
    /// The record's position in the archive, from 0.
    pub index: u64,
    /// The original record's name.
    pub name: String,
    /// What differs.
    pub what: String,
}

/// The outcome of [`verify`].
#[derive(Clone, Debug)]
#[must_use]
pub struct VerifyReport {
    /// The policy checked against.
    pub policy: Policy,
    /// Header lines of the original that the archive lacks or changed.
    pub header_differences: Vec<String>,
    /// Records compared.
    pub records: u64,
    /// Records that differ, or are missing from or extra in the archive.
    pub differences: u64,
    /// The first few differences.
    pub described: Vec<Difference>,
}

impl VerifyReport {
    /// Whether the header kept the original's lines and every record matched.
    pub fn is_match(&self) -> bool {
        self.header_differences.is_empty() && self.differences == 0
    }

    pub(crate) fn new(policy: Policy, header_differences: Vec<String>) -> Self {
        Self {
            policy,
            header_differences,
            records: 0,
            differences: 0,
            described: Vec::new(),
        }
    }

    pub(crate) fn record(&mut self, index: u64, name: impl Into<String>, what: impl Into<String>) {
        self.differences += 1;
        if self.described.len() < MAX_DESCRIBED {
            self.described.push(Difference {
                index,
                name: name.into(),
                what: what.into(),
            });
        }
    }
}

/// The archive being read.
enum Archive {
    Cram(Box<cram::io::Reader<File>>),
    Bam(Box<bam::io::Reader<bgzf::io::MultithreadedReader<File>>>),
}

/// Compares `options.archive` with what its policy makes of `options.original`.
///
/// # Errors
///
/// Fails when either file cannot be read, the archive records no policy and none is given, the
/// two files have different reference sequences, or the reference is missing where it is needed.
/// Differences between records are not errors; they are in the report.
pub fn verify(options: &VerifyOptions) -> anyhow::Result<VerifyReport> {
    if options.format == Format::Fastq {
        return crate::fastq::verify(options);
    }
    let open =
        |path: &PathBuf| File::open(path).with_context(|| format!("opening {}", path.display()));
    let repository = options
        .reference
        .as_deref()
        .map(open_reference)
        .transpose()?;
    let (mut original, original_header) = Alignments::open(
        Box::new(open(&options.original)?),
        options.threads,
        repository.as_ref(),
    )
    .context("reading the original")?;

    let (mut archive, archive_header) = match options.format {
        Format::Cram => {
            let repository = repository
                .clone()
                .context("reading a CRAM archive needs the reference (-T)")?;
            let mut reader = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repository)
                .build_from_reader(open(&options.archive)?);
            let header = reader.read_header().context("reading the archive header")?;
            (Archive::Cram(Box::new(reader)), header)
        }
        Format::Fastq => unreachable!("FASTQ archives are verified by fastq::verify"),
        Format::Bam => {
            let mut reader =
                bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(
                    options.threads,
                    open(&options.archive)?,
                ));
            let header = reader.read_header().context("reading the archive header")?;
            (Archive::Bam(Box::new(reader)), header)
        }
    };
    let policy = match (options.policy.clone(), recorded_policy(&archive_header)?) {
        (Some(given), Some(recorded)) if given != recorded => bail!(
            "the archive records a different policy from the one given (-c or --platform); \
             leave them out to check against the recorded policy:\n  given:    {}\n  \
             recorded: {}",
            given.summary(),
            recorded.summary()
        ),
        (Some(policy), _) | (None, Some(policy)) => policy,
        (None, None) => bail!(
            "the archive does not record its policy (it was made by an older lrsqueeze or another \
             tool); pass the --platform (and --max or --methylation) it was packed with, or the \
             policy file"
        ),
    };
    policy.validate().map_err(anyhow::Error::msg)?;
    if names_and_lengths(&original_header) != names_and_lengths(&archive_header) {
        bail!("the original and the archive have different reference sequences (@SQ lines)");
    }

    let transformer = Transformer::new(&policy);
    let needs_reference =
        transformer.homopolymer_length().is_some() || transformer.needs_reference();
    if needs_reference && repository.is_none() {
        bail!("this policy uses the reference (-T)");
    }
    let mut references = ReferenceCache::new(
        transformer.homopolymer_length(),
        transformer.needs_reference(),
        repository,
    );
    let mut report = VerifyReport::new(
        policy,
        header_differences(&original_header, &archive_header),
    );

    thread::scope(|scope| -> anyhow::Result<()> {
        let (batches, decoded) = mpsc::sync_channel(4);
        let header = &archive_header;
        // Two CRAM files share the decoding threads, so that verify stays near -@ threads.
        let decoders = if original.is_cram() && options.format == Format::Cram {
            NonZero::new(options.threads.get().div_ceil(2)).expect("at least one thread")
        } else {
            options.threads
        };
        let decoder = scope.spawn(move || decode(&mut archive, header, decoders, &batches));
        let mut stream = Decoded {
            batches: decoded,
            ready: VecDeque::new(),
        };
        let compared = transform_in_parallel(
            original.records(&original_header, decoders),
            &original_header,
            &mut references,
            &transformer,
            options.threads,
            BATCH_BYTES,
            |record| {
                let index = report.records;
                report.records += 1;
                let actual = stream.next()?;
                // A BAM archive holds the transformed records byte for byte.
                if let Some(ArchiveRecord::Encoded(actual)) = &actual
                    && actual.as_ref() == record.as_ref()
                {
                    return Ok(());
                }
                let expected = RecordBuf::try_from_alignment_record(&original_header, record)
                    .context("decoding an original record")?;
                let name = record_name(&expected);
                match actual {
                    Some(actual) => {
                        if let Some(what) = difference(&expected, &actual.decode(header)?) {
                            report.record(index, name, what);
                        }
                    }
                    None => report.record(index, name, "missing from the archive"),
                }
                Ok(())
            },
        );
        let extra = match &compared {
            Ok(_) => report_extra(&mut stream, header, &mut report),
            Err(_) => Ok(()),
        };
        // Dropping the receiver stops the decoder if it is still sending.
        drop(stream);
        let decoded = decoder
            .join()
            .map_err(|_| anyhow::anyhow!("the archive decoder panicked"))?;
        compared?;
        extra?;
        decoded
    })?;
    Ok(report)
}

/// The reference sequences' names and lengths, which a CRAM writer keeps (it may add `M5`).
fn names_and_lengths(header: &sam::Header) -> Vec<(&[u8], usize)> {
    header
        .reference_sequences()
        .iter()
        .map(|(name, rs)| (name.as_ref(), rs.length().get()))
        .collect()
}

/// The original's header lines that the archive lacks or changed. lrsqueeze adds `@PG` and `@CO`
/// lines, and a CRAM writer adds `M5`/`UR` to `@SQ` lines (whose names and lengths are checked
/// separately), so only the original's own lines are looked for.
fn header_differences(original: &sam::Header, archive: &sam::Header) -> Vec<String> {
    let mut differences = Vec::new();
    if original.header() != archive.header() {
        differences.push("@HD".to_owned());
    }
    for (id, read_group) in original.read_groups() {
        if archive.read_groups().get(id) != Some(read_group) {
            differences.push(format!("@RG ID:{id}"));
        }
    }
    for (id, program) in original.programs().as_ref() {
        if archive.programs().as_ref().get(id) != Some(program) {
            differences.push(format!("@PG ID:{id}"));
        }
    }
    for comment in original.comments() {
        if !archive.comments().contains(comment) {
            differences.push(format!("@CO {}", describe_text(&comment.to_string())));
        }
    }
    differences
}

/// The policy recorded in an archive's header by `lrsqueeze pack`, if any: the last one, since
/// packing an archive again adds its policy after the earlier one.
fn recorded_policy(header: &sam::Header) -> anyhow::Result<Option<Policy>> {
    header
        .comments()
        .iter()
        .rev()
        .find_map(|c| Policy::from_header_line(&String::from_utf8_lossy(c)))
        .transpose()
        .context("parsing the policy recorded in the archive")
}

/// An archive record as read. BAM records stay encoded, to be compared byte for byte first.
enum ArchiveRecord {
    Encoded(bam::Record),
    Decoded(RecordBuf),
}

impl ArchiveRecord {
    fn decode(self, header: &sam::Header) -> anyhow::Result<RecordBuf> {
        match self {
            Self::Encoded(record) => RecordBuf::try_from_alignment_record(header, &record)
                .context("decoding an archive record"),
            Self::Decoded(record) => Ok(record),
        }
    }
}

/// A batch of archive records, or why reading stopped.
type ArchiveBatch = anyhow::Result<Vec<ArchiveRecord>>;

/// Reads the archive's records and sends them in batches, until the end or until the receiver
/// is gone.
fn decode(
    archive: &mut Archive,
    header: &sam::Header,
    threads: NonZero<usize>,
    batches: &SyncSender<ArchiveBatch>,
) -> anyhow::Result<()> {
    let records: Box<dyn Iterator<Item = anyhow::Result<ArchiveRecord>> + '_> = match archive {
        Archive::Cram(reader) => Box::new(reader.multithreaded_records(header, threads).map(|r| {
            r.map(ArchiveRecord::Decoded)
                .context("reading an archive record")
        })),
        Archive::Bam(reader) => Box::new(reader.records().map(|r| {
            r.map(ArchiveRecord::Encoded)
                .context("reading an archive record")
        })),
    };
    let mut batch = Vec::with_capacity(DECODED_BATCH);
    for record in records {
        match record {
            Ok(record) => batch.push(record),
            Err(e) => {
                // Pass the error on; the comparison reports it.
                let _ = batches.send(Err(e));
                return Ok(());
            }
        }
        if batch.len() == DECODED_BATCH {
            let full = mem::replace(&mut batch, Vec::with_capacity(DECODED_BATCH));
            if batches.send(Ok(full)).is_err() {
                return Ok(());
            }
        }
    }
    if !batch.is_empty() {
        // A failed send means the comparison stopped; it reports why.
        let _ = batches.send(Ok(batch));
    }
    Ok(())
}

/// The archive's records, in order.
struct Decoded {
    batches: Receiver<ArchiveBatch>,
    ready: VecDeque<ArchiveRecord>,
}

impl Decoded {
    fn next(&mut self) -> anyhow::Result<Option<ArchiveRecord>> {
        while self.ready.is_empty() {
            match self.batches.recv() {
                Ok(batch) => self.ready.extend(batch?),
                Err(_) => return Ok(None),
            }
        }
        Ok(self.ready.pop_front())
    }
}

/// Counts the archive records left after the original ran out, one at a time.
fn report_extra(
    stream: &mut Decoded,
    header: &sam::Header,
    report: &mut VerifyReport,
) -> anyhow::Result<()> {
    let mut index = report.records;
    while let Some(actual) = stream.next()? {
        let name = record_name(&actual.decode(header)?);
        report.record(index, name, "extra in the archive");
        index += 1;
    }
    Ok(())
}

fn record_name(record: &RecordBuf) -> String {
    record.name().map(ToString::to_string).unwrap_or_default()
}

/// The first field in which `actual` is not `expected`, allowing what CRAM changes.
fn difference(expected: &RecordBuf, actual: &RecordBuf) -> Option<String> {
    macro_rules! check {
        ($field:ident, $label:literal) => {
            if expected.$field() != actual.$field() {
                return Some(format!(
                    "{}: expected {:?}, archive has {:?}",
                    $label,
                    expected.$field(),
                    actual.$field()
                ));
            }
        };
    }
    check!(name, "name");
    check!(flags, "flags");
    check!(reference_sequence_id, "reference sequence");
    check!(alignment_start, "position");
    check!(mapping_quality, "mapping quality");
    if cigar_as_m(expected) != cigar_as_m(actual) {
        return Some("CIGAR".to_owned());
    }
    check!(mate_reference_sequence_id, "mate reference sequence");
    check!(mate_alignment_start, "mate position");
    check!(template_length, "template length");
    if expected.sequence() != actual.sequence() {
        return Some("sequence".to_owned());
    }
    if expected.quality_scores() != actual.quality_scores() {
        return Some("quality scores".to_owned());
    }
    let (want, have) = (
        sorted_fields(expected, None),
        sorted_fields(actual, Some(expected)),
    );
    let tags: BTreeSet<Tag> = want.iter().chain(&have).map(|(tag, _)| *tag).collect();
    let (tag, a, b) = tags.into_iter().find_map(|tag| {
        let (a, b) = (value(&want, tag), value(&have, tag));
        (!same_value(a, b)).then_some((tag, a, b))
    })?;
    let [t0, t1] = <[u8; 2]>::from(tag);
    Some(format!(
        "tag {}{}: expected {}, archive has {}",
        char::from(t0),
        char::from(t1),
        describe(a),
        describe(b)
    ))
}

fn value<'v>(fields: &[(Tag, &'v Value)], tag: Tag) -> Option<&'v Value> {
    fields.iter().find(|(t, _)| *t == tag).map(|(_, v)| *v)
}

/// Whether two tag values are equal, taking integers of any width as SAM does (a BAM writer may
/// store `NM:i:209` in one byte, a CRAM decoder in four), and floats as stored.
fn same_value(a: Option<&Value>, b: Option<&Value>) -> bool {
    let same_floats = |x: &[f32], y: &[f32]| {
        x.len() == y.len() && x.iter().zip(y).all(|(x, y)| x.to_bits() == y.to_bits())
    };
    match (a, b) {
        (Some(a), Some(b)) => match (a, b) {
            // By bits, so that a stored NaN matches itself.
            (Value::Float(x), Value::Float(y)) => x.to_bits() == y.to_bits(),
            (Value::Array(Array::Float(x)), Value::Array(Array::Float(y))) => same_floats(x, y),
            _ => match (a.as_int(), b.as_int()) {
                (Some(x), Some(y)) => x == y,
                _ => a == b,
            },
        },
        (a, b) => a == b,
    }
}

fn describe(value: Option<&Value>) -> String {
    value.map_or_else(|| "none".to_owned(), |v| describe_text(&format!("{v:?}")))
}

/// `text`, cut to 60 characters.
fn describe_text(text: &str) -> String {
    match text.char_indices().nth(60) {
        Some((end, _)) => format!("{}…", &text[..end]),
        None => text.to_owned(),
    }
}

/// The CIGAR with `=`/`X` as `M` and adjacent operations of one kind merged, as CRAM decodes it.
fn cigar_as_m(record: &RecordBuf) -> Vec<(Kind, usize)> {
    let mut ops: Vec<(Kind, usize)> = Vec::new();
    for op in record.cigar().as_ref() {
        let kind = match op.kind() {
            Kind::SequenceMatch | Kind::SequenceMismatch => Kind::Match,
            kind => kind,
        };
        match ops.last_mut() {
            Some((last, len)) if *last == kind => *len += op.len(),
            _ => ops.push((kind, op.len())),
        }
    }
    ops
}

/// The record's tags sorted by tag. With `original`, `MD`/`NM` are left out unless the original
/// has them, since decoding CRAM adds them.
fn sorted_fields<'r>(record: &'r RecordBuf, original: Option<&RecordBuf>) -> Vec<(Tag, &'r Value)> {
    let added_by_decoding = |tag: Tag| {
        original.is_some_and(|o| {
            matches!(&<[u8; 2]>::from(tag), b"MD" | b"NM") && o.data().get(&tag).is_none()
        })
    };
    let mut fields: Vec<_> = record
        .data()
        .iter()
        .filter(|(tag, _)| !added_by_decoding(*tag))
        .collect();
    fields.sort_by_key(|(tag, _)| *tag);
    fields
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integers_compare_across_widths() {
        let (byte, int) = (Value::UInt8(209), Value::Int32(209));
        assert!(same_value(Some(&byte), Some(&int)));
        assert!(!same_value(Some(&byte), Some(&Value::Int32(210))));
        assert!(!same_value(Some(&byte), None));
        assert!(!same_value(
            Some(&Value::from("209")),
            Some(&Value::Int32(209))
        ));
    }

    #[test]
    fn stored_nans_match_themselves() {
        let nan = Value::Float(f32::NAN);
        assert!(same_value(Some(&nan), Some(&nan)));
        let array = Value::Array(Array::Float(vec![1.0, f32::NAN]));
        assert!(same_value(Some(&array), Some(&array)));
        assert!(!same_value(Some(&nan), Some(&Value::Float(0.0))));
    }

    fn header(text: &str) -> sam::Header {
        sam::io::Reader::new(text.as_bytes())
            .read_header()
            .expect("valid header")
    }

    #[test]
    fn recorded_policy_is_the_last_one() {
        let first = Policy::preset(crate::policy::Platform::Ont, crate::policy::Mode::Safe);
        let last = Policy::preset(crate::policy::Platform::Ont, crate::policy::Mode::Max);
        let text = format!(
            "@HD\tVN:1.6\n@CO\t{}\n@CO\t{}\n",
            first.to_header_line(),
            last.to_header_line()
        );
        assert_eq!(recorded_policy(&header(&text)).expect("valid"), Some(last));
    }

    #[test]
    fn header_lines_must_survive() {
        let original = header(
            "@HD\tVN:1.6\tSO:coordinate\n@RG\tID:rg1\tSM:HG002\n@PG\tID:minimap2\tPN:minimap2\n\
             @CO\tbasecalled\n",
        );
        let kept = header(
            "@HD\tVN:1.6\tSO:coordinate\n@RG\tID:rg1\tSM:HG002\n@PG\tID:minimap2\tPN:minimap2\n\
             @PG\tID:lrsqueeze\tPN:lrsqueeze\tPP:minimap2\n@CO\tbasecalled\n@CO\tlrsqueeze lossy policy: x\n",
        );
        assert!(header_differences(&original, &kept).is_empty());
        let lost = header("@HD\tVN:1.6\tSO:coordinate\n@RG\tID:rg1\tSM:other\n");
        assert_eq!(
            header_differences(&original, &lost),
            ["@RG ID:rg1", "@PG ID:minimap2", "@CO basecalled"]
        );
    }
}
