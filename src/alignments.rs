//! Reads aligned records from BAM or CRAM, as encoded BAM records.

use std::{
    io::{self, BufReader, Read},
    mem,
    num::NonZero,
};

use anyhow::Context;
use pasta::{
    bam, bgzf,
    cram::{self, io::reader::MdNm},
    fasta,
    sam::{self, alignment::io::Write},
};

/// The start of every CRAM file (CRAM specification, section 6).
pub(crate) const CRAM_MAGIC: &[u8] = b"CRAM";

/// The input source.
pub(crate) type Source = Box<dyn Read + Send>;

/// The first `n` bytes of `source` (fewer if it ends sooner), and the whole of it again. A pipe
/// may deliver fewer bytes per read than a format check needs.
///
/// # Errors
///
/// Fails if `source` cannot be read.
pub(crate) fn peek(mut source: Source, n: usize) -> io::Result<(Vec<u8>, Source)> {
    let mut head = Vec::with_capacity(n);
    (&mut source).take(n as u64).read_to_end(&mut head)?;
    let again: Source = Box::new(io::Cursor::new(head.clone()).chain(source));
    Ok((head, again))
}

/// An input of aligned records: BAM, or CRAM decoded against the reference.
#[expect(
    clippy::large_enum_variant,
    reason = "one per run; boxing the BAM reader saves nothing"
)]
pub(crate) enum Alignments {
    Bam(bam::io::Reader<bgzf::io::MultithreadedReader<BufReader<Source>>>),
    Cram(Box<cram::io::Reader<BufReader<Source>>>),
}

impl Alignments {
    /// Opens BAM or CRAM, whichever `source` starts with, and reads its header. BAM is
    /// decompressed on `threads` workers; CRAM is decoded against `repository`, with `MD` and
    /// `NM` as stored, so records keep the tags they had.
    ///
    /// # Errors
    ///
    /// Fails if the input cannot be read, is CRAM without a reference, or has a malformed header.
    pub(crate) fn open(
        source: Source,
        threads: NonZero<usize>,
        repository: Option<&fasta::Repository>,
    ) -> anyhow::Result<(Self, sam::Header)> {
        let (head, source) = peek(source, CRAM_MAGIC.len())?;
        let source = BufReader::with_capacity(1 << 20, source);
        let is_cram = head == CRAM_MAGIC;
        if is_cram {
            let repository = repository.context("reading CRAM input needs the reference (-T)")?;
            let mut reader = cram::io::reader::Builder::default()
                .set_reference_sequence_repository(repository.clone())
                .set_md_nm(MdNm::AsStored)
                .build_from_reader(source);
            let header = reader.read_header().context("reading the CRAM header")?;
            Ok((Self::Cram(Box::new(reader)), header))
        } else {
            let mut reader = bam::io::Reader::from(
                bgzf::io::MultithreadedReader::with_worker_count(threads, source),
            );
            let header = reader.read_header().context("reading the BAM header")?;
            Ok((Self::Bam(reader), header))
        }
    }

    /// Whether the input is CRAM.
    pub(crate) fn is_cram(&self) -> bool {
        matches!(self, Self::Cram(_))
    }

    /// The records, each encoded as BAM (without its length prefix). CRAM slices are decoded on
    /// `threads` workers, and records encoded on the calling thread.
    pub(crate) fn records<'a>(
        &'a mut self,
        header: &'a sam::Header,
        threads: NonZero<usize>,
    ) -> Box<dyn Iterator<Item = anyhow::Result<Vec<u8>>> + Send + 'a> {
        match self {
            Self::Bam(reader) => {
                let mut record = bam::Record::default();
                Box::new(std::iter::from_fn(move || {
                    match reader.read_record(&mut record) {
                        Ok(0) => None,
                        Ok(_) => Some(Ok(Vec::from(mem::take(&mut record)))),
                        Err(e) => Some(Err(anyhow::Error::new(e).context("reading a BAM record"))),
                    }
                }))
            }
            Self::Cram(reader) => {
                let mut encoder = bam::io::Writer::from(Vec::new());
                Box::new(
                    reader
                        .multithreaded_records(header, threads)
                        .map(move |record| {
                            let record = record.context("reading a CRAM record")?;
                            encoder
                                .write_alignment_record(header, &record)
                                .context("encoding a CRAM record as BAM")?;
                            // The writer puts each record's length first.
                            let encoded = encoder.get_ref()[4..].to_vec();
                            encoder.get_mut().clear();
                            Ok(encoded)
                        }),
                )
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Gives out one byte per read, as a slow pipe may.
    struct Trickle(io::Cursor<Vec<u8>>);

    impl Read for Trickle {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let n = buf.len().min(1);
            self.0.read(&mut buf[..n])
        }
    }

    #[test]
    fn peek_reads_past_short_reads_and_gives_everything_back() -> io::Result<()> {
        let source = Box::new(Trickle(io::Cursor::new(b"CRAM\x03\x01rest".to_vec())));
        let (head, mut again) = peek(source, 4)?;
        assert_eq!(head, CRAM_MAGIC);
        let mut all = Vec::new();
        again.read_to_end(&mut all)?;
        assert_eq!(all, b"CRAM\x03\x01rest");

        let (head, _) = peek(Box::new(io::Cursor::new(b"CR".to_vec())), 4)?;
        assert_eq!(head, b"CR", "a source shorter than asked for");
        Ok(())
    }
}
