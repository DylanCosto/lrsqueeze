//! Applies a [`Policy`] to encoded BAM records.

use std::io;

use crate::{
    homopolymer::Runs,
    methylation,
    policy::{HomopolymerSource, ModificationCalls, Policy, Quality, ReadQualitySource, Secondary},
    quality::{accuracy_to_phred, mean_error_phred, mean_phred},
    raw::{self, Layout, MISSING_QUALITY},
};

/// What a record's transform may need from its reference sequence.
#[derive(Clone, Copy, Debug, Default)]
pub struct RecordContext<'a> {
    /// Homopolymer runs, when [`Transformer::homopolymer_length`] is set.
    pub runs: Option<&'a Runs>,
    /// The whole reference sequence, when [`Transformer::needs_reference`] is set.
    pub reference: Option<&'a [u8]>,
}

/// Whether a transformed record is written.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Outcome {
    /// Write the (possibly edited) record.
    Write,
    /// Leave the record out.
    Skip,
}

/// Counts of what the transformer saw and changed.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Stats {
    /// Records read.
    pub records_in: u64,
    /// Records written.
    pub records_out: u64,
    /// Secondary alignments left out.
    pub secondary_dropped: u64,
    /// Quality scores seen on written records.
    pub quality_bases: u64,
    /// Quality scores kept native by homopolymer protection.
    pub quality_preserved: u64,
    /// Quality scores whose value changed.
    pub quality_changed: u64,
    /// `ML` probabilities quantized.
    pub ml_values: u64,
    /// Records whose `MD`/`NM` were removed because they have no `SEQ`.
    pub md_stripped: u64,
    /// Records whose `MM`/`ML` were restricted or removed.
    pub modifications_rewritten: u64,
    /// Written records that are unmapped.
    pub records_unmapped: u64,
}

impl std::ops::AddAssign for Stats {
    fn add_assign(&mut self, other: Self) {
        // Destructured without `..`, so a new field fails to compile until it is summed here.
        let Self {
            records_in,
            records_out,
            secondary_dropped,
            quality_bases,
            quality_preserved,
            quality_changed,
            ml_values,
            md_stripped,
            modifications_rewritten,
            records_unmapped,
        } = other;
        self.records_in += records_in;
        self.records_out += records_out;
        self.secondary_dropped += secondary_dropped;
        self.quality_bases += quality_bases;
        self.quality_preserved += quality_preserved;
        self.quality_changed += quality_changed;
        self.ml_values += ml_values;
        self.md_stripped += md_stripped;
        self.modifications_rewritten += modifications_rewritten;
        self.records_unmapped += records_unmapped;
    }
}

/// How the quality of each base is rewritten (BAM records and FASTQ alike).
#[derive(Clone, Debug)]
pub(crate) enum QualityOp {
    Keep,
    Bins {
        /// Phred score to binned Phred score.
        lut: Box<[u8; 256]>,
        /// The shortest protected homopolymer, and where homopolymers are found.
        preserve: Option<(usize, HomopolymerSource)>,
    },
    PerRead {
        source: ReadQualitySource,
        max: u8,
    },
    Drop,
}

impl QualityOp {
    pub(crate) fn new(quality: &Quality) -> Self {
        match quality {
            Quality::Keep => Self::Keep,
            Quality::Bins {
                bins,
                preserve_homopolymers,
                homopolymers_in,
            } => Self::Bins {
                lut: Box::new(bins.lut()),
                preserve: preserve_homopolymers.map(|n| (n, *homopolymers_in)),
            },
            Quality::PerRead { source, max } => Self::PerRead {
                source: *source,
                max: *max,
            },
            Quality::Drop => Self::Drop,
        }
    }
}

/// Bins each Phred score in `quals` with `lut`, except where `keep` is set.
pub(crate) fn bin_qualities(quals: &mut [u8], keep: &[bool], lut: &[u8; 256], stats: &mut Stats) {
    let mut changed = 0u64;
    let mut preserved = 0u64;
    for (q, &k) in quals.iter_mut().zip(keep) {
        if k {
            preserved += 1;
        } else {
            let v = lut[usize::from(*q)];
            changed += u64::from(v != *q);
            *q = v;
        }
    }
    stats.quality_changed += changed;
    stats.quality_preserved += preserved;
}

/// A read's single Phred score from its qualities `quals` and its `rq` accuracy, if it has one,
/// capped at `max`. Under [`ReadQualitySource::Rq`], reads without `rq` take the Phred score of
/// their mean per-base error.
pub(crate) fn read_quality(
    source: ReadQualitySource,
    quals: &[u8],
    rq: Option<f32>,
    max: u8,
) -> u8 {
    match (source, rq) {
        (ReadQualitySource::Rq, Some(accuracy)) => accuracy_to_phred(accuracy, max),
        (ReadQualitySource::Rq | ReadQualitySource::MeanError, _) => mean_error_phred(quals, max),
        (ReadQualitySource::MeanPhred, _) => mean_phred(quals, max),
    }
}

/// Gives every base of a read the Phred score `q`.
pub(crate) fn fill_quality(quals: &mut [u8], q: u8, stats: &mut Stats) {
    stats.quality_changed += quals.iter().filter(|&&v| v != q).count() as u64;
    quals.fill(q);
}

/// Applies one policy to records, reusing its buffers between them.
#[derive(Clone, Debug)]
pub struct Transformer {
    quality: QualityOp,
    ml_lut: Option<Box<[u8; 256]>>,
    calls: ModificationCalls,
    drop_secondary: bool,
    strip_md_without_sequence: bool,
    keep: Vec<bool>,
    bases: Vec<u8>,
    stats: Stats,
}

/// A lookup table quantizing `ML` probabilities (0–255) to the midpoints of `levels` equal bins.
pub fn ml_lut(levels: u16) -> [u8; 256] {
    let levels = usize::from(levels.clamp(2, 256));
    let width = 256 / levels;
    std::array::from_fn(|v| {
        let bin = (v / width).min(levels - 1);
        (bin * width + width / 2).min(255) as u8
    })
}

impl Transformer {
    /// A transformer for `policy` (which should be [validated](Policy::validate)).
    pub fn new(policy: &Policy) -> Self {
        Self {
            quality: QualityOp::new(&policy.quality),
            ml_lut: policy.modifications.ml_levels.map(|n| Box::new(ml_lut(n))),
            calls: policy.modifications.calls,
            drop_secondary: policy.alignments.secondary == Secondary::Drop,
            strip_md_without_sequence: policy.alignments.strip_md_without_sequence,
            keep: Vec::new(),
            bases: Vec::new(),
            stats: Stats::default(),
        }
    }

    /// The homopolymer length whose reference runs [`apply`](Self::apply) needs, if any.
    pub fn homopolymer_length(&self) -> Option<usize> {
        match self.quality {
            QualityOp::Bins {
                preserve:
                    Some((n, HomopolymerSource::Reference | HomopolymerSource::ReferenceOrRead)),
                ..
            } => Some(n),
            _ => None,
        }
    }

    /// Whether [`apply`](Self::apply) needs the reference sequence (for CpG context).
    pub fn needs_reference(&self) -> bool {
        matches!(
            self.calls,
            ModificationCalls::Cpg | ModificationCalls::CpgOnly
        )
    }

    /// The counts so far.
    pub fn stats(&self) -> Stats {
        self.stats
    }

    /// Transforms the encoded record in `buf`, with what `context` provides for its reference
    /// sequence.
    ///
    /// # Errors
    ///
    /// Fails on a malformed record or malformed `MM`/`ML`.
    pub fn apply(&mut self, buf: &mut Vec<u8>, context: RecordContext<'_>) -> io::Result<Outcome> {
        self.stats.records_in += 1;
        // Checks the record is whole before its fields are read.
        Layout::new(buf)?;
        if self.drop_secondary && raw::is_secondary(buf) {
            self.stats.secondary_dropped += 1;
            return Ok(Outcome::Skip);
        }
        if self.calls != ModificationCalls::All {
            self.rewrite_modifications(buf, context.reference)?;
        }
        let runs = context.runs;
        let layout = Layout::new(buf)?;
        let fields = {
            let data = layout.data(buf.len());
            raw::fields(&buf[data.clone()], data.start)?
        };

        self.rewrite_quality(buf, &layout, &fields, runs)?;

        if let Some(lut) = &self.ml_lut
            && let Some(ml) = fields
                .iter()
                .find(|f| f.tag == *b"ML" && f.subtype == Some(b'C'))
        {
            for v in &mut buf[ml.value.clone()] {
                *v = lut[usize::from(*v)];
            }
            self.stats.ml_values += ml.value.len() as u64;
        }

        if self.strip_md_without_sequence && layout.sequence_length() == 0 {
            let spans: Vec<_> = fields
                .iter()
                .filter(|f| matches!(&f.tag, b"MD" | b"NM"))
                .map(|f| f.span.clone())
                .collect();
            if !spans.is_empty() {
                raw::remove_fields(buf, &spans);
                self.stats.md_stripped += 1;
            }
        }

        self.stats.records_out += 1;
        self.stats.records_unmapped += u64::from(raw::is_unmapped(buf));
        Ok(Outcome::Write)
    }

    /// Restricts or removes `MM`/`ML` (and `MN`) as `self.calls` says. Records without `SEQ`
    /// or without both tags are left alone by the CpG modes, and so are records whose `MN` is not
    /// their `SEQ` length (such as a hard-clipped supplementary alignment that kept the primary's
    /// tags): their calls are not about these bases, and readers ignore them.
    fn rewrite_modifications(
        &mut self,
        buf: &mut Vec<u8>,
        reference: Option<&[u8]>,
    ) -> io::Result<()> {
        let layout = Layout::new(buf)?;
        let data = layout.data(buf.len());
        let fields = raw::fields(&buf[data.clone()], data.start)?;
        let find = |tag: &[u8; 2]| fields.iter().find(|f| f.tag == *tag);
        let span = |tag: &[u8; 2]| find(tag).map(|f| f.span.clone());
        let mut remove: Vec<_> = Vec::new();
        let mut replacement = None;
        match self.calls {
            ModificationCalls::All => return Ok(()),
            ModificationCalls::None => remove.extend(
                [span(b"MM"), span(b"ML"), span(b"MN")]
                    .into_iter()
                    .flatten(),
            ),
            ModificationCalls::Cpg | ModificationCalls::CpgOnly => {
                let (Some(mm), Some(ml)) = (find(b"MM"), find(b"ML")) else {
                    return Ok(());
                };
                if layout.sequence_length() == 0 || mm.kind != b'Z' || ml.subtype != Some(b'C') {
                    return Ok(());
                }
                if let Some(mn) = find(b"MN").and_then(|f| raw::int_value(buf, f))
                    && usize::try_from(mn).ok() != Some(layout.sequence_length())
                {
                    return Ok(());
                }
                let reverse = raw::is_reverse(buf);
                layout.bases_into(buf, &mut self.bases);
                let seq = methylation::original_orientation(&self.bases, reverse);
                let ref_cpg = match (reference, raw::position(buf)) {
                    (Some(reference), Some(start)) if !raw::is_unmapped(buf) => {
                        methylation::reference_cpg(
                            start,
                            layout.cigar(buf, &fields),
                            reference,
                            seq.len(),
                            reverse,
                        )
                    }
                    _ => vec![false; seq.len()],
                };
                let mm_text = std::str::from_utf8(&buf[mm.value.clone()])
                    .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "MM is not text"))?;
                let (new_mm, new_ml) = methylation::restrict_to_cpg(
                    mm_text,
                    &buf[ml.value.clone()],
                    &seq,
                    |p| seq.get(p + 1) == Some(&b'G') || ref_cpg[p],
                    self.calls == ModificationCalls::CpgOnly,
                )?;
                remove.extend([mm.span.clone(), ml.span.clone()]);
                if new_mm.is_empty() {
                    remove.extend(span(b"MN"));
                } else {
                    replacement = Some((new_mm, new_ml));
                }
            }
        }
        if remove.is_empty() {
            return Ok(());
        }
        remove.sort_by_key(|r| r.start);
        raw::remove_fields(buf, &remove);
        if let Some((mm, ml)) = replacement {
            raw::push_string(buf, *b"MM", &mm);
            raw::push_u8_array(buf, *b"ML", &ml);
        }
        self.stats.modifications_rewritten += 1;
        Ok(())
    }

    fn rewrite_quality(
        &mut self,
        buf: &mut [u8],
        layout: &Layout,
        fields: &[raw::Field],
        runs: Option<&Runs>,
    ) -> io::Result<()> {
        let range = layout.quality();
        if range.is_empty() || buf[range.start] == MISSING_QUALITY {
            return Ok(());
        }
        self.stats.quality_bases += range.len() as u64;
        match &self.quality {
            QualityOp::Keep => {}
            QualityOp::Bins { lut, preserve } => {
                let keep = &mut self.keep;
                keep.clear();
                keep.resize(range.len(), false);
                if let Some((n, source)) = preserve {
                    let unmapped = raw::is_unmapped(buf);
                    let in_read = match source {
                        HomopolymerSource::Reference => false,
                        HomopolymerSource::Read => true,
                        HomopolymerSource::ReferenceOrRead => unmapped,
                    };
                    if in_read {
                        layout.bases_into(buf, &mut self.bases);
                        Runs::find(&self.bases, *n).mark(keep);
                    } else if let (Some(runs), Some(start)) = (runs, raw::position(buf))
                        && !unmapped
                    {
                        let covered = raw::query_length(layout.cigar(buf, fields));
                        if covered != range.len() as u64 {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!(
                                    "the CIGAR covers {covered} read bases, but the record has {}",
                                    range.len()
                                ),
                            ));
                        }
                        runs.project(start, layout.cigar(buf, fields), keep);
                    }
                }
                bin_qualities(&mut buf[range], keep, lut, &mut self.stats);
            }
            QualityOp::PerRead { source, max } => {
                let rq = match source {
                    ReadQualitySource::Rq => fields
                        .iter()
                        .find(|f| f.tag == *b"rq")
                        .and_then(|f| raw::f32_value(buf, f)),
                    ReadQualitySource::MeanError | ReadQualitySource::MeanPhred => None,
                };
                let q = read_quality(*source, &buf[range.clone()], rq, *max);
                fill_quality(&mut buf[range], q, &mut self.stats);
            }
            QualityOp::Drop => {
                self.stats.quality_changed += range.len() as u64;
                buf[range].fill(MISSING_QUALITY);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        policy::{Mode, Modifications, Platform, Policy},
        raw::tests::{cg_field, record},
    };

    fn quals(buf: &[u8]) -> Vec<u8> {
        buf[Layout::new(buf).expect("valid").quality()].to_vec()
    }

    #[test]
    fn ml_lut_matches_prototype_midpoints() {
        let lut = ml_lut(16);
        assert_eq!(
            (lut[0], lut[15], lut[16], lut[204], lut[255]),
            (8, 8, 24, 200, 248)
        );
        let lut = ml_lut(10);
        assert_eq!(lut[255], 237, "the remainder joins the top bin");
    }

    #[test]
    fn ont_preset_bins_outside_homopolymers_and_quantizes_ml() {
        let mut t = Transformer::new(&Policy::preset(Platform::Ont, Mode::Max));
        let runs = Runs::find(b"ACGTAAAAACGT", 5); // run at 4..9
        let aux = [b"MLBC".as_slice(), &3u32.to_le_bytes(), &[0, 100, 255]].concat();
        // 6M starting at reference 2: reference 2..8, of which 4..8 are in the run.
        let mut buf = record(0, 2, &[(0, 6)], &[3, 12, 40, 41, 42, 25], &aux);
        assert_eq!(
            t.apply(
                &mut buf,
                RecordContext {
                    runs: Some(&runs),
                    reference: None
                }
            )
            .expect("valid"),
            Outcome::Write
        );
        assert_eq!(quals(&buf), vec![7, 15, 40, 41, 42, 25]);
        assert!(buf.ends_with(&[8, 104, 248]));
        let s = t.stats();
        assert_eq!(
            (s.quality_preserved, s.quality_changed, s.ml_values),
            (4, 2, 3)
        );
    }

    #[test]
    fn homopolymer_protection_follows_a_long_cigar_in_cg() {
        let mut t = Transformer::new(&Policy::preset(Platform::Ont, Mode::Max));
        let runs = Runs::find(b"ACGTAAAAACGT", 5); // run at 4..9
        let aux = cg_field(&[(0, 6)]);
        let mut buf = record(0, 2, &[(4, 6), (3, 6)], &[3, 12, 40, 41, 42, 25], &aux);
        t.apply(
            &mut buf,
            RecordContext {
                runs: Some(&runs),
                reference: None,
            },
        )
        .expect("valid");
        assert_eq!(quals(&buf), vec![7, 15, 40, 41, 42, 25]);
    }

    #[test]
    fn cpg_context_follows_a_long_cigar_in_cg() {
        let mut policy = Policy::lossless();
        policy.modifications.calls = ModificationCalls::Cpg;
        let mut t = Transformer::new(&policy);
        // Read AACAAA: its C is followed by A in the read, by G in the reference.
        let aux = [
            cg_field(&[(0, 6)]).as_slice(),
            b"MMZC+m?,0;\0",
            b"MLBC",
            &1u32.to_le_bytes(),
            &[200],
        ]
        .concat();
        let mut buf = record(0, 0, &[(4, 6), (3, 6)], &[30; 6], &aux);
        let sequence = Layout::new(&buf).expect("valid").sequence();
        buf[sequence].copy_from_slice(&[0x11, 0x21, 0x11]);
        t.apply(
            &mut buf,
            RecordContext {
                runs: None,
                reference: Some(b"AACGAA"),
            },
        )
        .expect("valid");
        let layout = Layout::new(&buf).expect("valid");
        let data = layout.data(buf.len());
        let fields = raw::fields(&buf[data.clone()], data.start).expect("valid");
        let ml = fields.iter().find(|f| f.tag == *b"ML").expect("call kept");
        assert_eq!(&buf[ml.value.clone()], &[200]);
    }

    #[test]
    fn malformed_records_are_errors() {
        let mut t = Transformer::new(&Policy::preset(Platform::Ont, Mode::Max));
        let runs = Runs::find(b"AAAAAAAAAA", 5);
        let context = RecordContext {
            runs: Some(&runs),
            reference: None,
        };
        // A CIGAR covering 6 bases on a 4-base read.
        let mut buf = record(0, 0, &[(0, 6)], &[30; 4], b"");
        let error = t
            .apply(&mut buf, context)
            .expect_err("CIGAR longer than SEQ");
        assert!(error.to_string().contains("covers 6 read bases"), "{error}");
        // Shorter than a record's fixed fields.
        let mut short = vec![0; 10];
        t.apply(&mut short, context).expect_err("too short");
    }

    #[test]
    fn cpg_modes_leave_calls_made_on_another_sequence() {
        let mut policy = Policy::lossless();
        policy.modifications.calls = ModificationCalls::CpgOnly;
        let mut t = Transformer::new(&policy);
        // The read's bases are AAA, with no C for the call to be on: MM was computed on the
        // 100-base read this record was clipped from.
        let calls = [
            b"MMZC+m?,5;\0".as_slice(),
            b"MLBC",
            &1u32.to_le_bytes(),
            &[200],
        ]
        .concat();
        let with_mn = [calls.as_slice(), b"MNC", &[100]].concat();
        let mut buf = record(0x800, 0, &[(5, 97), (0, 3)], &[30; 3], &with_mn);
        let before = buf.clone();
        t.apply(&mut buf, RecordContext::default())
            .expect("left alone");
        assert_eq!(buf, before);

        let mut buf = record(0x800, 0, &[(5, 97), (0, 3)], &[30; 3], &calls);
        t.apply(&mut buf, RecordContext::default())
            .expect_err("without MN, the call must be on this read");
    }

    #[test]
    fn read_homopolymers_need_no_alignment() {
        let mut policy = Policy::preset(Platform::Ont, Mode::Safe);
        let Quality::Bins {
            homopolymers_in, ..
        } = &mut policy.quality
        else {
            panic!("ONT bins qualities");
        };
        *homopolymers_in = HomopolymerSource::Read;
        policy.validate().expect("valid");
        let mut t = Transformer::new(&policy);
        assert_eq!(t.homopolymer_length(), None, "no reference runs are needed");
        // The test records' bases are all A: an unmapped 6-base read is one run of 6, and keeps
        // its native qualities; a 4-base read is too short a run, and is binned.
        let mut long = record(0x4, -1, &[], &[13; 6], b"");
        t.apply(&mut long, RecordContext::default()).expect("valid");
        assert_eq!(quals(&long), vec![13; 6]);
        let mut short = record(0x4, -1, &[], &[13; 4], b"");
        t.apply(&mut short, RecordContext::default())
            .expect("valid");
        assert_eq!(quals(&short), vec![12; 4]);
    }

    #[test]
    fn ont_presets_protect_unmapped_reads_by_their_own_bases() {
        let runs = Runs::find(b"ACGTAAAAACGT", 5); // run at 4..9
        let context = RecordContext {
            runs: Some(&runs),
            reference: None,
        };
        for mode in [Mode::Safe, Mode::Max] {
            let mut t = Transformer::new(&Policy::preset(Platform::Ont, mode));
            assert_eq!(
                t.homopolymer_length(),
                Some(5),
                "mapped reads use the reference"
            );
            // An unmapped read of six As is a run of six in its own bases: native qualities.
            let mut unmapped = record(0x4, -1, &[], &[13; 6], b"");
            t.apply(&mut unmapped, context).expect("valid");
            assert_eq!(quals(&unmapped), vec![13; 6], "{mode:?}");
            // A mapped read of the same bases at reference 0..6 overlaps the run at 4..6 only.
            let mut mapped = record(0, 0, &[(0, 6)], &[13; 6], b"");
            t.apply(&mut mapped, context).expect("valid");
            let binned = quals(&mapped)[0];
            assert_ne!(binned, 13, "{mode:?}");
            assert_eq!(quals(&mapped), vec![binned, binned, binned, binned, 13, 13]);
        }
    }

    #[test]
    fn drops_secondaries_and_strips_md_without_sequence() {
        let mut t = Transformer::new(&Policy::preset(Platform::Pacbio, Mode::Max));
        let mut secondary = record(0x100, 0, &[(0, 2)], &[30, 30], b"");
        assert_eq!(
            t.apply(&mut secondary, RecordContext::default())
                .expect("valid"),
            Outcome::Skip
        );

        let mut keep = Policy::preset(Platform::Pacbio, Mode::Safe);
        keep.alignments.secondary = Secondary::Keep;
        let mut t = Transformer::new(&keep);
        let aux = [
            b"NMi".as_slice(),
            &1i32.to_le_bytes(),
            b"MDZ1A0\0",
            b"rqf",
            &0.99f32.to_le_bytes(),
        ]
        .concat();
        let mut no_seq = record(0x100, 0, &[(0, 2)], &[], &aux);
        assert_eq!(
            t.apply(&mut no_seq, RecordContext::default())
                .expect("valid"),
            Outcome::Write
        );
        let layout = Layout::new(&no_seq).expect("valid");
        let data = layout.data(no_seq.len());
        let tags: Vec<_> = raw::fields(&no_seq[data.clone()], data.start)
            .expect("valid")
            .into_iter()
            .map(|f| f.tag)
            .collect();
        assert_eq!(tags, vec![*b"rq"]);
        assert_eq!(t.stats().md_stripped, 1);
    }

    #[test]
    fn per_read_quality_uses_rq_or_mean_error() {
        let mut policy = Policy::preset(Platform::Pacbio, Mode::Max);
        policy.quality = Quality::PerRead {
            source: ReadQualitySource::Rq,
            max: 40,
        };
        let mut t = Transformer::new(&policy);
        let aux = [b"rqf".as_slice(), &0.99f32.to_le_bytes()].concat();
        let mut buf = record(0, 0, &[(0, 3)], &[40, 3, 22], &aux);
        t.apply(&mut buf, RecordContext::default()).expect("valid");
        assert_eq!(quals(&buf), vec![20, 20, 20]);
        let mut buf = record(0, 0, &[(0, 2)], &[20, 20], b"");
        t.apply(&mut buf, RecordContext::default()).expect("valid");
        assert_eq!(quals(&buf), vec![20, 20], "no rq: mean per-base error");
    }

    #[test]
    fn pacbio_max_gives_each_read_its_mean_phred_score() {
        let mut t = Transformer::new(&Policy::preset(Platform::Pacbio, Mode::Max));
        // `rq` is ignored: the mean of 40, 3 and 23 is 22.
        let aux = [b"rqf".as_slice(), &0.99f32.to_le_bytes()].concat();
        let mut buf = record(0, 0, &[(0, 3)], &[40, 3, 23], &aux);
        t.apply(&mut buf, RecordContext::default()).expect("valid");
        assert_eq!(quals(&buf), vec![22, 22, 22]);
    }

    #[test]
    fn lossless_changes_nothing_and_missing_quality_stays_missing() {
        let policy = Policy::lossless();
        assert_eq!(policy.modifications, Modifications::default());
        let mut t = Transformer::new(&policy);
        let aux = [b"MLBC".as_slice(), &1u32.to_le_bytes(), &[77]].concat();
        let original = record(0x100, 0, &[(0, 2)], &[13, 37], &aux);
        let mut buf = original.clone();
        t.apply(&mut buf, RecordContext::default()).expect("valid");
        assert_eq!(buf, original);

        let mut t = Transformer::new(&Policy::preset(Platform::Ont, Mode::Max));
        let original = record(0, 0, &[(0, 2)], &[MISSING_QUALITY; 2], b"");
        let mut buf = original.clone();
        t.apply(&mut buf, RecordContext::default()).expect("valid");
        assert_eq!(buf, original);
    }
}
