//! End to end: SAM → BAM → `lrsqueeze::pack::pack` → CRAM or BAM, read back with pasta.

use std::{fs, io, num::NonZero, path::Path};

use lrsqueeze::{
    pack::{self, Format, IndexOutcome, Input, PackOptions},
    policy::{Mode, ModificationCalls, Overrides, Platform, Policy},
    verify::{self, VerifyOptions, VerifyReport},
};
use pasta::{
    bam, cram, fasta,
    sam::{
        self,
        alignment::{
            RecordBuf,
            record::data::field::Tag,
            record_buf::data::field::{Value, value::Array},
        },
    },
};

/// `ACGTACGTAC`, a 7 bp A run at 10..17, then 23 more bases: 40 bp.
const REFERENCE: &str = "ACGTACGTACAAAAAAACGTACGTACGTACGTACGTACGT";

/// Q13 (`.`) everywhere: the safe ONT bins store it as 12, the max ONT bins as 15 (except on
/// the homopolymer, where it stays 13).
const SAM: &str = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chrT\tLN:40
@RG\tID:rg1\tPL:ONT
r_run\t0\tchrT\t9\t60\t12M\t*\t0\t0\tACAAAAAAACGT\t............\tRG:Z:rg1\tMM:Z:C+m?,0;\tML:B:C,100
r_secondary\t256\tchrT\t20\t0\t4M\t*\t0\t0\t*\t*\tRG:Z:rg1\tMD:Z:4\tNM:i:0
r_unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t....\tRG:Z:rg1
";

fn write_inputs(dir: &Path, sam: &str) -> io::Result<()> {
    fs::write(dir.join("ref.fa"), format!(">chrT\n{REFERENCE}\n"))?;
    fs::write(dir.join("ref.fa.fai"), "chrT\t40\t6\t40\t41\n")?;

    let mut reader = sam::io::Reader::new(sam.as_bytes());
    let header = reader.read_header()?;
    let mut writer = bam::io::Writer::new(fs::File::create(dir.join("in.bam"))?);
    writer.write_header(&header)?;
    for result in reader.record_bufs(&header) {
        use sam::alignment::io::Write;
        writer.write_alignment_record(&header, &result?)?;
    }
    writer.try_finish()
}

fn options(dir: &Path, mode: Mode, format: Format) -> PackOptions {
    let output = match format {
        Format::Cram => dir.join("out.cram"),
        Format::Bam => dir.join("out.bam"),
        Format::Fastq => dir.join("out.fastq.gz"),
    };
    PackOptions {
        input: Input::Path(dir.join("in.bam")),
        output,
        format,
        reference: Some(dir.join("ref.fa")),
        mode,
        platform: None,
        custom_policy: None,
        methylation: None,
        keep_ml: false,
        threads: NonZero::new(2).expect("2 is nonzero"),
        write_index: true,
        command_line: "test".to_owned(),
    }
}

fn assert_policy_recorded(header: &sam::Header) {
    let comments = header.comments();
    assert!(
        comments
            .iter()
            .any(|c| c.starts_with(b"lrsqueeze lossy policy:")),
        "the policy is described in the header"
    );
    let recorded = comments
        .iter()
        .filter_map(|c| Policy::from_header_line(&String::from_utf8_lossy(c)))
        .next()
        .expect("the policy is recorded in the header");
    recorded.expect("the recorded policy parses");
}

fn read_cram(dir: &Path) -> io::Result<Vec<RecordBuf>> {
    let reader =
        fasta::io::indexed_reader::Builder::default().build_from_path(dir.join("ref.fa"))?;
    let repository =
        fasta::Repository::new(fasta::repository::adapters::IndexedReader::new(reader));
    let mut reader = cram::io::reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .build_from_path(dir.join("out.cram"))?;
    let header = reader.read_header()?;
    assert_policy_recorded(&header);
    reader
        .records(&header)
        .map(|r| RecordBuf::try_from_alignment_record(&header, &r?))
        .collect()
}

fn read_bam(dir: &Path) -> io::Result<Vec<RecordBuf>> {
    let mut reader = bam::io::reader::Builder.build_from_path(dir.join("out.bam"))?;
    let header = reader.read_header()?;
    assert_policy_recorded(&header);
    reader.record_bufs(&header).collect()
}

/// The CRAM major and minor version from the file definition.
fn cram_version(dir: &Path) -> io::Result<(u8, u8)> {
    let bytes = fs::read(dir.join("out.cram"))?;
    assert_eq!(&bytes[..4], b"CRAM");
    Ok((bytes[4], bytes[5]))
}

fn quals(record: &RecordBuf) -> Vec<u8> {
    record.quality_scores().as_ref().to_vec()
}

fn ml(record: &RecordBuf) -> Option<&Value> {
    record.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES)
}

fn names(records: &[RecordBuf]) -> Vec<String> {
    records
        .iter()
        .map(|r| r.name().map(|n| n.to_string()).unwrap_or_default())
        .collect()
}

#[test]
fn default_is_safe_cram_3_0_keeping_every_record() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), SAM)?;
    let report = pack::pack(&options(dir.path(), Mode::Safe, Format::Cram))?;
    assert_eq!(report.platform, Some(Platform::Ont), "detected from @RG PL");
    assert_eq!(report.stats.secondary_dropped, 0);
    assert!(matches!(report.index, IndexOutcome::Written(_)));
    assert_eq!(cram_version(dir.path())?, (3, 0));

    let records = read_cram(dir.path())?;
    assert_eq!(names(&records), ["r_run", "r_secondary", "r_unmapped"]);
    // 12M from reference 8: query 2..=8 lie on the A run (reference 10..17) and keep Q13; the
    // rest is in the 10-14 bin, stored as 12.
    let mut want = vec![12; 12];
    want[2..=8].fill(13);
    assert_eq!(quals(&records[0]), want);
    assert_eq!(
        ml(&records[0]),
        Some(&Value::Array(Array::UInt8(vec![104])))
    );
    assert_eq!(
        quals(&records[2]),
        vec![12; 4],
        "unmapped reads are never protected"
    );
    Ok(())
}

#[test]
fn max_is_cram_3_1_with_homopolymer_protection() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), SAM)?;
    let report = pack::pack(&options(dir.path(), Mode::Max, Format::Cram))?;
    assert_eq!(report.stats.secondary_dropped, 1);
    assert_eq!(cram_version(dir.path())?, (3, 1));

    let records = read_cram(dir.path())?;
    assert_eq!(names(&records), ["r_run", "r_unmapped"]);
    // 12M from reference 8: query 2..=8 lie on the A run (reference 10..17) and keep Q13.
    let mut want = vec![15; 12];
    want[2..=8].fill(13);
    assert_eq!(quals(&records[0]), want);
    assert_eq!(
        quals(&records[1]),
        vec![15; 4],
        "unmapped reads are never protected"
    );
    Ok(())
}

#[test]
fn bam_output_needs_no_reference_without_homopolymer_protection() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), SAM)?;
    let mut opts = options(dir.path(), Mode::Safe, Format::Bam);
    opts.reference = None;
    opts.platform = Some(Platform::Pacbio);
    // The header says ONT; a platform given explicitly takes precedence over detection.
    let report = pack::pack(&opts)?;
    assert_eq!(
        report.index,
        IndexOutcome::Written(dir.path().join("out.bam.bai"))
    );
    let records = read_bam(dir.path())?;
    assert_eq!(records.len(), 3);
    assert_eq!(
        quals(&records[0]),
        vec![10; 12],
        "Q13 in the PacBio 0-17 bin"
    );
    assert_eq!(
        ml(&records[0]),
        Some(&Value::Array(Array::UInt8(vec![104])))
    );
    Ok(())
}

#[test]
fn reference_is_required_where_it_is_used() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), SAM)?;
    // ONT protects homopolymers in both modes; CRAM always needs the reference.
    for (mode, format) in [
        (Mode::Safe, Format::Bam),
        (Mode::Max, Format::Bam),
        (Mode::Safe, Format::Cram),
    ] {
        let mut opts = options(dir.path(), mode, format);
        opts.reference = None;
        let err = pack::pack(&opts).expect_err("needs the reference");
        assert!(err.to_string().contains("-T"), "{mode:?} {format:?}: {err}");
    }
    Ok(())
}

#[test]
fn platform_comes_from_the_header_or_the_option() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), &SAM.replace("\tPL:ONT", ""))?;
    let err = pack::pack(&options(dir.path(), Mode::Safe, Format::Cram))
        .expect_err("no platform to detect");
    assert!(err.to_string().contains("--platform"), "{err}");

    let mut opts = options(dir.path(), Mode::Safe, Format::Cram);
    opts.platform = Some(Platform::Pacbio);
    let report = pack::pack(&opts)?;
    assert_eq!(report.platform, Some(Platform::Pacbio));
    let records = read_cram(dir.path())?;
    assert_eq!(
        quals(&records[0]),
        vec![10; 12],
        "Q13 in the PacBio 0-17 bin"
    );
    Ok(())
}

#[test]
fn custom_policy_replaces_the_built_in_one() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), SAM)?;
    let mut opts = options(dir.path(), Mode::Safe, Format::Cram);
    opts.custom_policy = Some(Policy::lossless());
    let report = pack::pack(&opts)?;
    assert_eq!(report.platform, None);
    let records = read_cram(dir.path())?;
    assert_eq!(records.len(), 3);
    assert_eq!(quals(&records[0]), vec![13; 12]);
    assert_eq!(
        ml(&records[0]),
        Some(&Value::Array(Array::UInt8(vec![100])))
    );
    Ok(())
}

#[test]
fn format_follows_the_extension() {
    assert_eq!(Format::from_path(Path::new("a/x.bam")), Format::Bam);
    assert_eq!(Format::from_path(Path::new("x.BAM")), Format::Bam);
    assert_eq!(Format::from_path(Path::new("x.cram")), Format::Cram);
    assert_eq!(Format::from_path(Path::new("x")), Format::Cram);
}

/// One read with Cs at query 1 (reference CA) and 9 (reference CG at 17), an implicit 5mC part
/// calling both, and a 6mA call on its third A (query 3).
const SAM_MODS: &str = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chrT\tLN:40
@RG\tID:rg1\tPL:PACBIO
r_mods\t0\tchrT\t9\t60\t12M\t*\t0\t0\tACAAAAAAACGT\t............\tRG:Z:rg1\tMM:Z:C+m.,0,0;A+a.,2;\tML:B:C,50,220,7
";

fn mods(record: &RecordBuf) -> (Option<String>, Option<Vec<u8>>) {
    let mm = match record.data().get(&Tag::BASE_MODIFICATIONS) {
        Some(Value::String(s)) => Some(s.to_string()),
        _ => None,
    };
    let ml = match ml(record) {
        Some(Value::Array(Array::UInt8(v))) => Some(v.clone()),
        _ => None,
    };
    (mm, ml)
}

#[test]
fn methylation_modes_keep_cpg_calls() -> anyhow::Result<()> {
    // 16-level ML: 220 -> 216, 7 -> 8.
    let cases = [
        (
            ModificationCalls::Cpg,
            Some("C+m?,1;A+a.,2;"),
            Some(vec![216, 8]),
        ),
        (ModificationCalls::CpgOnly, Some("C+m?,1;"), Some(vec![216])),
        (ModificationCalls::None, None, None),
    ];
    for (calls, want_mm, want_ml) in cases {
        let dir = tempfile::tempdir()?;
        write_inputs(dir.path(), SAM_MODS)?;
        let mut opts = options(dir.path(), Mode::Safe, Format::Cram);
        opts.methylation = Some(calls);
        let report = pack::pack(&opts)?;
        assert_eq!(report.stats.modifications_rewritten, 1, "{calls}");
        let records = read_cram(dir.path())?;
        assert_eq!(
            mods(&records[0]),
            (want_mm.map(str::to_owned), want_ml),
            "{calls}"
        );
    }
    Ok(())
}

#[test]
fn keep_ml_keeps_probabilities_exactly() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), SAM_MODS)?;
    let mut opts = options(dir.path(), Mode::Safe, Format::Cram);
    opts.keep_ml = true;
    opts.methylation = Some(ModificationCalls::Cpg);
    pack::pack(&opts)?;
    let records = read_cram(dir.path())?;
    // As in methylation_modes_keep_cpg_calls, without 220 -> 216 and 7 -> 8.
    assert_eq!(
        mods(&records[0]),
        (Some("C+m?,1;A+a.,2;".to_owned()), Some(vec![220, 7]))
    );
    Ok(())
}

#[test]
fn cpg_methylation_needs_the_reference() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), SAM_MODS)?;
    let mut opts = options(dir.path(), Mode::Safe, Format::Bam);
    opts.reference = None;
    opts.methylation = Some(ModificationCalls::CpgOnly);
    let err = pack::pack(&opts).expect_err("CpG context needs the reference");
    assert!(err.to_string().contains("-T"), "{err}");
    // Without it, PacBio BAM output with all calls kept works.
    opts.methylation = None;
    pack::pack(&opts)?;
    Ok(())
}

fn verify_archive(
    dir: &Path,
    format: Format,
    policy: Option<Policy>,
) -> anyhow::Result<VerifyReport> {
    let archive = match format {
        Format::Cram => "out.cram",
        Format::Bam => "out.bam",
        Format::Fastq => "out.fastq.gz",
    };
    verify::verify(&VerifyOptions {
        original: dir.join("in.bam"),
        archive: dir.join(archive),
        format,
        reference: Some(dir.join("ref.fa")),
        policy,
        threads: NonZero::new(2).expect("2 is nonzero"),
    })
}

#[test]
fn verify_accepts_every_mode_and_format() -> anyhow::Result<()> {
    for (mode, format) in [
        (Mode::Safe, Format::Cram),
        (Mode::Max, Format::Cram),
        (Mode::Safe, Format::Bam),
    ] {
        let dir = tempfile::tempdir()?;
        write_inputs(dir.path(), SAM)?;
        pack::pack(&options(dir.path(), mode, format))?;
        let report = verify_archive(dir.path(), format, None)?;
        assert!(
            report.is_match(),
            "{mode:?} {format:?}: {:?}",
            report.described
        );
        let written = if mode == Mode::Max { 2 } else { 3 };
        assert_eq!(report.records, written, "{mode:?} {format:?}");
    }
    Ok(())
}

#[test]
fn verify_finds_records_that_differ_or_are_missing() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), SAM)?;
    pack::pack(&options(dir.path(), Mode::Safe, Format::Cram))?;

    // A policy that differs from the one the archive records is refused, not checked.
    let error = verify_archive(dir.path(), Format::Cram, Some(Policy::lossless()))
        .expect_err("the archive records another policy");
    assert!(
        format!("{error:#}").contains("records a different policy"),
        "{error:#}"
    );

    // Unbinned qualities (a plain copy) differ from what the ONT policy makes of them.
    fs::copy(dir.path().join("in.bam"), dir.path().join("out.bam"))?;
    let report = verify_archive(
        dir.path(),
        Format::Bam,
        Some(Policy::preset(Platform::Ont, Mode::Safe)),
    )?;
    assert!(!report.is_match());
    assert_eq!(report.described[0].name, "r_run");
    assert_eq!(report.described[0].what, "quality scores");

    // An original with one more record than was packed.
    let extra = "r_extra\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t....\tRG:Z:rg1\n";
    write_inputs(dir.path(), &format!("{SAM}{extra}"))?;
    let report = verify_archive(dir.path(), Format::Cram, None)?;
    assert_eq!(report.differences, 1);
    assert_eq!(report.described[0].name, "r_extra");
    assert_eq!(report.described[0].what, "missing from the archive");
    Ok(())
}

#[test]
fn verify_needs_a_recorded_or_given_policy() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), SAM)?;
    // An archive another tool made: a plain BAM copy, with no lrsqueeze policy in its header.
    fs::copy(dir.path().join("in.bam"), dir.path().join("out.bam"))?;
    let error = verify_archive(dir.path(), Format::Bam, None).expect_err("no policy");
    assert!(
        format!("{error:#}").contains("does not record its policy"),
        "{error:#}"
    );
    let report = verify_archive(dir.path(), Format::Bam, Some(Policy::lossless()))?;
    assert!(report.is_match(), "{:?}", report.described);
    Ok(())
}

/// A read with `=`/`X` CIGAR operations and `NM`/`MD`: reference 20..28 is `ACGTACGT`, the read
/// has `A` for the `T` at 23.
const SAM_EQX: &str = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chrT\tLN:40
@RG\tID:rg1\tPL:ONT
r_run\t0\tchrT\t9\t60\t12M\t*\t0\t0\tACAAAAAAACGT\t............\tRG:Z:rg1
r_eqx\t0\tchrT\t21\t60\t3=1X4=\t*\t0\t0\tACGAACGT\t........\tRG:Z:rg1\tNM:i:1\tMD:Z:3T4
";

#[test]
fn verify_allows_what_cram_changes() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), SAM_EQX)?;
    pack::pack(&options(dir.path(), Mode::Safe, Format::Cram))?;
    let records = read_cram(dir.path())?;
    // CRAM gives the CIGAR back as M; verify must accept that and the decoded NM.
    let cigar: Vec<_> = records[1]
        .cigar()
        .as_ref()
        .iter()
        .map(|op| (op.kind(), op.len()))
        .collect();
    assert_eq!(
        cigar,
        [(pasta::sam::alignment::record::cigar::op::Kind::Match, 8)]
    );
    let report = verify_archive(dir.path(), Format::Cram, None)?;
    assert!(report.is_match(), "{:?}", report.described);
    assert_eq!(report.records, 2);
    Ok(())
}

#[test]
fn verify_finds_records_extra_in_the_archive() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), SAM_EQX)?;
    pack::pack(&options(dir.path(), Mode::Safe, Format::Cram))?;
    // An original with the last record left out.
    let fewer = SAM_EQX.lines().take(4).collect::<Vec<_>>().join("\n") + "\n";
    write_inputs(dir.path(), &fewer)?;
    let report = verify_archive(dir.path(), Format::Cram, None)?;
    assert_eq!((report.records, report.differences), (1, 1));
    assert_eq!(report.described[0].index, 1);
    assert_eq!(report.described[0].name, "r_eqx");
    assert_eq!(report.described[0].what, "extra in the archive");
    Ok(())
}

#[test]
fn bam_index_needs_a_declared_sort_order() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    // In order, but the header does not say so (as after `samtools reset`).
    write_inputs(dir.path(), &SAM.replace("\tSO:coordinate", ""))?;
    let report = pack::pack(&options(dir.path(), Mode::Safe, Format::Bam))?;
    assert_eq!(report.index, IndexOutcome::SkippedUndeclaredOrder);
    assert_eq!(read_bam(dir.path())?.len(), 3);
    Ok(())
}

/// Two ONT reads as `samtools fastq -T MM,ML` writes them; `r1` has a 6 bp A run.
const FASTQ: &str = "@r1\tMM:Z:C+m?,0;\tML:B:C,100
CAAAAAACGT
+
..........
@r2 runid=x
ACGTACGT
+
IIIIIIII
";

fn fastq_options(dir: &Path, platform: Option<Platform>) -> PackOptions {
    fs::write(dir.join("in.fastq"), FASTQ).expect("writable");
    PackOptions {
        input: Input::Path(dir.join("in.fastq")),
        reference: None,
        platform,
        ..options(dir, Mode::Safe, Format::Fastq)
    }
}

fn read_gz(path: &Path) -> anyhow::Result<String> {
    let mut text = String::new();
    io::Read::read_to_string(
        &mut flate2::read::MultiGzDecoder::new(fs::File::open(path)?),
        &mut text,
    )?;
    Ok(text)
}

#[test]
fn fastq_packs_to_fastq_gz_and_verifies() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    let report = pack::pack(&fastq_options(dir.path(), Some(Platform::Ont)))?;
    assert_eq!(report.index, IndexOutcome::NotRequested);
    // Q13 ('.') is binned to 12 ('-') outside the read's A run; Q40 ('I') stays 40; ML 100 -> 104.
    assert_eq!(
        read_gz(&dir.path().join("out.fastq.gz"))?,
        "@r1\tMM:Z:C+m?,0;\tML:B:C,104\nCAAAAAACGT\n+\n-......---\n@r2 runid=x\nACGTACGT\n+\nIIIIIIII\n"
    );
    let verify = |policy| {
        verify::verify(&VerifyOptions {
            original: dir.path().join("in.fastq"),
            archive: dir.path().join("out.fastq.gz"),
            format: Format::Fastq,
            reference: None,
            policy,
            threads: NonZero::new(2).expect("2 is nonzero"),
        })
    };
    let report = verify(Some(Policy::preset(Platform::Ont, Mode::Safe)))?;
    assert!(report.is_match(), "{:?}", report.described);
    assert_eq!(report.records, 2);
    let report = verify(Some(Policy::preset(Platform::Ont, Mode::Max)))?;
    assert_eq!(report.described[0].what, "quality scores");
    assert!(verify(None).is_err(), "a FASTQ does not record its policy");
    Ok(())
}

#[test]
fn fastq_needs_a_platform_and_fastq_output() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    let error = pack::pack(&fastq_options(dir.path(), None)).expect_err("no platform");
    assert!(format!("{error:#}").contains("--platform"), "{error:#}");
    let to_cram = PackOptions {
        format: Format::Cram,
        output: dir.path().join("out.cram"),
        ..fastq_options(dir.path(), Some(Platform::Ont))
    };
    let error = pack::pack(&to_cram).expect_err("FASTQ to CRAM");
    assert!(format!("{error:#}").contains(".fastq.gz"), "{error:#}");
    Ok(())
}

#[test]
fn fastq_packs_pacbio_reads() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    for (mode, quals) in [
        // 3 bins: Q13 ('.') to 10 ('+'), Q40 ('I') stays; no homopolymer protection.
        (Mode::Safe, ["++++++++++", "IIIIIIII"]),
        // One score per read, from its mean error.
        (Mode::Max, ["..........", "IIIIIIII"]),
    ] {
        let options = PackOptions {
            mode,
            ..fastq_options(dir.path(), Some(Platform::Pacbio))
        };
        pack::pack(&options)?;
        let text = read_gz(&options.output)?;
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!([lines[3], lines[7]], quals, "{mode:?}");
        let report = verify::verify(&VerifyOptions {
            original: dir.path().join("in.fastq"),
            archive: options.output.clone(),
            format: Format::Fastq,
            reference: None,
            policy: Some(Policy::preset(Platform::Pacbio, mode)),
            threads: NonZero::new(2).expect("2 is nonzero"),
        })?;
        assert!(report.is_match(), "{mode:?}: {:?}", report.described);
    }
    Ok(())
}

#[test]
fn fastq_verify_finds_missing_and_extra_records() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    pack::pack(&fastq_options(dir.path(), Some(Platform::Ont)))?;
    let verify = |original: &str| -> anyhow::Result<VerifyReport> {
        fs::write(dir.path().join("in.fastq"), original)?;
        verify::verify(&VerifyOptions {
            original: dir.path().join("in.fastq"),
            archive: dir.path().join("out.fastq.gz"),
            format: Format::Fastq,
            reference: None,
            policy: Some(Policy::preset(Platform::Ont, Mode::Safe)),
            threads: NonZero::new(2).expect("2 is nonzero"),
        })
    };
    let report = verify(&format!("{FASTQ}@r3\nAC\n+\nII\n"))?;
    assert_eq!((report.records, report.differences), (3, 1));
    assert_eq!(
        (report.described[0].index, report.described[0].name.as_str()),
        (2, "r3")
    );
    assert_eq!(report.described[0].what, "missing from the archive");

    let first: String = FASTQ.lines().take(4).map(|l| format!("{l}\n")).collect();
    let report = verify(&first)?;
    assert_eq!((report.records, report.differences), (1, 1));
    assert_eq!(
        (report.described[0].index, report.described[0].name.as_str()),
        (1, "r2")
    );
    assert_eq!(report.described[0].what, "extra in the archive");
    Ok(())
}

#[test]
fn fastq_keep_ml_keeps_probabilities_and_verifies() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    let options = PackOptions {
        keep_ml: true,
        ..fastq_options(dir.path(), Some(Platform::Ont))
    };
    pack::pack(&options)?;
    let text = read_gz(&options.output)?;
    // ML 100 stays 100 instead of becoming 104.
    assert!(
        text.starts_with("@r1\tMM:Z:C+m?,0;\tML:B:C,100\n"),
        "{text}"
    );
    let verify = |policy| {
        verify::verify(&VerifyOptions {
            original: dir.path().join("in.fastq"),
            archive: options.output.clone(),
            format: Format::Fastq,
            reference: None,
            policy: Some(policy),
            threads: NonZero::new(2).expect("2 is nonzero"),
        })
    };
    let keep_ml = Overrides {
        methylation: None,
        keep_ml: true,
    };
    let report = verify(Policy::preset(Platform::Ont, Mode::Safe).with_overrides(keep_ml))?;
    assert!(report.is_match(), "{:?}", report.described);
    let report = verify(Policy::preset(Platform::Ont, Mode::Safe))?;
    assert_eq!(report.described[0].what, "header line", "16 levels differ");
    Ok(())
}

/// The reference in `dir` as a repository.
fn repository(dir: &Path) -> io::Result<fasta::Repository> {
    let reader =
        fasta::io::indexed_reader::Builder::default().build_from_path(dir.join("ref.fa"))?;
    Ok(fasta::Repository::new(
        fasta::repository::adapters::IndexedReader::new(reader),
    ))
}

/// Writes `sam` as `in.cram` in `dir` (whose reference [`write_inputs`] wrote), unchanged.
fn write_cram_input(dir: &Path, sam: &str) -> anyhow::Result<()> {
    let mut reader = sam::io::Reader::new(sam.as_bytes());
    let header = reader.read_header()?;
    let mut writer = cram::io::writer::Builder::default()
        .set_reference_sequence_repository(repository(dir)?)
        .build_from_path(dir.join("in.cram"))?;
    writer.write_header(&header)?;
    for result in reader.record_bufs(&header) {
        use sam::alignment::io::Write;
        writer.write_alignment_record(&header, &result?)?;
    }
    writer.try_finish(&header)?;
    Ok(())
}

/// `records` with each one's tags in sorted order.
fn tags_sorted(mut records: Vec<RecordBuf>) -> Vec<RecordBuf> {
    for record in &mut records {
        let mut fields: Vec<(Tag, Value)> = record
            .data()
            .iter()
            .map(|(tag, value)| (tag, value.clone()))
            .collect();
        fields.sort_by_key(|(tag, _)| *tag.as_ref());
        *record.data_mut() = fields.into_iter().collect();
    }
    records
}

#[test]
fn cram_input_packs_and_verifies_like_bam() -> anyhow::Result<()> {
    for (mode, format) in [
        (Mode::Safe, Format::Cram),
        (Mode::Max, Format::Cram),
        (Mode::Safe, Format::Bam),
    ] {
        let (from_bam, from_cram) = (tempfile::tempdir()?, tempfile::tempdir()?);
        write_inputs(from_bam.path(), SAM)?;
        pack::pack(&options(from_bam.path(), mode, format))?;
        write_inputs(from_cram.path(), SAM)?;
        write_cram_input(from_cram.path(), SAM)?;
        let options = PackOptions {
            input: Input::Path(from_cram.path().join("in.cram")),
            ..options(from_cram.path(), mode, format)
        };
        let report = pack::pack(&options)?;
        assert_eq!(report.stats.records_in, 3, "{mode:?} {format:?}");

        let read = |dir: &Path| match format {
            Format::Cram => read_cram(dir),
            _ => read_bam(dir),
        };
        // Decoding CRAM puts `RG` after the other tags, as HTSlib does.
        assert_eq!(
            tags_sorted(read(from_cram.path())?),
            tags_sorted(read(from_bam.path())?),
            "{mode:?} {format:?}"
        );

        let archive = if format == Format::Cram {
            "out.cram"
        } else {
            "out.bam"
        };
        let report = verify::verify(&VerifyOptions {
            original: from_cram.path().join("in.cram"),
            archive: from_cram.path().join(archive),
            format,
            reference: Some(from_cram.path().join("ref.fa")),
            policy: None,
            threads: NonZero::new(2).expect("2 is nonzero"),
        })?;
        assert!(
            report.is_match(),
            "{mode:?} {format:?}: {:?}",
            report.described
        );
    }
    Ok(())
}

#[test]
fn cram_input_needs_the_reference() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), SAM)?;
    write_cram_input(dir.path(), SAM)?;
    let options = PackOptions {
        input: Input::Path(dir.path().join("in.cram")),
        reference: None,
        custom_policy: Some(Policy::lossless()),
        ..options(dir.path(), Mode::Safe, Format::Bam)
    };
    let error = pack::pack(&options).expect_err("CRAM is decoded against the reference");
    assert!(
        format!("{error:#}").contains("reading CRAM input needs the reference (-T)"),
        "{error:#}"
    );
    Ok(())
}

#[test]
fn a_cram_cut_short_is_found() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), SAM)?;
    write_cram_input(dir.path(), SAM)?;
    let path = dir.path().join("in.cram");
    assert!(!pack::missing_eof_marker(&path)?, "a complete CRAM");
    let cram = fs::read(&path)?;
    // Cut off the end-of-file container: the CRAM still reads to its last whole container.
    fs::write(&path, &cram[..cram.len() - 38])?;
    assert!(
        pack::missing_eof_marker(&path)?,
        "a CRAM cut short between containers"
    );
    Ok(())
}

#[test]
fn packing_onto_the_input_is_refused() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), SAM)?;
    let before = fs::read(dir.path().join("in.bam"))?;
    for output in ["in.bam", "./in.bam"] {
        let options = PackOptions {
            output: dir.path().join(output),
            ..options(dir.path(), Mode::Safe, Format::Bam)
        };
        let error = pack::pack(&options).expect_err("the output is the input");
        assert!(format!("{error:#}").contains("is the input"), "{error:#}");
    }
    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(dir.path().join("in.bam"), dir.path().join("link.bam"))?;
        let options = PackOptions {
            output: dir.path().join("link.bam"),
            ..options(dir.path(), Mode::Safe, Format::Bam)
        };
        pack::pack(&options).expect_err("the output links to the input");
    }
    assert_eq!(
        fs::read(dir.path().join("in.bam"))?,
        before,
        "the input is untouched"
    );
    Ok(())
}

#[test]
fn indexes_of_an_earlier_output_are_removed() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    write_inputs(dir.path(), SAM)?;
    for (format, indexes) in [
        (Format::Cram, &["out.cram.crai"][..]),
        (Format::Bam, &["out.bam.bai", "out.bam.csi"][..]),
    ] {
        for index in indexes {
            fs::write(dir.path().join(index), "an index of another file")?;
        }
        let options = PackOptions {
            write_index: false,
            ..options(dir.path(), Mode::Safe, format)
        };
        pack::pack(&options)?;
        for index in indexes {
            assert!(!dir.path().join(index).exists(), "{index} is removed");
        }
    }
    Ok(())
}
