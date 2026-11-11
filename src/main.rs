use std::{
    env, fs,
    io::{self, Write},
    num::NonZero,
    path::{Path, PathBuf},
    process::ExitCode,
};

use anyhow::Context;
use clap::{Parser, Subcommand, ValueEnum};
use lrsqueeze::{
    pack::{self, Format, IndexOutcome, Input, PackOptions},
    policy::{Mode, ModificationCalls, Overrides, Platform, Policy},
    verify::{self, VerifyOptions},
};

/// Lossy archives for ONT and PacBio long reads that pipelines read like the original.
#[derive(Debug, Parser)]
#[command(version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Pack an aligned BAM or CRAM into a lossy CRAM (or BAM) archive, or a FASTQ into a lossy
    /// FASTQ.gz.
    ///
    /// By default the changes are ones measured as negligible for variant callers, methylation
    /// pileups, and read filters, written as CRAM 3.0 that any htslib-based tool reads. `--max`
    /// trades a small measured accuracy cost for the smallest files.
    Pack(PackArgs),
    /// Check that an archive holds exactly what its policy makes of the original.
    ///
    /// Run it before deleting the original. Every record is compared: names, flags, alignments,
    /// sequence, qualities, and tags must be what the policy makes of the original, allowing only
    /// what CRAM itself changes (`=`/`X` CIGAR operations read back as `M`, and `MD`/`NM` added
    /// on decoding). Exits with an error if any record differs.
    Verify(VerifyArgs),
    /// Print the built-in policy for a platform as TOML, to edit and pass with `pack --config`.
    Policy {
        /// The platform.
        #[arg(long, value_enum)]
        platform: PlatformArg,
        /// The --max policy instead of the default.
        #[arg(long)]
        max: bool,
    },
}

#[derive(Debug, clap::Args)]
struct PackArgs {
    /// Aligned BAM or CRAM input, or FASTQ (.fastq or .fq, optionally .gz); `-` for standard
    /// input.
    input: String,
    /// Output archive: CRAM, BAM when the name ends in .bam, or FASTQ.gz for FASTQ input.
    #[arg(short, long)]
    output: PathBuf,
    /// Indexed FASTA the reads are aligned to (needed for CRAM input and output).
    #[arg(short = 'T', long)]
    reference: Option<PathBuf>,
    /// Smallest files: more aggressive changes and CRAM 3.1, which needs htslib 1.12+ to read.
    #[arg(long)]
    max: bool,
    /// The platform: needed for FASTQ, and for a BAM whose @RG PL does not say ONT or PACBIO.
    #[arg(long, value_enum)]
    platform: Option<PlatformArg>,
    /// Output format, if not chosen by the output name.
    #[arg(short = 'O', long, value_enum)]
    output_format: Option<FormatArg>,
    /// Base-modification calls to keep: all (default); cpg (5mC/5hmC only at CpG, others such
    /// as 6mA kept); cpg-only (CpG 5mC/5hmC and nothing else); none.
    #[arg(long, value_name = "CALLS")]
    methylation: Option<ModificationCalls>,
    /// Keep base-modification probabilities (`ML`) exactly instead of in 16 levels, for
    /// methylation calls at a fixed threshold that is not a multiple of 1/16. About 7% larger on
    /// PacBio Revio.
    #[arg(long)]
    keep_ml: bool,
    /// A custom policy (TOML; start from `lrsqueeze policy`), e.g. to change the quality bins.
    #[arg(short, long)]
    config: Option<PathBuf>,
    /// Threads: the size of each worker pool (BAM decompression, record rewriting, and
    /// compression; CRAM input is decoded on half as many). The pools mostly take turns, so CPU
    /// use stays within about this many cores.
    #[arg(short = '@', long, default_value_t = NonZero::new(4).expect("4 is nonzero"))]
    threads: NonZero<usize>,
    /// Do not index the output.
    #[arg(long)]
    no_index: bool,
}

#[derive(Debug, clap::Args)]
struct VerifyArgs {
    /// The BAM, CRAM, or FASTQ that was packed.
    original: PathBuf,
    /// The archive: CRAM, BAM when the name ends in .bam, or FASTQ.gz.
    archive: PathBuf,
    /// Indexed FASTA the reads are aligned to (needed for CRAM and for ONT archives).
    #[arg(short = 'T', long)]
    reference: Option<PathBuf>,
    /// The archive's format, if not chosen by its name.
    #[arg(short = 'O', long, value_enum)]
    archive_format: Option<FormatArg>,
    /// The policy the archive was made with (TOML), for archives that do not record it.
    #[arg(short, long, conflicts_with = "platform")]
    config: Option<PathBuf>,
    /// For archives that do not record their policy (FASTQ), the --platform they were packed
    /// with.
    #[arg(long, value_enum)]
    platform: Option<PlatformArg>,
    /// With --platform: the archive was packed with --max.
    #[arg(long, requires = "platform")]
    max: bool,
    /// With --platform: the --methylation the archive was packed with.
    #[arg(long, value_name = "CALLS", requires = "platform")]
    methylation: Option<ModificationCalls>,
    /// With --platform: the archive was packed with --keep-ml.
    #[arg(long, requires = "platform")]
    keep_ml: bool,
    /// Threads for each worker pool; two CRAM files are decoded on half as many each.
    #[arg(short = '@', long, default_value_t = NonZero::new(4).expect("4 is nonzero"))]
    threads: NonZero<usize>,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum PlatformArg {
    Ont,
    Pacbio,
}

impl From<PlatformArg> for Platform {
    fn from(arg: PlatformArg) -> Self {
        match arg {
            PlatformArg::Ont => Self::Ont,
            PlatformArg::Pacbio => Self::Pacbio,
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum FormatArg {
    Cram,
    Bam,
    Fastq,
}

fn mode(max: bool) -> Mode {
    if max { Mode::Max } else { Mode::Safe }
}

fn read_policy(path: &PathBuf) -> anyhow::Result<Policy> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    toml::from_str(&text).with_context(|| format!("parsing {}", path.display()))
}

fn format(arg: Option<FormatArg>, path: &Path) -> Format {
    match arg {
        Some(FormatArg::Cram) => Format::Cram,
        Some(FormatArg::Bam) => Format::Bam,
        Some(FormatArg::Fastq) => Format::Fastq,
        None => Format::from_path(path),
    }
}

/// Warns when `path` is BGZF or CRAM without its end-of-file marker, which a copy cut short
/// leaves.
fn warn_if_truncated(path: &Path, done: &str) -> anyhow::Result<()> {
    if pack::missing_eof_marker(path).with_context(|| format!("reading {}", path.display()))? {
        eprintln!(
            "lrsqueeze: warning: {} does not end with an end-of-file marker, so it may have been \
             cut short; its records are {done} as they are. Check that it is complete before \
             deleting it",
            path.display()
        );
    }
    Ok(())
}

fn run_pack(args: PackArgs) -> anyhow::Result<()> {
    let custom_policy = args.config.as_ref().map(read_policy).transpose()?;
    let format = format(args.output_format, &args.output);
    let options = PackOptions {
        input: if args.input == "-" {
            Input::Stdin
        } else {
            Input::Path(PathBuf::from(&args.input))
        },
        format,
        output: args.output,
        reference: args.reference,
        mode: mode(args.max),
        platform: args.platform.map(Platform::from),
        custom_policy,
        methylation: args.methylation,
        keep_ml: args.keep_ml,
        threads: args.threads,
        write_index: !args.no_index,
        command_line: env::args().collect::<Vec<_>>().join(" "),
    };

    if let Input::Path(path) = &options.input {
        warn_if_truncated(path, "packed")?;
    }
    let report = pack::pack(&options)?;
    let s = report.stats;
    let pct = |n: u64, d: u64| {
        if d == 0 {
            0.0
        } else {
            100.0 * n as f64 / d as f64
        }
    };
    let source = match report.platform {
        Some(platform) if args.max => format!("{platform} --max"),
        Some(platform) => platform.to_string(),
        None => "custom".to_owned(),
    };
    eprintln!("lrsqueeze: policy ({source}): {}", report.policy.summary());
    eprintln!(
        "lrsqueeze: {} records in, {} written ({} secondary dropped, {} with MM/ML restricted)",
        s.records_in, s.records_out, s.secondary_dropped, s.modifications_rewritten
    );
    eprintln!(
        "lrsqueeze: {} quality scores, {:.2}% changed; {} ML values quantized",
        s.quality_bases,
        pct(s.quality_changed, s.quality_bases),
        s.ml_values
    );
    if s.records_unmapped * 2 > s.records_out {
        let mut costs = Vec::new();
        if report.policy.reference_homopolymer_length().is_some()
            && !report.policy.protects_unmapped_homopolymers()
        {
            costs.push("get no homopolymer protection");
        }
        if format == Format::Cram {
            costs.push("have their bases stored in full instead of against the reference");
        }
        if !costs.is_empty() {
            eprintln!(
                "lrsqueeze: warning: {:.0}% of the records are unmapped, and unmapped reads {}. \
                 Align the reads first (e.g. with minimap2 or dorado aligner)",
                pct(s.records_unmapped, s.records_out),
                costs.join(" and ")
            );
        }
    }
    match report.index {
        IndexOutcome::SkippedUnsorted => eprintln!(
            "lrsqueeze: warning: records are not in coordinate order (reference sequences in header \
             order, positions ascending), so no index was written; sort the input to index it"
        ),
        IndexOutcome::SkippedUndeclaredOrder => eprintln!(
            "lrsqueeze: warning: the header does not declare SO:coordinate, which a BAM index needs, \
             so no index was written; sort the input (samtools sort) to index it"
        ),
        IndexOutcome::Written(_) | IndexOutcome::NotRequested => {}
    }
    Ok(())
}

fn run_verify(args: VerifyArgs) -> anyhow::Result<()> {
    let options = VerifyOptions {
        format: format(args.archive_format, &args.archive),
        original: args.original,
        archive: args.archive,
        reference: args.reference,
        policy: match (&args.config, args.platform) {
            (Some(path), _) => Some(read_policy(path)?),
            (None, Some(platform)) => Some(
                Policy::preset(platform.into(), mode(args.max)).with_overrides(Overrides {
                    methylation: args.methylation,
                    keep_ml: args.keep_ml,
                }),
            ),
            (None, None) => None,
        },
        threads: args.threads,
    };
    warn_if_truncated(&options.original, "checked")?;
    let report = verify::verify(&options)?;
    eprintln!("lrsqueeze: policy: {}", report.policy.summary());
    for line in &report.header_differences {
        eprintln!("lrsqueeze: header: the archive lacks or changed the original's {line}");
    }
    for d in &report.described {
        eprintln!("lrsqueeze: record {} ({}): {}", d.index + 1, d.name, d.what);
    }
    if !report.is_match() {
        anyhow::bail!(
            "{} header lines and {} of {} records differ from what the policy makes of the \
             original",
            report.header_differences.len(),
            report.differences,
            report.records
        );
    }
    eprintln!(
        "lrsqueeze: all {} records are what the policy makes of the original",
        report.records
    );
    Ok(())
}

fn print_policy(platform: Platform, max: bool) -> anyhow::Result<()> {
    let flag = if max { " --max" } else { "" };
    let text =
        toml::to_string(&Policy::preset(platform, mode(max))).context("serializing the policy")?;
    let written = write!(
        io::stdout().lock(),
        "# lrsqueeze built-in policy: {platform}{flag}\n{text}"
    );
    match written {
        // A reader that stops early (`| head`) is not an error.
        Err(e) if e.kind() != io::ErrorKind::BrokenPipe => Err(e).context("writing the policy"),
        _ => Ok(()),
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Command::Pack(args) => run_pack(args),
        Command::Verify(args) => run_verify(args),
        Command::Policy { platform, max } => print_policy(platform.into(), max),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("lrsqueeze: error: {e:#}");
            ExitCode::FAILURE
        }
    }
}
