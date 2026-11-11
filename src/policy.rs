//! What `lrsqueeze` changes in each record: built-in policies per platform and mode, or TOML.

use std::{fmt, str::FromStr};

use serde::{Deserialize, Serialize};

use crate::quality::{Bin, Bins};

/// The sequencing platform, which picks the built-in policy.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Platform {
    /// Oxford Nanopore (`@RG` `PL:ONT`).
    Ont,
    /// PacBio HiFi (`@RG` `PL:PACBIO`).
    Pacbio,
}

impl Platform {
    /// Both platforms.
    pub const ALL: [Self; 2] = [Self::Ont, Self::Pacbio];

    /// The platform's name on the command line.
    pub fn name(self) -> &'static str {
        match self {
            Self::Ont => "ont",
            Self::Pacbio => "pacbio",
        }
    }
}

impl fmt::Display for Platform {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

impl FromStr for Platform {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::ALL
            .into_iter()
            .find(|p| p.name().eq_ignore_ascii_case(s))
            .ok_or_else(|| format!("unknown platform {s:?} (expected ont or pacbio)"))
    }
}

/// How hard to compress: the default keeps archives nearly indistinguishable to pipelines.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum Mode {
    /// Changes that measured as negligible for callers and filters, in CRAM 3.0 that every
    /// htslib-based reader decodes.
    #[default]
    Safe,
    /// The smallest files with a small measured accuracy cost, in CRAM 3.1 (archive profile),
    /// which needs htslib 1.12 or newer (or pasta) to read.
    Max,
}

impl Policy {
    /// The built-in policy for a platform and mode (see the README for the evidence).
    pub fn preset(platform: Platform, mode: Mode) -> Self {
        let bins = |spec: &str| spec.parse::<Bins>().expect("built-in bins are valid");
        let modifications = Modifications {
            ml_levels: Some(16),
            calls: ModificationCalls::All,
        };
        let (quality, secondary) = match (platform, mode) {
            (Platform::Ont, Mode::Safe) => (
                Quality::Bins {
                    bins: bins(ONT_Q8),
                    preserve_homopolymers: Some(5),
                    homopolymers_in: HomopolymerSource::ReferenceOrRead,
                },
                Secondary::Keep,
            ),
            (Platform::Ont, Mode::Max) => (
                Quality::Bins {
                    bins: bins(ONT_Q4),
                    preserve_homopolymers: Some(5),
                    homopolymers_in: HomopolymerSource::ReferenceOrRead,
                },
                Secondary::Drop,
            ),
            (Platform::Pacbio, Mode::Safe) => (
                Quality::Bins {
                    bins: bins(PACBIO_Q3),
                    preserve_homopolymers: None,
                    homopolymers_in: HomopolymerSource::Reference,
                },
                Secondary::Keep,
            ),
            (Platform::Pacbio, Mode::Max) => (
                Quality::PerRead {
                    source: ReadQualitySource::MeanPhred,
                    max: 40,
                },
                Secondary::Drop,
            ),
        };
        Self {
            quality,
            modifications,
            alignments: Alignments {
                secondary,
                strip_md_without_sequence: true,
            },
        }
    }

    /// A policy that changes nothing, for baselines.
    pub fn lossless() -> Self {
        Self {
            quality: Quality::Keep,
            modifications: Modifications::default(),
            alignments: Alignments {
                secondary: Secondary::Keep,
                strip_md_without_sequence: false,
            },
        }
    }
}

/// ONT four bins with within-bin median representatives (HG002 HAC v5.2 calibration).
pub const ONT_Q4: &str = "0-9:7,10-19:15,20-29:27,30-93:38";
/// ONT eight bins with within-bin median representatives (HG002 HAC v5.2 calibration).
pub const ONT_Q8: &str = "0-4:4,5-9:7,10-14:12,15-19:17,20-24:23,25-29:28,30-34:32,35-93:40";
/// PacBio Revio three bins over its native seven (Q3/10/17/22/27/35/40).
pub const PACBIO_Q3: &str = "0-17:10,18-29:25,30-93:40";

/// Everything `lrsqueeze` changes in a record.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Policy {
    /// Per-base quality scores.
    pub quality: Quality,
    /// Base-modification probabilities.
    #[serde(default)]
    pub modifications: Modifications,
    /// Which records and alignment tags are kept.
    #[serde(default)]
    pub alignments: Alignments,
}

/// How per-base quality scores are stored.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "mode", rename_all = "kebab-case", deny_unknown_fields)]
pub enum Quality {
    /// Unchanged.
    Keep,
    /// Each score replaced by its bin's value.
    Bins {
        /// The bins, covering Phred 0–93.
        bins: Bins,
        /// Keep native scores in homopolymers at least this long.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        preserve_homopolymers: Option<usize>,
        /// Where those homopolymers are found.
        #[serde(default, skip_serializing_if = "is_default")]
        homopolymers_in: HomopolymerSource,
    },
    /// One score for every base of a read.
    PerRead {
        /// Where the read's score comes from.
        source: ReadQualitySource,
        /// The highest score written.
        max: u8,
    },
    /// Removed (`QUAL` = `*`). Some callers require qualities; prefer `per-read`.
    Drop,
}

/// Changes the command line makes to a built-in or custom policy (`--methylation`,
/// `--keep-ml`), the same for packing and for verifying.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Overrides {
    /// Which base-modification calls to keep.
    pub methylation: Option<ModificationCalls>,
    /// Keep `ML` probabilities exactly.
    pub keep_ml: bool,
}

/// Where protected homopolymers are found.
///
/// The default, for a policy that does not say, is [`Self::Reference`]: archives made by 0.2.1
/// and earlier record ONT policies without `homopolymers_in`, and did not protect unmapped reads.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum HomopolymerSource {
    /// Runs in the reference: bases aligned to them, and insertions anchored in them. Needs the
    /// alignment and the reference; unmapped reads and soft clips get no protection.
    #[default]
    Reference,
    /// Runs in the read's own bases, for reads without an alignment (such as FASTQ).
    Read,
    /// The reference's runs for mapped reads (soft clips get no protection), and the read's own
    /// for unmapped reads, such as Dorado's unaligned BAM. The built-in ONT policies use it.
    ReferenceOrRead,
}

fn is_default<T: Default + PartialEq>(value: &T) -> bool {
    *value == T::default()
}

/// The source of a read's single quality score.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ReadQualitySource {
    /// PacBio's predicted read accuracy tag `rq:f`; reads without it use `mean-error`.
    Rq,
    /// The Phred score of the read's mean per-base error probability.
    MeanError,
    /// The mean of the read's Phred scores.
    MeanPhred,
}

/// How base-modification calls (`MM`/`ML`) are stored.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Modifications {
    /// Quantize each `ML` probability to the midpoint of one of this many equal bins (2–256).
    /// Unset keeps them. Low-probability calls are never dropped: that rewrites `MM` deltas,
    /// which costs more than it saves once `ML` is quantized.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ml_levels: Option<u16>,
    /// Which calls are kept. Everything but `all` loses information some analyses need.
    #[serde(default, skip_serializing_if = "ModificationCalls::is_all")]
    pub calls: ModificationCalls,
}

/// Which base-modification calls are kept.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ModificationCalls {
    /// Every call, as basecalled.
    #[default]
    All,
    /// C modifications (5mC, 5hmC) only at CpG sites; other modifications (e.g. 6mA) kept.
    Cpg,
    /// C modifications only at CpG sites, and no other modifications.
    CpgOnly,
    /// No base-modification tags (`MM`, `ML`, `MN` removed).
    None,
}

impl ModificationCalls {
    /// Every mode, in display order.
    pub const ALL: [Self; 4] = [Self::All, Self::Cpg, Self::CpgOnly, Self::None];

    /// The mode's name on the command line and in TOML.
    pub fn name(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::Cpg => "cpg",
            Self::CpgOnly => "cpg-only",
            Self::None => "none",
        }
    }

    fn is_all(&self) -> bool {
        *self == Self::All
    }
}

impl fmt::Display for ModificationCalls {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

impl FromStr for ModificationCalls {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::ALL
            .into_iter()
            .find(|m| m.name() == s)
            .ok_or_else(|| {
                format!("unknown methylation mode {s:?} (expected all, cpg, cpg-only, or none)")
            })
    }
}

/// What happens to secondary alignments.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Secondary {
    /// Written out.
    #[default]
    Keep,
    /// Left out. Variant, SV, and modification callers ignore them, but read counts and
    /// `flagstat` change.
    Drop,
}

/// Which records and alignment tags are kept.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Alignments {
    /// Secondary alignments.
    #[serde(default)]
    pub secondary: Secondary,
    /// Remove `MD`/`NM` from records without `SEQ`, where CRAM cannot regenerate them and
    /// stores them verbatim. Records with `SEQ` always get them regenerated on decoding.
    #[serde(default = "default_true")]
    pub strip_md_without_sequence: bool,
}

fn default_true() -> bool {
    true
}

impl Default for Alignments {
    fn default() -> Self {
        Self {
            secondary: Secondary::Keep,
            strip_md_without_sequence: true,
        }
    }
}

/// Starts the archive header's `@CO` line that holds the policy as one line of TOML.
pub const HEADER_LINE_PREFIX: &str = "lrsqueeze policy toml: ";

/// The same line as written before the tool was renamed lrsqueeze (0.1.0 to 0.2.1).
const LEGACY_HEADER_LINE_PREFIXES: [&str; 1] = ["lrslim policy toml: "];

impl Policy {
    /// The `@CO` header line recording this policy: [`HEADER_LINE_PREFIX`], then the TOML with
    /// backslashes, newlines and tabs escaped.
    pub fn to_header_line(&self) -> String {
        let text = toml::to_string(self).expect("policies serialize to TOML");
        let mut line = String::from(HEADER_LINE_PREFIX);
        for c in text.trim_end().chars() {
            match c {
                '\\' => line.push_str("\\\\"),
                '\n' => line.push_str("\\n"),
                '\t' => line.push_str("\\t"),
                c => line.push(c),
            }
        }
        line
    }

    /// The policy in a header line written by [`to_header_line`](Self::to_header_line) (or by
    /// lrsqueeze under its earlier name, lrslim), or `None` when `line` is another kind of comment.
    ///
    /// # Errors
    ///
    /// Fails when the recorded TOML does not parse as a policy.
    pub fn from_header_line(line: &str) -> Option<Result<Self, toml::de::Error>> {
        let escaped = std::iter::once(HEADER_LINE_PREFIX)
            .chain(LEGACY_HEADER_LINE_PREFIXES)
            .find_map(|prefix| line.strip_prefix(prefix))?;
        let mut text = String::with_capacity(escaped.len());
        let mut chars = escaped.chars();
        while let Some(c) = chars.next() {
            if c != '\\' {
                text.push(c);
                continue;
            }
            match chars.next() {
                Some('n') => text.push('\n'),
                Some('t') => text.push('\t'),
                Some(other) => text.push(other),
                None => text.push('\\'),
            }
        }
        Some(toml::from_str(&text))
    }

    /// This policy with the command line's overrides applied.
    #[must_use]
    pub fn with_overrides(mut self, overrides: Overrides) -> Self {
        if let Some(calls) = overrides.methylation {
            self.modifications.calls = calls;
        }
        if overrides.keep_ml {
            self.modifications.ml_levels = None;
        }
        self
    }

    /// This policy as it applies to reads without an alignment (FASTQ): protected homopolymers
    /// are found in each read's own bases, and one score per read is the mean of its Phred
    /// scores, since a FASTQ has no `rq`. On PacBio Revio HG002 that matched the original's
    /// calls (indel F1 −0.005 pp), where the mean error probability lost 0.83 pp.
    #[must_use]
    pub fn for_fastq(mut self) -> Self {
        match &mut self.quality {
            Quality::Bins {
                preserve_homopolymers: Some(_),
                homopolymers_in,
                ..
            } => *homopolymers_in = HomopolymerSource::Read,
            Quality::PerRead {
                source: source @ ReadQualitySource::Rq,
                ..
            } => *source = ReadQualitySource::MeanPhred,
            _ => {}
        }
        self
    }

    /// The shortest reference homopolymer whose bases keep their native qualities, if reference
    /// homopolymers are protected.
    pub fn reference_homopolymer_length(&self) -> Option<usize> {
        match &self.quality {
            Quality::Bins {
                preserve_homopolymers,
                homopolymers_in: HomopolymerSource::Reference | HomopolymerSource::ReferenceOrRead,
                ..
            } => *preserve_homopolymers,
            _ => None,
        }
    }

    /// Whether unmapped reads keep native qualities in homopolymers (found in their own bases).
    pub fn protects_unmapped_homopolymers(&self) -> bool {
        matches!(
            self.quality,
            Quality::Bins {
                preserve_homopolymers: Some(_),
                homopolymers_in: HomopolymerSource::Read | HomopolymerSource::ReferenceOrRead,
                ..
            }
        )
    }

    /// Checks values that the types alone cannot.
    ///
    /// # Errors
    ///
    /// Describes the first invalid setting.
    pub fn validate(&self) -> Result<(), String> {
        if let Some(levels) = self.modifications.ml_levels
            && !(2..=256).contains(&levels)
        {
            return Err(format!("ml_levels must be between 2 and 256, not {levels}"));
        }
        if let Quality::Bins {
            preserve_homopolymers: Some(n),
            ..
        } = self.quality
            && n < 2
        {
            return Err(format!("preserve_homopolymers must be at least 2, not {n}"));
        }
        if let Quality::Bins {
            preserve_homopolymers: None,
            homopolymers_in: HomopolymerSource::Read | HomopolymerSource::ReferenceOrRead,
            ..
        } = self.quality
        {
            return Err("homopolymers_in needs preserve_homopolymers".to_owned());
        }
        if let Quality::PerRead { max, .. } = self.quality
            && max > crate::quality::MAX_PHRED
        {
            return Err(format!(
                "per-read max must be at most {}, not {max}",
                crate::quality::MAX_PHRED
            ));
        }
        Ok(())
    }

    /// A one-line summary for the output header's `@CO` record.
    pub fn summary(&self) -> String {
        let quality = match &self.quality {
            Quality::Keep => "keep".to_owned(),
            Quality::Bins {
                bins,
                preserve_homopolymers,
                homopolymers_in,
            } => match preserve_homopolymers {
                Some(n) => {
                    let source = match homopolymers_in {
                        HomopolymerSource::Reference => "reference",
                        HomopolymerSource::Read => "read",
                        HomopolymerSource::ReferenceOrRead => "reference (unmapped reads: read)",
                    };
                    format!("bins {bins}, native in {source} homopolymers >={n}bp")
                }
                None => format!("bins {bins}"),
            },
            Quality::PerRead { source, max } => format!("per-read from {source:?} (max {max})"),
            Quality::Drop => "dropped".to_owned(),
        };
        let ml = match self.modifications.ml_levels {
            Some(n) => format!("{n} levels"),
            None => "keep".to_owned(),
        };
        let ml = match self.modifications.calls {
            ModificationCalls::All => ml,
            calls => format!("{ml}, calls: {calls}"),
        };
        format!(
            "quality: {quality}; ML: {ml}; secondary: {:?}; strip MD/NM without SEQ: {}",
            self.alignments.secondary, self.alignments.strip_md_without_sequence
        )
    }

    /// The bins used, if any (for reporting).
    pub fn bins(&self) -> Option<&[Bin]> {
        match &self.quality {
            Quality::Bins { bins, .. } => Some(bins.as_slice()),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn presets_validate_and_round_trip_through_toml() {
        let presets = Platform::ALL
            .into_iter()
            .flat_map(|p| [Mode::Safe, Mode::Max].map(|m| Policy::preset(p, m)));
        for policy in presets.chain([Policy::lossless()]) {
            policy.validate().expect("built-in policies are valid");
            let text = toml::to_string(&policy).expect("serializes");
            let back: Policy = toml::from_str(&text).expect("parses its own output");
            assert_eq!(back, policy, "{text}");
        }
    }

    #[test]
    fn safe_presets_keep_secondaries_and_protect_ont_homopolymers() {
        for platform in Platform::ALL {
            let policy = Policy::preset(platform, Mode::Safe);
            assert_eq!(policy.alignments.secondary, Secondary::Keep, "{platform}");
            let Quality::Bins {
                preserve_homopolymers,
                ..
            } = policy.quality
            else {
                panic!("{platform}: safe mode bins qualities");
            };
            // Genome-wide, ONT indel accuracy depends on it; PacBio's was unchanged without it.
            let want = match platform {
                Platform::Ont => Some(5),
                Platform::Pacbio => None,
            };
            assert_eq!(preserve_homopolymers, want, "{platform}");
        }
    }

    #[test]
    fn policies_recorded_under_the_earlier_name_still_parse() {
        let policy = Policy::preset(Platform::Ont, Mode::Safe);
        let line = policy.to_header_line();
        let legacy = line.replacen(HEADER_LINE_PREFIX, "lrslim policy toml: ", 1);
        assert_eq!(
            Policy::from_header_line(&legacy).map(|p| p.expect("parses")),
            Some(policy)
        );
        assert!(Policy::from_header_line("lrsqueeze lossy policy: bins").is_none());
    }

    #[test]
    fn platform_names_parse_case_insensitively() {
        assert_eq!("ONT".parse::<Platform>(), Ok(Platform::Ont));
        assert_eq!("pacbio".parse::<Platform>(), Ok(Platform::Pacbio));
        assert!("illumina".parse::<Platform>().is_err());
    }

    #[test]
    fn policies_without_a_homopolymer_source_are_reference_only() {
        // An ONT policy as 0.2.1 and earlier recorded it, without `homopolymers_in`.
        let mut old = Policy::preset(Platform::Ont, Mode::Safe);
        if let Quality::Bins {
            homopolymers_in, ..
        } = &mut old.quality
        {
            *homopolymers_in = HomopolymerSource::Reference;
        }
        let text = toml::to_string(&old).expect("serializes");
        assert!(!text.contains("homopolymers_in"), "{text}");
        let back: Policy = toml::from_str(&text).expect("parses");
        assert_eq!(back, old);
        assert!(!back.protects_unmapped_homopolymers());

        // The built-in policy records its source, so verify tells the two apart.
        let new = Policy::preset(Platform::Ont, Mode::Safe);
        let text = toml::to_string(&new).expect("serializes");
        assert!(
            text.contains(r#"homopolymers_in = "reference-or-read""#),
            "{text}"
        );
        assert!(new.protects_unmapped_homopolymers());
        assert_ne!(new, old);
    }

    #[test]
    fn edited_toml_changes_bins() {
        let text = r#"
            [quality]
            mode = "bins"
            bins = [{ min = 0, max = 19, value = 12 }, { min = 20, max = 93, value = 35 }]
            preserve_homopolymers = 6

            [modifications]
            ml_levels = 32
        "#;
        let policy: Policy = toml::from_str(text).expect("valid");
        assert_eq!(policy.bins().map(<[Bin]>::len), Some(2));
        assert_eq!(
            policy.alignments,
            Alignments::default(),
            "omitted sections take defaults"
        );
        assert!(toml::from_str::<Policy>(&text.replace("min = 20", "min = 21")).is_err());
    }

    #[test]
    fn policies_round_trip_through_a_header_line() {
        let presets = Platform::ALL
            .into_iter()
            .flat_map(|p| [Mode::Safe, Mode::Max].map(|m| Policy::preset(p, m)));
        for policy in presets.chain([Policy::lossless()]) {
            let line = policy.to_header_line();
            assert!(!line.contains(['\n', '\t']), "{line}");
            let back = Policy::from_header_line(&line)
                .expect("a policy line")
                .expect("parses");
            assert_eq!(back, policy, "{line}");
        }
        assert!(Policy::from_header_line("lrsqueeze lossy policy: quality: keep").is_none());
    }

    #[test]
    fn rejects_out_of_range_settings() {
        let mut policy = Policy::preset(Platform::Ont, Mode::Max);
        policy.modifications.ml_levels = Some(1);
        assert!(policy.validate().is_err());

        let mut policy = Policy::preset(Platform::Pacbio, Mode::Safe);
        if let Quality::Bins {
            homopolymers_in, ..
        } = &mut policy.quality
        {
            *homopolymers_in = HomopolymerSource::Read;
        }
        assert!(
            policy.validate().is_err(),
            "read homopolymers need a length"
        );
    }

    #[test]
    fn fastq_policies_use_the_read() {
        for platform in [Platform::Ont, Platform::Pacbio] {
            for mode in [Mode::Safe, Mode::Max] {
                let policy = Policy::preset(platform, mode).for_fastq();
                policy.validate().expect("every preset applies to FASTQ");
                match policy.quality {
                    Quality::Bins {
                        preserve_homopolymers: Some(_),
                        homopolymers_in,
                        ..
                    } => assert_eq!(homopolymers_in, HomopolymerSource::Read),
                    Quality::Bins {
                        preserve_homopolymers: None,
                        homopolymers_in,
                        ..
                    } => assert_eq!(homopolymers_in, HomopolymerSource::Reference),
                    Quality::PerRead { source, .. } => {
                        assert_eq!(source, ReadQualitySource::MeanPhred);
                    }
                    _ => panic!("presets bin qualities or give one per read"),
                }
            }
        }
    }
}
