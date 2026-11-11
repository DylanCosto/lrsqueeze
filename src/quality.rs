//! Quality-score bins: contiguous Phred ranges, each stored as one representative value.

use std::{fmt, str::FromStr};

use serde::{Deserialize, Serialize};

/// The highest Phred score SAM can represent (`~` is 93 + 33).
pub const MAX_PHRED: u8 = 93;

/// A range of input Phred scores, inclusive, stored as `value`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Bin {
    /// The lowest input score in the bin.
    pub min: u8,
    /// The highest input score in the bin.
    pub max: u8,
    /// The score written for every input in the bin.
    pub value: u8,
}

/// Bins covering every Phred score from 0 to [`MAX_PHRED`], in order, without gaps or overlaps.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(try_from = "Vec<Bin>", into = "Vec<Bin>")]
pub struct Bins(Vec<Bin>);

/// Why a set of bins was rejected.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BinsError {
    /// No bins were given.
    Empty,
    /// A bin starts somewhere other than just after the previous one (or at 0).
    Gap { expected_min: u8, found_min: u8 },
    /// A bin's `min` is above its `max`.
    Inverted(Bin),
    /// A `max` or `value` above [`MAX_PHRED`].
    OutOfRange(Bin),
    /// The last bin ends below [`MAX_PHRED`].
    Incomplete { last_max: u8 },
    /// A bin in a `MIN-MAX:VALUE` list could not be parsed.
    Syntax(String),
}

impl fmt::Display for BinsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "no quality bins given"),
            Self::Gap {
                expected_min,
                found_min,
            } => write!(
                f,
                "quality bins must be contiguous: expected a bin starting at {expected_min}, found {found_min}"
            ),
            Self::Inverted(b) => write!(f, "quality bin {}-{} has min above max", b.min, b.max),
            Self::OutOfRange(b) => write!(
                f,
                "quality bin {}-{}:{} exceeds Phred {MAX_PHRED}",
                b.min, b.max, b.value
            ),
            Self::Incomplete { last_max } => write!(
                f,
                "quality bins must reach Phred {MAX_PHRED}; the last ends at {last_max}"
            ),
            Self::Syntax(s) => write!(f, "invalid quality bin {s:?}: expected MIN-MAX:VALUE"),
        }
    }
}

impl std::error::Error for BinsError {}

impl Bins {
    /// Validates bins that must start at 0, be contiguous, and end at [`MAX_PHRED`].
    ///
    /// # Errors
    ///
    /// Returns a [`BinsError`] describing the first violation.
    pub fn new(bins: Vec<Bin>) -> Result<Self, BinsError> {
        let mut expected_min = 0u8;
        for &bin in &bins {
            if bin.min != expected_min {
                return Err(BinsError::Gap {
                    expected_min,
                    found_min: bin.min,
                });
            }
            if bin.min > bin.max {
                return Err(BinsError::Inverted(bin));
            }
            if bin.max > MAX_PHRED || bin.value > MAX_PHRED {
                return Err(BinsError::OutOfRange(bin));
            }
            expected_min = bin.max.saturating_add(1);
        }
        match bins.last() {
            None => Err(BinsError::Empty),
            Some(last) if last.max < MAX_PHRED => Err(BinsError::Incomplete { last_max: last.max }),
            Some(_) => Ok(Self(bins)),
        }
    }

    /// The bins, lowest first.
    pub fn as_slice(&self) -> &[Bin] {
        &self.0
    }

    /// The number of bins.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Always false: validated bins are never empty.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// A lookup table from input score to output score. Scores above [`MAX_PHRED`] (outside
    /// SAM's range) map to the last bin's value.
    pub fn lut(&self) -> [u8; 256] {
        let last = self.0.last().expect("validated bins are not empty").value;
        let mut lut = [last; 256];
        for bin in &self.0 {
            lut[usize::from(bin.min)..=usize::from(bin.max)].fill(bin.value);
        }
        lut
    }
}

impl TryFrom<Vec<Bin>> for Bins {
    type Error = BinsError;

    fn try_from(bins: Vec<Bin>) -> Result<Self, Self::Error> {
        Self::new(bins)
    }
}

impl From<Bins> for Vec<Bin> {
    fn from(bins: Bins) -> Self {
        bins.0
    }
}

/// Parses `MIN-MAX:VALUE,...`, e.g. `0-9:7,10-19:15,20-29:27,30-93:38`.
impl FromStr for Bins {
    type Err = BinsError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parse_bin = |item: &str| -> Option<Bin> {
            let (range, value) = item.trim().split_once(':')?;
            let (min, max) = range.split_once('-')?;
            Some(Bin {
                min: min.trim().parse().ok()?,
                max: max.trim().parse().ok()?,
                value: value.trim().parse().ok()?,
            })
        };
        let bins = s
            .split(',')
            .map(|item| parse_bin(item).ok_or_else(|| BinsError::Syntax(item.to_owned())))
            .collect::<Result<Vec<_>, _>>()?;
        Self::new(bins)
    }
}

impl fmt::Display for Bins {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, b) in self.0.iter().enumerate() {
            if i > 0 {
                f.write_str(",")?;
            }
            write!(f, "{}-{}:{}", b.min, b.max, b.value)?;
        }
        Ok(())
    }
}

/// The Phred score of a read accuracy such as PacBio's `rq`, capped at `max`.
///
/// Rounds half to even, matching the Python prototype the presets were validated with.
pub fn accuracy_to_phred(accuracy: f32, max: u8) -> u8 {
    /// Error rates below this (Q40) are clamped before taking the logarithm.
    const MIN_ERROR: f64 = 1e-4;
    let error = (1.0 - f64::from(accuracy)).max(MIN_ERROR);
    let q = (-10.0 * error.log10()).round_ties_even();
    // `q` is in 0..=40 by the clamp above, and `max` caps it further.
    q.clamp(0.0, f64::from(max)) as u8
}

/// The mean of a read's Phred scores, rounded and capped at `max`.
pub fn mean_phred(quals: &[u8], max: u8) -> u8 {
    if quals.is_empty() {
        return max;
    }
    let sum: u64 = quals.iter().map(|&q| u64::from(q)).sum();
    let q = (sum as f64 / quals.len() as f64).round_ties_even();
    q.min(f64::from(max)) as u8
}

/// The Phred score of a read's mean per-base error probability, capped at `max`.
pub fn mean_error_phred(quals: &[u8], max: u8) -> u8 {
    if quals.is_empty() {
        return max;
    }
    let mean: f64 = quals
        .iter()
        .map(|&q| 10f64.powf(-f64::from(q) / 10.0))
        .sum::<f64>()
        / quals.len() as f64;
    let q = (-10.0 * mean.log10()).round_ties_even();
    q.clamp(0.0, f64::from(max)) as u8
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_validates_and_round_trips() {
        let spec = "0-9:7,10-19:15,20-29:27,30-93:38";
        let bins: Bins = spec.parse().expect("valid spec");
        assert_eq!(bins.len(), 4);
        assert_eq!(bins.to_string(), spec);
        let lut = bins.lut();
        assert_eq!(
            (lut[0], lut[9], lut[10], lut[29], lut[30], lut[93]),
            (7, 7, 15, 27, 38, 38)
        );
        assert_eq!(lut[200], 38, "out-of-range scores fall in the last bin");
    }

    #[test]
    fn rejects_gaps_overlaps_and_short_coverage() {
        assert!(matches!(
            "0-9:7,11-93:20".parse::<Bins>(),
            Err(BinsError::Gap { .. })
        ));
        assert!(matches!(
            "0-9:7,9-93:20".parse::<Bins>(),
            Err(BinsError::Gap { .. })
        ));
        assert!(matches!(
            "0-50:7".parse::<Bins>(),
            Err(BinsError::Incomplete { .. })
        ));
        assert!(matches!(
            "1-93:7".parse::<Bins>(),
            Err(BinsError::Gap { .. })
        ));
        assert!(matches!(
            "0-93:94".parse::<Bins>(),
            Err(BinsError::OutOfRange(_))
        ));
        assert!(matches!("0-93".parse::<Bins>(), Err(BinsError::Syntax(_))));
    }

    #[test]
    fn accuracy_phred_matches_prototype() {
        // Python: min(40, round(-10*log10(max(1-rq, 1e-4)))) on f32 rq values.
        assert_eq!(accuracy_to_phred(0.999_854, 40), 38);
        assert_eq!(accuracy_to_phred(0.99, 40), 20);
        assert_eq!(accuracy_to_phred(1.0, 40), 40);
        assert_eq!(accuracy_to_phred(0.5, 40), 3);
    }

    #[test]
    fn mean_phred_averages_scores() {
        assert_eq!(mean_phred(&[20, 40], 60), 30);
        // One low score barely moves it, unlike the mean error probability.
        assert_eq!(mean_phred(&[40, 40, 40, 3], 60), 31);
        assert_eq!(mean_error_phred(&[40, 40, 40, 3], 60), 9);
        assert_eq!(mean_phred(&[93, 93], 40), 40, "capped");
        assert_eq!(mean_phred(&[], 40), 40, "no scores");
    }

    #[test]
    fn mean_error_phred_weights_errors() {
        assert_eq!(mean_error_phred(&[20, 20], 60), 20);
        // One Q10 base dominates a Q40 base: mean error (0.1 + 0.0001) / 2.
        assert_eq!(mean_error_phred(&[10, 40], 60), 13);
    }
}
