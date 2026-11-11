//! Homopolymer runs: in the reference, projected onto read positions, or in a read itself.
//!
//! The protected positions are those of the validated `hybrid4` policy: query bases aligned
//! (`M`/`=`/`X`) to a reference A/C/G/T run of at least `min_length` (case-insensitive, no
//! flanks), and inserted bases whose reference anchor on either side is in such a run. Soft
//! clips and unmapped reads are never protected.

/// Half-open reference intervals of homopolymer runs, sorted and disjoint.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Runs {
    starts: Vec<u32>,
    ends: Vec<u32>,
}

impl Runs {
    /// Finds the runs of one base (A, C, G, or T, either case) at least `min_length` long.
    pub fn find(sequence: &[u8], min_length: usize) -> Self {
        let mut runs = Self::default();
        let mut i = 0;
        while i < sequence.len() {
            let base = sequence[i].to_ascii_uppercase();
            let mut j = i + 1;
            while j < sequence.len() && sequence[j].to_ascii_uppercase() == base {
                j += 1;
            }
            if matches!(base, b'A' | b'C' | b'G' | b'T') && j - i >= min_length {
                runs.starts
                    .push(u32::try_from(i).expect("reference sequences fit in u32"));
                runs.ends
                    .push(u32::try_from(j).expect("reference sequences fit in u32"));
            }
            i = j;
        }
        runs
    }

    /// The number of runs.
    pub fn len(&self) -> usize {
        self.starts.len()
    }

    /// Whether there are no runs.
    pub fn is_empty(&self) -> bool {
        self.starts.is_empty()
    }

    /// Whether reference position `pos` lies in a run.
    pub fn contains(&self, pos: usize) -> bool {
        // The first run ending after `pos` is the only one that can contain it.
        let i = self.ends.partition_point(|&end| end as usize <= pos);
        i < self.starts.len() && self.starts[i] as usize <= pos
    }

    /// Marks in `keep` every position inside a run, for runs found in a read's own bases.
    pub fn mark(&self, keep: &mut [bool]) {
        for (&start, &end) in self.starts.iter().zip(&self.ends) {
            keep[start as usize..end as usize].fill(true);
        }
    }

    /// Marks in `keep` (one entry per query base) the bases of an alignment starting at
    /// `start` with `cigar` (`(op, len)` BAM codes) that the policy protects.
    pub fn project(&self, start: usize, cigar: impl Iterator<Item = (u8, u32)>, keep: &mut [bool]) {
        let (mut rp, mut qp) = (start, 0usize);
        for (op, n) in cigar {
            let n = n as usize;
            match op {
                // M, =, X: aligned bases.
                0 | 7 | 8 => {
                    let stop = rp + n;
                    let mut j = self.ends.partition_point(|&end| end as usize <= rp);
                    while j < self.starts.len() && (self.starts[j] as usize) < stop {
                        let lo = rp.max(self.starts[j] as usize);
                        let hi = stop.min(self.ends[j] as usize);
                        keep[qp + lo - rp..qp + hi - rp].fill(true);
                        j += 1;
                    }
                    rp = stop;
                    qp += n;
                }
                // I: protected when either reference anchor is in a run.
                1 => {
                    if (rp > 0 && self.contains(rp - 1)) || self.contains(rp) {
                        keep[qp..qp + n].fill(true);
                    }
                    qp += n;
                }
                // D, N: reference only.
                2 | 3 => rp += n,
                // S: query only, never protected.
                4 => qp += n,
                // H, P: neither.
                _ => {}
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_runs_case_insensitively_and_skips_n() {
        let runs = Runs::find(b"ACaaaAAGTTTTNNNNNNcccccc", 5);
        // aaaAA at 2..7, NNNNNN ignored, cccccc at 18..24; TTTT is too short.
        assert_eq!(
            (runs.starts.clone(), runs.ends.clone()),
            (vec![2, 18], vec![7, 24])
        );
        assert!(runs.contains(2) && runs.contains(6) && !runs.contains(7) && !runs.contains(1));
    }

    #[test]
    fn projects_matches_insertions_and_clips() {
        // Reference run at 10..15.
        let runs = Runs {
            starts: vec![10],
            ends: vec![15],
        };
        // 2S 4M 3I 4M 1D 2M starting at reference 8:
        // 4M covers ref 8..12 (query 2..6) -> ref 10,11 protected -> query 4,5.
        // 3I anchored between ref 11 and 12, both in the run -> query 6..9.
        // 4M covers ref 12..16 (query 9..13) -> ref 12..15 -> query 9,10,11.
        // 1D skips ref 16; 2M ref 17..19 not in a run.
        let cigar = [(4, 2), (0, 4), (1, 3), (0, 4), (2, 1), (0, 2)];
        let mut keep = vec![false; 15];
        runs.project(8, cigar.into_iter(), &mut keep);
        let protected: Vec<usize> = (0..15).filter(|&i| keep[i]).collect();
        assert_eq!(protected, vec![4, 5, 6, 7, 8, 9, 10, 11]);
    }

    #[test]
    fn insertion_at_run_edge_uses_either_anchor() {
        let runs = Runs {
            starts: vec![10],
            ends: vec![15],
        };
        // Insertion right after the run: left anchor 14 is in it.
        let mut keep = vec![false; 7];
        runs.project(10, [(0, 5), (1, 2)].into_iter(), &mut keep);
        assert!(keep[5] && keep[6]);
        // Insertion right before the run: right anchor 10 is in it.
        let mut keep = vec![false; 3];
        runs.project(8, [(0, 2), (1, 1)].into_iter(), &mut keep);
        assert!(keep[2]);
        // Insertion two bases before: neither anchor (8, 9) is in it.
        let mut keep = vec![false; 2];
        runs.project(8, [(0, 1), (1, 1)].into_iter(), &mut keep);
        assert!(!keep[1]);
    }
}
