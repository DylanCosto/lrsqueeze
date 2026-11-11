//! Keeping base-modification calls (`MM`/`ML`) only at CpG sites.
//!
//! `MM` counts positions among the bases of one kind in the read as sequenced (`SEQ`, reverse
//! complemented back for reverse-strand records), and `ML` holds one probability per call and
//! modification. A C call is kept when the read has CpG there or the base is aligned to a
//! reference CpG (the C of CG on the forward strand; the G, for reverse-strand records).
//!
//! Kept parts become explicit (`?`). Calls the basecaller left out of an implicit (`.`) part
//! mean "unmodified", so at CpGs they become explicit 0-probability calls: CpG pileups keep
//! their coverage, and non-CpG positions read as "no call", not "unmodified".

use std::{fmt::Write, io};

fn invalid(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

/// One `MM` part with its `ML` values.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Part<'a> {
    base: u8,
    strand: u8,
    /// The modification codes, as written (letters, or one ChEBI number).
    codes: &'a str,
    /// The mode character, if any (`.` or `?`).
    mode: Option<u8>,
    deltas: Vec<usize>,
    /// `codes.len()` values per call (one for a ChEBI number).
    values: &'a [u8],
}

impl Part<'_> {
    fn n_codes(&self) -> usize {
        if self.codes.bytes().all(|b| b.is_ascii_digit()) {
            1
        } else {
            self.codes.len()
        }
    }
}

fn parse<'a>(mm: &'a str, ml: &'a [u8]) -> io::Result<Vec<Part<'a>>> {
    let mut parts = Vec::new();
    let mut k = 0;
    for part in mm.split(';').filter(|p| !p.is_empty()) {
        let mut items = part.split(',');
        let head = items
            .next()
            .ok_or_else(|| invalid("empty MM part"))?
            .as_bytes();
        if head.len() < 3 {
            return Err(invalid("MM part head too short"));
        }
        let (mode, codes_end) = match head[head.len() - 1] {
            m @ (b'.' | b'?') => (Some(m), head.len() - 1),
            _ => (None, head.len()),
        };
        let codes = std::str::from_utf8(&head[2..codes_end]).map_err(|_| invalid("MM codes"))?;
        let deltas = items
            .map(|d| d.parse::<usize>().map_err(|_| invalid("MM delta")))
            .collect::<io::Result<Vec<_>>>()?;
        let mut p = Part {
            base: head[0],
            strand: head[1],
            codes,
            mode,
            deltas,
            values: &[],
        };
        let n = p.n_codes() * p.deltas.len();
        p.values = ml
            .get(k..k + n)
            .ok_or_else(|| invalid("ML shorter than MM"))?;
        k += n;
        parts.push(p);
    }
    if k != ml.len() {
        return Err(invalid("ML longer than MM"));
    }
    Ok(parts)
}

/// Rewrites `mm`/`ml` keeping C-modification calls only where `is_cpg(position)` holds
/// (positions in the original read orientation, of length `seq.len()`). With `drop_other`, the
/// parts for other bases or strands are left out; otherwise they are kept as they are.
///
/// # Errors
///
/// Fails on malformed `MM`/`ML`.
pub fn restrict_to_cpg(
    mm: &str,
    ml: &[u8],
    seq: &[u8],
    is_cpg: impl Fn(usize) -> bool,
    drop_other: bool,
) -> io::Result<(String, Vec<u8>)> {
    let mut out_mm = String::new();
    let mut out_ml = Vec::new();
    let c_positions: Vec<usize> = seq
        .iter()
        .enumerate()
        .filter(|(_, b)| b.eq_ignore_ascii_case(&b'C'))
        .map(|(i, _)| i)
        .collect();
    for part in parse(mm, ml)? {
        let n = part.n_codes();
        if part.base != b'C' || part.strand != b'+' {
            if !drop_other {
                out_mm.push(part.base as char);
                out_mm.push(part.strand as char);
                out_mm.push_str(part.codes);
                if let Some(m) = part.mode {
                    out_mm.push(m as char);
                }
                for d in &part.deltas {
                    write!(out_mm, ",{d}").expect("writing to a String cannot fail");
                }
                out_mm.push(';');
                out_ml.extend_from_slice(part.values);
            }
            continue;
        }
        // Calls by index among the read's Cs.
        let mut called = vec![None; c_positions.len()];
        let mut idx = 0usize;
        for (i, d) in part.deltas.iter().enumerate() {
            idx += d;
            let slot = called
                .get_mut(idx)
                .ok_or_else(|| invalid("MM delta past the read's Cs"))?;
            *slot = Some(&part.values[i * n..(i + 1) * n]);
            idx += 1;
        }
        // In an implicit part, an uncalled C is unmodified: written as an explicit 0 call.
        let implicit = part.mode != Some(b'?');
        let unmodified = vec![0u8; n];
        out_mm.push_str("C+");
        out_mm.push_str(part.codes);
        out_mm.push('?');
        let mut gap = 0usize;
        for (j, &pos) in c_positions.iter().enumerate() {
            let values = match called[j] {
                _ if !is_cpg(pos) => None,
                Some(v) => Some(v),
                None if implicit => Some(&unmodified[..]),
                None => None,
            };
            match values {
                Some(v) => {
                    write!(out_mm, ",{gap}").expect("writing to a String cannot fail");
                    out_ml.extend_from_slice(v);
                    gap = 0;
                }
                None => gap += 1,
            }
        }
        out_mm.push(';');
    }
    Ok((out_mm, out_ml))
}

/// The read's bases in their original orientation (reverse complemented for reverse-strand
/// records), uppercase.
pub fn original_orientation(seq: &[u8], reverse: bool) -> Vec<u8> {
    if !reverse {
        return seq.to_ascii_uppercase();
    }
    seq.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

/// Marks, in original-orientation read coordinates, the bases aligned to a reference CpG.
/// `reference` is the whole reference sequence; `cigar` is `(op, len)` BAM codes.
pub fn reference_cpg(
    start: usize,
    cigar: impl Iterator<Item = (u8, u32)>,
    reference: &[u8],
    read_len: usize,
    reverse: bool,
) -> Vec<bool> {
    let mut marks = vec![false; read_len];
    let up = |i: usize| reference.get(i).map(u8::to_ascii_uppercase);
    let (mut rp, mut qp) = (start, 0usize);
    for (op, n) in cigar {
        let n = n as usize;
        match op {
            0 | 7 | 8 => {
                for i in 0..n {
                    let (r, q) = (rp + i, qp + i);
                    let hit = if reverse {
                        r > 0 && up(r) == Some(b'G') && up(r - 1) == Some(b'C')
                    } else {
                        up(r) == Some(b'C') && up(r + 1) == Some(b'G')
                    };
                    if hit && q < read_len {
                        marks[if reverse { read_len - 1 - q } else { q }] = true;
                    }
                }
                rp += n;
                qp += n;
            }
            1 | 4 => qp += n,
            2 | 3 => rp += n,
            _ => {}
        }
    }
    marks
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Decodes calls as (part head without mode, read position, values).
    fn calls(mm: &str, ml: &[u8], seq: &[u8]) -> Vec<(String, usize, Vec<u8>)> {
        let mut out = Vec::new();
        for p in parse(mm, ml).expect("valid") {
            let n = p.n_codes();
            let positions: Vec<usize> = seq
                .iter()
                .enumerate()
                .filter(|(_, b)| **b == p.base)
                .map(|(i, _)| i)
                .collect();
            let mut idx = 0;
            for (i, d) in p.deltas.iter().enumerate() {
                idx += d;
                let head = format!("{}{}{}", p.base as char, p.strand as char, p.codes);
                out.push((head, positions[idx], p.values[i * n..(i + 1) * n].to_vec()));
                idx += 1;
            }
        }
        out
    }

    #[test]
    fn keeps_cpg_calls_and_makes_implicit_ones_explicit() {
        // Cs at 0 (CG), 3 (CA), 5 (CG), 8 (CG). Implicit part calls Cs 0 and 3 only.
        let seq = b"CGACACGACG";
        let mm = "C+hm.,0,0;A+a.,1;";
        let ml = [10, 200, 50, 60, 77];
        let (new_mm, new_ml) =
            restrict_to_cpg(mm, &ml, seq, |p| seq.get(p + 1) == Some(&b'G'), false).expect("valid");
        assert_eq!(new_mm, "C+hm?,0,1,0;A+a.,1;");
        assert_eq!(
            calls(&new_mm, &new_ml, seq),
            vec![
                ("C+hm".into(), 0, vec![10, 200]),
                ("C+hm".into(), 5, vec![0, 0]),
                ("C+hm".into(), 8, vec![0, 0]),
                ("A+a".into(), 4, vec![77]),
            ]
        );
    }

    #[test]
    fn explicit_parts_omit_uncalled_cpgs_and_drop_other_removes_6ma() {
        let seq = b"CGACACGACG";
        let mm = "C+m?,0,2;A+a.,0;";
        let ml = [240, 30, 99];
        let (new_mm, new_ml) =
            restrict_to_cpg(mm, &ml, seq, |p| seq.get(p + 1) == Some(&b'G'), true).expect("valid");
        // Calls at Cs 0 and 8 (deltas 0, 2); the CpG C at 5 was not called in an explicit
        // part, so it stays absent, and 6mA is dropped.
        assert_eq!(new_mm, "C+m?,0,2;");
        assert_eq!(
            calls(&new_mm, &new_ml, seq),
            vec![("C+m".into(), 0, vec![240]), ("C+m".into(), 8, vec![30])]
        );
    }

    #[test]
    fn rejects_inconsistent_tags() {
        assert!(restrict_to_cpg("C+m.,0,0;", &[1], b"CC", |_| true, false).is_err());
        assert!(restrict_to_cpg("C+m.,5;", &[1], b"CC", |_| true, false).is_err());
    }

    #[test]
    fn reference_cpg_marks_both_strands() {
        //             0123456789
        let reference = b"AACGTTACGA";
        // 6M at reference 1: read positions 0..6 align to reference 1..7.
        let fwd = reference_cpg(1, [(0u8, 6u32)].into_iter(), reference, 6, false);
        // Reference C at 2 (CG) -> read 1.
        assert_eq!(fwd, vec![false, true, false, false, false, false]);
        // Reverse: the reference G of CG at 3 -> SEQ position 2 -> original position 6-1-2 = 3.
        let rev = reference_cpg(1, [(0u8, 6u32)].into_iter(), reference, 6, true);
        assert_eq!(rev, vec![false, false, false, true, false, false]);
    }

    #[test]
    fn original_orientation_reverse_complements() {
        assert_eq!(original_orientation(b"aCGTn", true), b"NACGT");
        assert_eq!(original_orientation(b"aCGTn", false), b"ACGTN");
    }
}
