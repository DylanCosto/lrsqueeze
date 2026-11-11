//! Reading and editing encoded BAM records (without `block_size`) in place.
//!
//! Quality scores and `ML` values keep their length when transformed, so they are rewritten in
//! the record's own buffer; only removing tags rebuilds the auxiliary-data section.

use std::{io, ops::Range};

/// Offsets of the fixed-length fields (SAMv1 §4.2).
const NAME_LENGTH: usize = 8;
const CIGAR_OP_COUNT: Range<usize> = 12..14;
const FLAGS: Range<usize> = 14..16;
const SEQUENCE_LENGTH: Range<usize> = 16..20;
const POSITION: Range<usize> = 4..8;
const NAME: usize = 32;

/// The flag of secondary alignments.
const SECONDARY: u16 = 0x100;
/// The flag of unmapped reads.
const UNMAPPED: u16 = 0x4;

/// SAM letters for BAM's 4-bit base codes.
const BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

/// The flag of reverse-strand alignments.
const REVERSE: u16 = 0x10;

/// A quality score byte meaning the record has no qualities (when it fills `QUAL`).
pub const MISSING_QUALITY: u8 = 0xff;

/// The section boundaries of an encoded record.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Layout {
    cigar_start: usize,
    cigar_ops: usize,
    sequence_length: usize,
    quality_start: usize,
    data_start: usize,
}

fn u16_at(buf: &[u8], range: Range<usize>) -> u16 {
    u16::from_le_bytes(buf[range].try_into().expect("two bytes"))
}

fn u32_at(buf: &[u8], range: Range<usize>) -> u32 {
    u32::from_le_bytes(buf[range].try_into().expect("four bytes"))
}

fn invalid(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

impl Layout {
    /// Finds the sections of `buf`, checking that they fit.
    ///
    /// # Errors
    ///
    /// Fails when the declared lengths run past the end of the buffer.
    pub fn new(buf: &[u8]) -> io::Result<Self> {
        if buf.len() < NAME {
            return Err(invalid("BAM record shorter than its fixed fields"));
        }
        let cigar_start = NAME + usize::from(buf[NAME_LENGTH]);
        let cigar_ops = usize::from(u16_at(buf, CIGAR_OP_COUNT));
        let sequence_length = u32_at(buf, SEQUENCE_LENGTH) as usize;
        let sequence_start = cigar_start + 4 * cigar_ops;
        let quality_start = sequence_start + sequence_length.div_ceil(2);
        let data_start = quality_start + sequence_length;
        if data_start > buf.len() {
            return Err(invalid("BAM record sections run past its end"));
        }
        Ok(Self {
            cigar_start,
            cigar_ops,
            sequence_length,
            quality_start,
            data_start,
        })
    }

    /// The read length (`l_seq`); 0 when `SEQ` is `*`.
    pub fn sequence_length(&self) -> usize {
        self.sequence_length
    }

    /// The range of the 4-bit packed `SEQ` in the buffer.
    pub fn sequence(&self) -> Range<usize> {
        let start = self.cigar_start + 4 * self.cigar_ops;
        start..start + self.sequence_length.div_ceil(2)
    }

    /// Writes the bases of `SEQ` into `out` (cleared first), as SAM letters.
    pub fn bases_into(&self, buf: &[u8], out: &mut Vec<u8>) {
        let packed = &buf[self.sequence()];
        out.clear();
        out.extend((0..self.sequence_length).map(|i| {
            let b = packed[i / 2];
            BASES[usize::from(if i % 2 == 0 { b >> 4 } else { b & 0xf })]
        }));
    }

    /// The bases of `SEQ`, as SAM letters.
    pub fn bases(&self, buf: &[u8]) -> Vec<u8> {
        let mut bases = Vec::with_capacity(self.sequence_length);
        self.bases_into(buf, &mut bases);
        bases
    }

    /// The range of `QUAL` in the buffer.
    pub fn quality(&self) -> Range<usize> {
        self.quality_start..self.data_start
    }

    /// The range of the auxiliary data in a buffer of length `len`.
    pub fn data(&self, len: usize) -> Range<usize> {
        self.data_start..len
    }

    /// The CIGAR operations as `(op, length)`, with `op` the BAM code (0 = `M` ... 8 = `X`).
    ///
    /// `fields` are the record's auxiliary fields, where a CIGAR too long for the record's own
    /// field is kept (see [`Self::long_cigar`]).
    pub fn cigar<'a>(
        &self,
        buf: &'a [u8],
        fields: &[Field],
    ) -> impl Iterator<Item = (u8, u32)> + 'a {
        let range = self
            .long_cigar(buf, fields)
            .unwrap_or(self.cigar_start..self.cigar_start + 4 * self.cigar_ops);
        let (words, _) = buf[range].as_chunks::<4>();
        words.iter().map(|&w| {
            let v = u32::from_le_bytes(w);
            ((v & 0xf) as u8, v >> 4)
        })
    }

    /// The range of the real CIGAR when it is in a `CG` tag.
    ///
    /// A record with more than 65535 operations has `<l_seq>S<reference length>N` as its CIGAR
    /// and the real one in `CG:B:I` (SAMv1 §4.2.2).
    fn long_cigar(&self, buf: &[u8], fields: &[Field]) -> Option<Range<usize>> {
        const SOFT_CLIP: u32 = 4;
        const SKIP: u32 = 3;
        if self.cigar_ops != 2 || self.sequence_length == 0 {
            return None;
        }
        let word = |i: usize| {
            u32_at(
                buf,
                self.cigar_start + 4 * i..self.cigar_start + 4 * (i + 1),
            )
        };
        if word(0) != ((self.sequence_length as u32) << 4 | SOFT_CLIP) || word(1) & 0xf != SKIP {
            return None;
        }
        fields
            .iter()
            .find(|f| f.tag == *b"CG" && matches!(f.subtype, Some(b'I' | b'i')))
            .map(|f| f.value.clone())
    }
}

/// The record's flags.
pub fn flags(buf: &[u8]) -> u16 {
    u16_at(buf, FLAGS)
}

/// Whether the record is a secondary alignment.
pub fn is_secondary(buf: &[u8]) -> bool {
    flags(buf) & SECONDARY != 0
}

/// Whether the record is aligned to the reverse strand.
pub fn is_reverse(buf: &[u8]) -> bool {
    flags(buf) & REVERSE != 0
}

/// Whether the record is unmapped.
pub fn is_unmapped(buf: &[u8]) -> bool {
    flags(buf) & UNMAPPED != 0
}

/// The reference sequence ID, if any (and `None` for a record too short to have one).
pub fn reference_sequence_id(buf: &[u8]) -> Option<usize> {
    usize::try_from(i32::from_le_bytes(buf.get(0..4)?.try_into().ok()?)).ok()
}

/// The 0-based alignment start, if any (and `None` for a record too short to have one).
pub fn position(buf: &[u8]) -> Option<usize> {
    usize::try_from(i32::from_le_bytes(buf.get(POSITION)?.try_into().ok()?)).ok()
}

/// The number of read bases a CIGAR (`(op, length)` BAM codes) covers: `M`, `I`, `S`, `=`, `X`.
pub fn query_length(cigar: impl Iterator<Item = (u8, u32)>) -> u64 {
    cigar
        .filter(|&(op, _)| matches!(op, 0 | 1 | 4 | 7 | 8))
        .map(|(_, n)| u64::from(n))
        .sum()
}

/// One auxiliary field: its tag, type code, and where its value and whole field lie.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Field {
    /// The two-character tag.
    pub tag: [u8; 2],
    /// The SAM/BAM type code (`A`, `c`, `C`, `s`, `S`, `i`, `I`, `f`, `Z`, `H`, `B`).
    pub kind: u8,
    /// The value bytes. For `B` arrays, the elements only (after subtype and count).
    pub value: Range<usize>,
    /// For `B` arrays, the element subtype.
    pub subtype: Option<u8>,
    /// The whole field, tag included.
    pub span: Range<usize>,
}

fn scalar_size(kind: u8) -> Option<usize> {
    match kind {
        b'A' | b'c' | b'C' => Some(1),
        b's' | b'S' => Some(2),
        b'i' | b'I' | b'f' => Some(4),
        _ => None,
    }
}

/// The auxiliary fields of `data` (offsets relative to `base`, the data section's start).
///
/// # Errors
///
/// Fails on an unknown type code or a field that runs past the end.
pub fn fields(data: &[u8], base: usize) -> io::Result<Vec<Field>> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < data.len() {
        let start = i;
        let header = data
            .get(i..i + 3)
            .ok_or_else(|| invalid("truncated auxiliary field"))?;
        let (tag, kind) = ([header[0], header[1]], header[2]);
        i += 3;
        let (value, subtype) = match kind {
            b'Z' | b'H' => {
                let len = data[i..]
                    .iter()
                    .position(|&b| b == 0)
                    .ok_or_else(|| invalid("unterminated string field"))?;
                let value = i..i + len;
                i += len + 1;
                (value, None)
            }
            b'B' => {
                let sub = *data
                    .get(i)
                    .ok_or_else(|| invalid("truncated array field"))?;
                let size = scalar_size(sub).ok_or_else(|| invalid("invalid array subtype"))?;
                let count = data
                    .get(i + 1..i + 5)
                    .map(|c| u32::from_le_bytes(c.try_into().expect("four bytes")) as usize)
                    .ok_or_else(|| invalid("truncated array field"))?;
                i += 5;
                let value = i..i + size * count;
                i = value.end;
                (value, Some(sub))
            }
            _ => {
                let size = scalar_size(kind).ok_or_else(|| invalid("invalid field type"))?;
                let value = i..i + size;
                i += size;
                (value, None)
            }
        };
        if i > data.len() {
            return Err(invalid("auxiliary field runs past the record"));
        }
        out.push(Field {
            tag,
            kind,
            value: base + value.start..base + value.end,
            subtype,
            span: base + start..base + i,
        });
    }
    Ok(out)
}

/// Removes the fields whose spans are given (sorted, from [`fields`]) from `buf`.
pub fn remove_fields(buf: &mut Vec<u8>, spans: &[Range<usize>]) {
    for span in spans.iter().rev() {
        buf.drain(span.clone());
    }
}

/// Appends a `Z` (string) field.
pub fn push_string(buf: &mut Vec<u8>, tag: [u8; 2], value: &str) {
    buf.extend_from_slice(&tag);
    buf.push(b'Z');
    buf.extend_from_slice(value.as_bytes());
    buf.push(0);
}

/// Appends a `B:C` (unsigned byte array) field.
pub fn push_u8_array(buf: &mut Vec<u8>, tag: [u8; 2], values: &[u8]) {
    buf.extend_from_slice(&tag);
    buf.extend_from_slice(b"BC");
    buf.extend_from_slice(&(values.len() as u32).to_le_bytes());
    buf.extend_from_slice(values);
}

/// An integer field's value, whatever its width.
pub fn int_value(buf: &[u8], field: &Field) -> Option<i64> {
    let v = &buf[field.value.clone()];
    Some(match field.kind {
        b'c' => i64::from(v[0] as i8),
        b'C' => i64::from(v[0]),
        b's' => i64::from(i16::from_le_bytes(v.try_into().ok()?)),
        b'S' => i64::from(u16::from_le_bytes(v.try_into().ok()?)),
        b'i' => i64::from(i32::from_le_bytes(v.try_into().ok()?)),
        b'I' => i64::from(u32::from_le_bytes(v.try_into().ok()?)),
        _ => return None,
    })
}

/// A float field's value.
pub fn f32_value(buf: &[u8], field: &Field) -> Option<f32> {
    (field.kind == b'f')
        .then(|| f32::from_le_bytes(buf[field.value.clone()].try_into().expect("four bytes")))
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// Encodes a minimal record: name `r`, the given CIGAR, sequence of `len` `A`s, `quals`,
    /// and raw auxiliary bytes.
    pub(crate) fn record(
        flags: u16,
        pos: i32,
        cigar: &[(u8, u32)],
        quals: &[u8],
        aux: &[u8],
    ) -> Vec<u8> {
        let len = quals.len();
        let mut b = Vec::new();
        b.extend(0i32.to_le_bytes());
        b.extend(pos.to_le_bytes());
        b.push(2); // l_read_name, "r\0"
        b.push(60);
        b.extend(4680u16.to_le_bytes());
        b.extend((cigar.len() as u16).to_le_bytes());
        b.extend(flags.to_le_bytes());
        b.extend((len as u32).to_le_bytes());
        b.extend((-1i32).to_le_bytes());
        b.extend((-1i32).to_le_bytes());
        b.extend(0i32.to_le_bytes());
        b.extend(b"r\0");
        for &(op, n) in cigar {
            b.extend(((n << 4) | u32::from(op)).to_le_bytes());
        }
        b.extend(std::iter::repeat_n(0x11, len.div_ceil(2)));
        b.extend(quals);
        b.extend(aux);
        b
    }

    #[test]
    fn layout_and_fields() {
        let aux = [
            b"NMi".as_slice(),
            &3i32.to_le_bytes(),
            b"MDZ10\0",
            b"MLBC",
            &2u32.to_le_bytes(),
            &[7, 200],
        ]
        .concat();
        let buf = record(0, 100, &[(0, 4)], &[30, 31, 32, 33], &aux);
        let layout = Layout::new(&buf).expect("valid");
        assert_eq!(&buf[layout.quality()], &[30, 31, 32, 33]);
        assert_eq!(position(&buf), Some(100));
        let data = layout.data(buf.len());
        let f = fields(&buf[data.clone()], data.start).expect("valid aux");
        assert_eq!(layout.cigar(&buf, &f).collect::<Vec<_>>(), vec![(0, 4)]);
        assert_eq!(
            f.iter().map(|f| f.tag).collect::<Vec<_>>(),
            vec![*b"NM", *b"MD", *b"ML"]
        );
        assert_eq!(&buf[f[1].value.clone()], b"10");
        assert_eq!(&buf[f[2].value.clone()], &[7, 200]);
        assert_eq!(f[2].subtype, Some(b'C'));

        let mut edited = buf.clone();
        remove_fields(&mut edited, &[f[0].span.clone(), f[1].span.clone()]);
        let data = layout.data(edited.len());
        let f = fields(&edited[data.clone()], data.start).expect("valid aux");
        assert_eq!(f.iter().map(|f| f.tag).collect::<Vec<_>>(), vec![*b"ML"]);
    }

    /// The `CG:B:I` field holding `cigar`.
    pub(crate) fn cg_field(cigar: &[(u8, u32)]) -> Vec<u8> {
        let mut aux = b"CGBI".to_vec();
        aux.extend((cigar.len() as u32).to_le_bytes());
        for &(op, n) in cigar {
            aux.extend(((n << 4) | u32::from(op)).to_le_bytes());
        }
        aux
    }

    fn cigar_of(buf: &[u8]) -> Vec<(u8, u32)> {
        let layout = Layout::new(buf).expect("valid");
        let data = layout.data(buf.len());
        let f = fields(&buf[data.clone()], data.start).expect("valid aux");
        layout.cigar(buf, &f).collect()
    }

    #[test]
    fn long_cigar_comes_from_the_cg_tag() {
        let real = [(4, 1), (0, 2), (1, 1), (2, 3), (0, 2)];
        let placeholder = [(4, 6), (3, 7)];
        let buf = record(0, 0, &placeholder, &[30; 6], &cg_field(&real));
        assert_eq!(cigar_of(&buf), real);

        // A fully clipped read with no CG tag is what it says it is.
        let buf = record(0, 0, &placeholder, &[30; 6], b"");
        assert_eq!(cigar_of(&buf), placeholder);

        // A CG tag on a record whose CIGAR is not the placeholder is not read.
        for cigar in [[(4, 5), (3, 7)], [(4, 6), (2, 7)]] {
            let buf = record(0, 0, &cigar, &[30; 6], &cg_field(&real));
            assert_eq!(cigar_of(&buf), cigar);
        }
    }

    #[test]
    fn decodes_bases_and_appends_fields() {
        let mut buf = record(0x10, 5, &[(0, 3)], &[30, 31, 32], b"");
        let layout = Layout::new(&buf).expect("valid");
        assert_eq!(
            layout.bases(&buf),
            b"AAA",
            "the test record packs A (code 1)"
        );
        assert!(is_reverse(&buf));
        push_string(&mut buf, *b"MM", "C+m?,0;");
        push_u8_array(&mut buf, *b"ML", &[200]);
        let data = layout.data(buf.len());
        let f = fields(&buf[data.clone()], data.start).expect("valid aux");
        assert_eq!(&buf[f[0].value.clone()], b"C+m?,0;");
        assert_eq!(
            (f[1].subtype, &buf[f[1].value.clone()]),
            (Some(b'C'), &[200u8][..])
        );
    }

    #[test]
    fn rejects_truncated_records() {
        let buf = record(0, 0, &[(0, 4)], &[30, 31, 32, 33], b"");
        assert!(Layout::new(&buf[..buf.len() - 1]).is_err());
        assert!(fields(b"MLBC\x05\0\0\0\x01", 0).is_err());
    }
}
