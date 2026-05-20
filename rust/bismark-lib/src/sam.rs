use std::io::{self, BufRead};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SamError {
    #[error("malformed SAM record: {0}")]
    Malformed(String),
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
}

/// Read conversion strand (XR tag): CT or GA.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadConversion {
    CT,
    GA,
}

/// Genome conversion used for bisulfite index (XG tag): CT or GA.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenomeConversion {
    CT,
    GA,
}

/// Bismark strand orientation derived from XR + XG combination.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BismarkStrand {
    OT,   // Original Top    — XR=CT, XG=CT
    CTOT, // Comp. to OT     — XR=GA, XG=CT
    CTOB, // Comp. to OB     — XR=GA, XG=GA
    OB,   // Original Bottom — XR=CT, XG=GA
}

impl BismarkStrand {
    pub fn from_tags(xr: ReadConversion, xg: GenomeConversion) -> Self {
        match (xr, xg) {
            (ReadConversion::CT, GenomeConversion::CT) => BismarkStrand::OT,
            (ReadConversion::GA, GenomeConversion::CT) => BismarkStrand::CTOT,
            (ReadConversion::GA, GenomeConversion::GA) => BismarkStrand::CTOB,
            (ReadConversion::CT, GenomeConversion::GA) => BismarkStrand::OB,
        }
    }
}

/// Methylation counts from an XM string.
#[derive(Debug, Default, Clone)]
pub struct MethylationCounts {
    pub cpg_meth: u32,
    pub cpg_unmeth: u32,
    pub chg_meth: u32,
    pub chg_unmeth: u32,
    pub chh_meth: u32,
    pub chh_unmeth: u32,
    pub unknown_meth: u32,
    pub unknown_unmeth: u32,
}

impl MethylationCounts {
    pub fn total_non_cg_meth(&self) -> u32 {
        self.chg_meth + self.chh_meth
    }
    pub fn total_non_cg_unmeth(&self) -> u32 {
        self.chg_unmeth + self.chh_unmeth
    }
}

/// A parsed SAM record with Bismark-specific auxiliary tags.
#[derive(Debug, Clone)]
pub struct SamRecord {
    /// Raw line bytes (includes the trailing newline stripped)
    pub raw: Vec<u8>,
    pub qname: Vec<u8>,
    pub flag: u16,
    pub rname: Vec<u8>,
    /// 1-based alignment position
    pub pos: u32,
    pub mapq: u8,
    pub cigar: Vec<u8>,
    pub rnext: Vec<u8>,
    pub pnext: u32,
    pub tlen: i32,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    /// XM:Z methylation call string
    pub xm: Option<Vec<u8>>,
    pub xr: Option<ReadConversion>,
    pub xg: Option<GenomeConversion>,
    /// NM:i edit distance
    pub nm: Option<u32>,
}

impl SamRecord {
    /// Parse a SAM data line (not a header line).
    pub fn parse(line: &[u8]) -> Result<Self, SamError> {
        let fields: Vec<&[u8]> = line.splitn(12, |&b| b == b'\t').collect();
        if fields.len() < 11 {
            return Err(SamError::Malformed(
                String::from_utf8_lossy(line).into_owned(),
            ));
        }

        let flag = parse_u16(fields[1])?;
        let pos = parse_u32(fields[3])?;
        let mapq = parse_u8(fields[4])?;
        let pnext = parse_u32(fields[7])?;
        let tlen = parse_i32(fields[8])?;

        let mut xm: Option<Vec<u8>> = None;
        let mut xr: Option<ReadConversion> = None;
        let mut xg: Option<GenomeConversion> = None;
        let mut nm: Option<u32> = None;

        if fields.len() == 12 {
            for tag in fields[11].split(|&b| b == b'\t') {
                if tag.starts_with(b"XM:Z:") {
                    xm = Some(tag[5..].to_vec());
                } else if tag.starts_with(b"XR:Z:") {
                    xr = Some(match &tag[5..] {
                        b"CT" => ReadConversion::CT,
                        b"GA" => ReadConversion::GA,
                        _ => return Err(SamError::Malformed(format!(
                            "unknown XR value: {}",
                            String::from_utf8_lossy(&tag[5..])
                        ))),
                    });
                } else if tag.starts_with(b"XG:Z:") {
                    xg = Some(match &tag[5..] {
                        b"CT" => GenomeConversion::CT,
                        b"GA" => GenomeConversion::GA,
                        _ => return Err(SamError::Malformed(format!(
                            "unknown XG value: {}",
                            String::from_utf8_lossy(&tag[5..])
                        ))),
                    });
                } else if tag.starts_with(b"NM:i:") {
                    nm = parse_u32(&tag[5..]).ok();
                }
            }
        }

        Ok(SamRecord {
            raw: line.to_vec(),
            qname: fields[0].to_vec(),
            flag,
            rname: fields[2].to_vec(),
            pos,
            mapq,
            cigar: fields[5].to_vec(),
            rnext: fields[6].to_vec(),
            pnext,
            tlen,
            seq: fields[9].to_vec(),
            qual: fields[10].to_vec(),
            xm,
            xr,
            xg,
            nm,
        })
    }

    /// Determine Bismark strand from XR + XG tags.
    pub fn bismark_strand(&self) -> Option<BismarkStrand> {
        Some(BismarkStrand::from_tags(self.xr?, self.xg?))
    }

    /// Count methylation calls in each context from the XM string.
    pub fn count_methylation(&self) -> MethylationCounts {
        let mut counts = MethylationCounts::default();
        if let Some(xm) = &self.xm {
            for &b in xm.iter() {
                match b {
                    b'Z' => counts.cpg_meth += 1,
                    b'z' => counts.cpg_unmeth += 1,
                    b'X' => counts.chg_meth += 1,
                    b'x' => counts.chg_unmeth += 1,
                    b'H' => counts.chh_meth += 1,
                    b'h' => counts.chh_unmeth += 1,
                    b'U' => counts.unknown_meth += 1,
                    b'u' => counts.unknown_unmeth += 1,
                    _ => {}
                }
            }
        }
        counts
    }

    /// Compute the 1-based end position of the alignment using CIGAR.
    /// Operations that consume reference: M, D, N, =, X. (I and S do not.)
    pub fn end_pos(&self) -> u32 {
        let mut pos = self.pos;
        let cigar = &self.cigar;
        let mut num_buf = 0u32;
        for &b in cigar.iter() {
            if b.is_ascii_digit() {
                num_buf = num_buf * 10 + (b - b'0') as u32;
            } else {
                match b {
                    b'M' | b'D' | b'N' | b'=' | b'X' => pos += num_buf,
                    _ => {} // I, S, H, P do not consume reference
                }
                num_buf = 0;
            }
        }
        // pos is now one past the last base; last base is at pos-1
        pos.saturating_sub(1)
    }

    /// True if the read is unmapped (FLAG bit 0x4).
    pub fn is_unmapped(&self) -> bool {
        self.flag & 0x4 != 0
    }

    /// True if this is read 2 in a pair (FLAG bit 0x80).
    pub fn is_read2(&self) -> bool {
        self.flag & 0x80 != 0
    }

    /// True if reverse-complemented (FLAG bit 0x10).
    pub fn is_reverse(&self) -> bool {
        self.flag & 0x10 != 0
    }
}

/// Iterator over SAM records from a BufRead source, yielding header lines and
/// data records separately.
pub struct SamReader<R: BufRead> {
    inner: R,
    buf: Vec<u8>,
}

impl<R: BufRead> SamReader<R> {
    pub fn new(inner: R) -> Self {
        SamReader { inner, buf: Vec::with_capacity(4096) }
    }

    /// Read the next line. Returns `None` on EOF.
    pub fn next_line(&mut self) -> Option<Result<SamLine, SamError>> {
        self.buf.clear();
        match self.inner.read_until(b'\n', &mut self.buf) {
            Ok(0) => None,
            Ok(_) => {
                // strip trailing \r\n
                while self.buf.last() == Some(&b'\n') || self.buf.last() == Some(&b'\r') {
                    self.buf.pop();
                }
                if self.buf.starts_with(b"@") {
                    Some(Ok(SamLine::Header(self.buf.clone())))
                } else {
                    Some(SamRecord::parse(&self.buf).map(SamLine::Record))
                }
            }
            Err(e) => Some(Err(SamError::Io(e))),
        }
    }
}

pub enum SamLine {
    Header(Vec<u8>),
    Record(SamRecord),
}

// -- helpers --

fn parse_u8(s: &[u8]) -> Result<u8, SamError> {
    std::str::from_utf8(s)
        .ok()
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| SamError::Malformed(format!("expected u8: {}", String::from_utf8_lossy(s))))
}

fn parse_u16(s: &[u8]) -> Result<u16, SamError> {
    std::str::from_utf8(s)
        .ok()
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| SamError::Malformed(format!("expected u16: {}", String::from_utf8_lossy(s))))
}

fn parse_u32(s: &[u8]) -> Result<u32, SamError> {
    std::str::from_utf8(s)
        .ok()
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| SamError::Malformed(format!("expected u32: {}", String::from_utf8_lossy(s))))
}

fn parse_i32(s: &[u8]) -> Result<i32, SamError> {
    std::str::from_utf8(s)
        .ok()
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| SamError::Malformed(format!("expected i32: {}", String::from_utf8_lossy(s))))
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE_LINE: &[u8] =
        b"read1\t0\tchr1\t100\t255\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tXM:Z:Z.z.X.x.\tXR:Z:CT\tXG:Z:CT\tNM:i:2";

    #[test]
    fn test_parse_fields() {
        let r = SamRecord::parse(SAMPLE_LINE).unwrap();
        assert_eq!(r.flag, 0);
        assert_eq!(r.pos, 100);
        assert_eq!(r.xr, Some(ReadConversion::CT));
        assert_eq!(r.xg, Some(GenomeConversion::CT));
        assert_eq!(r.bismark_strand(), Some(BismarkStrand::OT));
    }

    #[test]
    fn test_count_methylation() {
        let r = SamRecord::parse(SAMPLE_LINE).unwrap();
        let c = r.count_methylation();
        assert_eq!(c.cpg_meth, 1);   // Z
        assert_eq!(c.cpg_unmeth, 1); // z
        assert_eq!(c.chg_meth, 1);   // X
        assert_eq!(c.chg_unmeth, 1); // x
    }

    #[test]
    fn test_end_pos_simple() {
        let mut r = SamRecord::parse(SAMPLE_LINE).unwrap();
        // 10M starting at 100 → end at 109
        assert_eq!(r.end_pos(), 109);

        r.cigar = b"5M2I3M".to_vec();
        r.pos = 1;
        // 5M + 2I(no ref) + 3M = 8 ref bases → end at 8
        assert_eq!(r.end_pos(), 8);

        r.cigar = b"5M2D3M".to_vec();
        r.pos = 1;
        // 5M + 2D + 3M = 10 ref bases → end at 10
        assert_eq!(r.end_pos(), 10);
    }

    #[test]
    fn test_bismark_strand_all() {
        assert_eq!(
            BismarkStrand::from_tags(ReadConversion::CT, GenomeConversion::CT),
            BismarkStrand::OT
        );
        assert_eq!(
            BismarkStrand::from_tags(ReadConversion::GA, GenomeConversion::CT),
            BismarkStrand::CTOT
        );
        assert_eq!(
            BismarkStrand::from_tags(ReadConversion::GA, GenomeConversion::GA),
            BismarkStrand::CTOB
        );
        assert_eq!(
            BismarkStrand::from_tags(ReadConversion::CT, GenomeConversion::GA),
            BismarkStrand::OB
        );
    }
}
