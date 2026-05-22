/// Convert all C → T in `seq` (bisulfite top-strand conversion).
pub fn ct_convert(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| if b == b'C' || b == b'c' { b'T' } else { b })
        .collect()
}

/// Convert all G → A in `seq` (bisulfite bottom-strand conversion).
pub fn ga_convert(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| if b == b'G' || b == b'g' { b'A' } else { b })
        .collect()
}

#[allow(dead_code)]
/// Reverse-complement a DNA sequence (uppercase).
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

#[inline(always)]
pub fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'N' | b'n' => b'N',
        _ => b'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ct() {
        assert_eq!(ct_convert(b"ACGT"), b"ATGT");
    }

    #[test]
    fn ga() {
        assert_eq!(ga_convert(b"ACGT"), b"ACAT");
    }

    #[test]
    fn rc() {
        assert_eq!(revcomp(b"ACGT"), b"ACGT");
        assert_eq!(revcomp(b"AACG"), b"CGTT");
    }
}
