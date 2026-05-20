use crate::sam::ReadConversion;

/// Produce the Bismark XM character for one base position.
///
/// For CT reads (forward strand): `observed` is the read base, `genomic` is a
/// slice of at least 3 bytes starting at this position (to allow +1 and +2
/// lookahead for context determination).
///
/// For GA reads (bottom strand): the genomic slice is shifted by 2 relative to
/// the read index (as in the Perl source) — caller must pass
/// `&genomic_slice[index..]` where `genomic_slice` has the 2-base prefix.
///
/// Returns one of: `.`, `Z`, `z`, `X`, `x`, `H`, `h`, `U`, `u`
#[inline(always)]
pub fn methylation_char(observed: u8, genomic: &[u8], conv: ReadConversion) -> u8 {
    debug_assert!(genomic.len() >= 3);
    match conv {
        ReadConversion::CT => {
            if genomic[0] == b'C' {
                if observed == b'C' {
                    // Methylated — protected from conversion
                    meth_context_ct(genomic[1], genomic[2], true)
                } else if observed == b'T' {
                    // Unmethylated — converted C→T
                    meth_context_ct(genomic[1], genomic[2], false)
                } else {
                    b'.'
                }
            } else {
                b'.'
            }
        }
        ReadConversion::GA => {
            // In GA reads, genomic[2] is the base at this read position (the G that
            // corresponds to a C on the opposite strand). genomic[1] is +1 upstream,
            // genomic[0] is +2 upstream (matching Perl's $genomic[$index], $genomic[$index+1],
            // $genomic[$index+2] with the 2-offset shift).
            if genomic[2] == b'G' {
                if observed == b'G' {
                    meth_context_ga(genomic[1], genomic[0], true)
                } else if observed == b'A' {
                    meth_context_ga(genomic[1], genomic[0], false)
                } else {
                    b'.'
                }
            } else {
                b'.'
            }
        }
    }
}

/// Determine context for CT-strand cytosine given downstream bases.
#[inline(always)]
fn meth_context_ct(next: u8, next2: u8, methylated: bool) -> u8 {
    if next == b'G' {
        if methylated { b'Z' } else { b'z' }
    } else if next == b'N' || next == b'X' {
        if methylated { b'U' } else { b'u' }
    } else if next2 == b'G' {
        if methylated { b'X' } else { b'x' }
    } else if next2 == b'N' || next2 == b'X' {
        if methylated { b'U' } else { b'u' }
    } else {
        if methylated { b'H' } else { b'h' }
    }
}

/// Determine context for GA-strand guanine given upstream bases.
/// `up1` = genomic[index+1] (one base upstream), `up2` = genomic[index] (two upstream).
#[inline(always)]
fn meth_context_ga(up1: u8, up2: u8, methylated: bool) -> u8 {
    if up1 == b'C' {
        if methylated { b'Z' } else { b'z' }
    } else if up1 == b'N' || up1 == b'X' {
        if methylated { b'U' } else { b'u' }
    } else if up2 == b'C' {
        if methylated { b'X' } else { b'x' }
    } else if up2 == b'N' || up2 == b'X' {
        if methylated { b'U' } else { b'u' }
    } else {
        if methylated { b'H' } else { b'h' }
    }
}

/// Build a complete XM methylation string from a read sequence, genomic
/// sequence, and read conversion direction.
///
/// `seq` and `genomic` must satisfy: `genomic.len() == seq.len() + 2`
/// (the 2-base lookahead / lookbehind padding used in Bismark).
pub fn build_xm_string(seq: &[u8], genomic: &[u8], conv: ReadConversion) -> Vec<u8> {
    assert_eq!(
        genomic.len(),
        seq.len() + 2,
        "genomic must be seq.len()+2 bytes"
    );
    let mut out = Vec::with_capacity(seq.len());
    match conv {
        ReadConversion::CT => {
            for i in 0..seq.len() {
                out.push(methylation_char(seq[i], &genomic[i..i + 3], conv));
            }
        }
        ReadConversion::GA => {
            for i in 0..seq.len() {
                // For GA reads, pass the slice starting at index i (which has
                // genomic[i], genomic[i+1], genomic[i+2]) and methylation_char
                // uses [2] as the base, [1] as upstream1, [0] as upstream2.
                out.push(methylation_char(seq[i], &genomic[i..i + 3], conv));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::ReadConversion;

    // CT strand tests
    #[test]
    fn ct_cpg_methylated() {
        // Genomic: CGx  (C followed by G → CpG context), observed C (protected)
        assert_eq!(methylation_char(b'C', b"CGT", ReadConversion::CT), b'Z');
    }

    #[test]
    fn ct_cpg_unmethylated() {
        assert_eq!(methylation_char(b'T', b"CGT", ReadConversion::CT), b'z');
    }

    #[test]
    fn ct_chg_methylated() {
        // Genomic: CAG → not CpG, second downstream = G → CHG
        assert_eq!(methylation_char(b'C', b"CAG", ReadConversion::CT), b'X');
    }

    #[test]
    fn ct_chg_unmethylated() {
        assert_eq!(methylation_char(b'T', b"CAG", ReadConversion::CT), b'x');
    }

    #[test]
    fn ct_chh_methylated() {
        // Genomic: CAT → not CpG, second downstream not G → CHH
        assert_eq!(methylation_char(b'C', b"CAT", ReadConversion::CT), b'H');
    }

    #[test]
    fn ct_chh_unmethylated() {
        assert_eq!(methylation_char(b'T', b"CAT", ReadConversion::CT), b'h');
    }

    #[test]
    fn ct_unknown_n() {
        // N downstream → unknown context
        assert_eq!(methylation_char(b'C', b"CNT", ReadConversion::CT), b'U');
        assert_eq!(methylation_char(b'T', b"CNT", ReadConversion::CT), b'u');
    }

    #[test]
    fn ct_non_cytosine() {
        assert_eq!(methylation_char(b'A', b"AGT", ReadConversion::CT), b'.');
        assert_eq!(methylation_char(b'G', b"GCG", ReadConversion::CT), b'.');
    }

    // GA strand tests
    #[test]
    fn ga_cpg_methylated() {
        // genomic[2]=G, genomic[1]=C → CpG context, observed G (protected)
        assert_eq!(methylation_char(b'G', b"ACG", ReadConversion::GA), b'Z');
    }

    #[test]
    fn ga_cpg_unmethylated() {
        assert_eq!(methylation_char(b'A', b"ACG", ReadConversion::GA), b'z');
    }

    #[test]
    fn ga_chg_methylated() {
        // genomic[2]=G, genomic[1]=T (not C), genomic[0]=C → CHG
        assert_eq!(methylation_char(b'G', b"CTG", ReadConversion::GA), b'X');
    }

    #[test]
    fn ga_chh_methylated() {
        // genomic[2]=G, genomic[1]=T, genomic[0]=A → CHH
        assert_eq!(methylation_char(b'G', b"ATG", ReadConversion::GA), b'H');
    }

    #[test]
    fn ga_non_guanine() {
        assert_eq!(methylation_char(b'A', b"ATA", ReadConversion::GA), b'.');
    }
}
