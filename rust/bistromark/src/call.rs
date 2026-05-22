use bismark_lib::fasta::Genome;
use bismark_lib::methylation::methylation_char;
use bismark_lib::sam::{GenomeConversion, ReadConversion};

use crate::align::{RawHit, StrandConfig};
use crate::convert::revcomp;

#[allow(dead_code)]
/// A fully resolved alignment ready for BAM output.
#[derive(Debug, Clone)]
pub struct AlignedRead {
    pub qname: Vec<u8>,
    pub flag: u16,
    pub rname: Vec<u8>,
    /// 1-based position.
    pub pos: u32,
    pub mapq: u8,
    pub cigar: Vec<u8>,
    pub rnext: Vec<u8>,
    pub pnext: u32,
    pub tlen: i32,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    /// XM methylation string (same length as seq after soft-clipping).
    pub xm: Vec<u8>,
    pub xr: ReadConversion,
    pub xg: GenomeConversion,
    pub nm: u32,
    pub as_score: i32,
    pub xs_score: Option<i32>,
    /// All extra tags from Bowtie2 (NM, MD, AS, XS, …).
    pub extra_tags: Vec<u8>,
    pub outcome: Outcome,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    Unique,
    Ambiguous,
    Unmapped,
}

fn has_unique_bowtie_best(hit: &RawHit) -> bool {
    match (hit.as_score, hit.xs_score) {
        (Some(as_score), Some(xs_score)) => xs_score < as_score,
        (Some(_), None) => true,
        _ => false,
    }
}

/// Select the best SE alignment across all strands and compute XM tags.
///
/// `strand_hits[s][i]` = alignment of read i from strand s (None = unmapped).
pub fn select_best_se(
    strand_hits: &[Vec<Option<RawHit>>],
    configs: &[StrandConfig],
    reads_seq: &[Vec<u8>],
    reads_qual: &[Vec<u8>],
    genome: &Genome,
) -> Vec<AlignedRead> {
    let n = strand_hits[0].len();
    let mut out = Vec::with_capacity(n);

    for i in 0..n {
        // Collect strands that produced a mapped hit.
        let mapped: Vec<(usize, &RawHit)> = strand_hits
            .iter()
            .enumerate()
            .filter_map(|(s, hits)| hits[i].as_ref().map(|h| (s, h)))
            .collect();

        if mapped.is_empty() {
            out.push(unmapped_read(&reads_seq[i]));
            continue;
        }

        // Find the maximum AS score.
        let best_score = mapped
            .iter()
            .filter_map(|(_, h)| h.as_score)
            .max()
            .unwrap_or(i32::MIN);

        let best_hits: Vec<(usize, &RawHit)> = mapped
            .iter()
            .filter(|(_, h)| h.as_score == Some(best_score))
            .copied()
            .collect();

        let outcome = if best_hits.len() == 1 && has_unique_bowtie_best(best_hits[0].1) {
            Outcome::Unique
        } else {
            Outcome::Ambiguous
        };

        let (strand_idx, hit) = best_hits[0];
        let cfg = &configs[strand_idx];

        let (seq, qual) = orient_read(&reads_seq[i], &reads_qual[i], hit.flag);
        let mut restored_hit = hit.clone();
        restored_hit.seq = seq.clone();
        restored_hit.qual = qual.clone();
        let xm = build_xm(
            &restored_hit,
            effective_conversion(cfg.read_conv, cfg.genome_conv),
            genome,
        );
        let (nm, md) = build_nm_md(&restored_hit, &seq, genome);

        out.push(AlignedRead {
            qname: hit.qname.clone(),
            flag: hit.flag,
            rname: hit.rname.clone(),
            pos: hit.pos,
            mapq: hit.mapq,
            cigar: hit.cigar.clone(),
            rnext: hit.rnext.clone(),
            pnext: hit.pnext,
            tlen: hit.tlen,
            seq,
            qual,
            xm,
            xr: cfg.read_conv,
            xg: cfg.genome_conv,
            nm,
            as_score: best_score,
            xs_score: hit.xs_score,
            extra_tags: format!("NM:i:{nm}\tMD:Z:{}", String::from_utf8_lossy(&md)).into_bytes(),
            outcome,
        });
    }

    out
}

/// Select the best PE alignment across all strands and compute XM tags.
///
/// `strand_hits[s][i]` = concordant pair for read-pair i from strand s.
pub fn select_best_pe(
    strand_hits: &[Vec<Option<(RawHit, RawHit)>>],
    configs: &[StrandConfig],
    r1_seqs: &[Vec<u8>],
    r2_seqs: &[Vec<u8>],
    r1_quals: &[Vec<u8>],
    r2_quals: &[Vec<u8>],
    genome: &Genome,
) -> Vec<(AlignedRead, AlignedRead)> {
    let n = strand_hits[0].len();
    let mut out = Vec::with_capacity(n);

    for i in 0..n {
        let mapped: Vec<(usize, &(RawHit, RawHit))> = strand_hits
            .iter()
            .enumerate()
            .filter_map(|(s, hits)| hits[i].as_ref().map(|p| (s, p)))
            .collect();

        if mapped.is_empty() {
            out.push((unmapped_read(&r1_seqs[i]), unmapped_read(&r2_seqs[i])));
            continue;
        }

        // Score a PE pair as the sum of R1+R2 AS scores.
        let pair_score = |r1: &RawHit, r2: &RawHit| -> i32 {
            r1.as_score.unwrap_or(i32::MIN / 2) + r2.as_score.unwrap_or(i32::MIN / 2)
        };

        let best_score = mapped
            .iter()
            .map(|(_, (r1, r2))| pair_score(r1, r2))
            .max()
            .unwrap_or(i32::MIN);

        let best_hits: Vec<(usize, &(RawHit, RawHit))> = mapped
            .iter()
            .filter(|(_, (r1, r2))| pair_score(r1, r2) == best_score)
            .copied()
            .collect();

        let unique_pair = if best_hits.len() == 1 {
            let (r1, r2) = best_hits[0].1;
            let r1_unique = has_unique_bowtie_best(r1);
            let r2_unique = has_unique_bowtie_best(r2);
            (r1_unique && r2_unique)
                || ((r1_unique || r2_unique) && r1.mapq >= 30 && r2.mapq >= 30)
        } else {
            false
        };
        let outcome = if unique_pair {
            Outcome::Unique
        } else {
            Outcome::Ambiguous
        };

        let (strand_idx, (r1, r2)) = best_hits[0];
        let cfg = &configs[strand_idx];

        let (seq1, qual1) = orient_read(&r1_seqs[i], &r1_quals[i], r1.flag);
        let (seq2, qual2) = orient_read(&r2_seqs[i], &r2_quals[i], r2.flag);
        let mut restored_r1 = r1.clone();
        restored_r1.seq = seq1.clone();
        restored_r1.qual = qual1.clone();
        let mut restored_r2 = r2.clone();
        restored_r2.seq = seq2.clone();
        restored_r2.qual = qual2.clone();

        let xm1 = build_xm(
            &restored_r1,
            effective_conversion(cfg.read_conv, cfg.genome_conv),
            genome,
        );
        let xm2 = build_xm(
            &restored_r2,
            effective_conversion(cfg.read_conv_r2, cfg.genome_conv),
            genome,
        );
        let (nm1, md1) = build_nm_md(&restored_r1, &seq1, genome);
        let (nm2, md2) = build_nm_md(&restored_r2, &seq2, genome);

        let make = |hit: &RawHit,
                    qname: &[u8],
                    seq: Vec<u8>,
                    qual: Vec<u8>,
                    xm: Vec<u8>,
                    xr: ReadConversion,
                    nm: u32,
                    md: Vec<u8>|
         -> AlignedRead {
            AlignedRead {
                qname: qname.to_vec(),
                flag: hit.flag,
                rname: hit.rname.clone(),
                pos: hit.pos,
                mapq: hit.mapq,
                cigar: hit.cigar.clone(),
                rnext: hit.rnext.clone(),
                pnext: hit.pnext,
                tlen: hit.tlen,
                seq,
                qual,
                xm,
                xr,
                xg: cfg.genome_conv,
                nm,
                as_score: hit.as_score.unwrap_or(0),
                xs_score: hit.xs_score,
                extra_tags: format!("NM:i:{nm}\tMD:Z:{}", String::from_utf8_lossy(&md))
                    .into_bytes(),
                outcome,
            }
        };

        out.push((
            make(r1, &r1.qname, seq1, qual1, xm1, cfg.read_conv, nm1, md1),
            make(r2, &r1.qname, seq2, qual2, xm2, cfg.read_conv_r2, nm2, md2),
        ));
    }

    out
}

// ── XM generation ────────────────────────────────────────────────────────────

fn orient_read(seq: &[u8], qual: &[u8], flag: u16) -> (Vec<u8>, Vec<u8>) {
    if flag & 0x10 != 0 {
        let mut q = qual.to_vec();
        q.reverse();
        (revcomp(seq), q)
    } else {
        (seq.to_vec(), qual.to_vec())
    }
}

fn effective_conversion(
    _read_conv: ReadConversion,
    genome_conv: GenomeConversion,
) -> ReadConversion {
    match genome_conv {
        GenomeConversion::CT => ReadConversion::CT,
        GenomeConversion::GA => ReadConversion::GA,
    }
}

fn build_nm_md(hit: &RawHit, seq: &[u8], genome: &Genome) -> (u32, Vec<u8>) {
    if hit.is_unmapped() || hit.rname == b"*" || hit.cigar == b"*" {
        return (0, b"*".to_vec());
    }

    let chr = match std::str::from_utf8(&hit.rname) {
        Ok(s) => s,
        Err(_) => return (0, b"*".to_vec()),
    };

    let ops = parse_cigar(&hit.cigar);
    let mut read_pos = 0usize;
    let mut genome_pos = hit.pos.saturating_sub(1) as usize;
    let mut nm = 0u32;
    let mut md = Vec::new();
    let mut matches = 0u32;

    for (len, op) in ops {
        let len = len as usize;
        match op {
            b'M' | b'=' | b'X' => {
                for _ in 0..len {
                    let rb = seq
                        .get(read_pos)
                        .copied()
                        .unwrap_or(b'N')
                        .to_ascii_uppercase();
                    let gb = genome.base(chr, genome_pos).to_ascii_uppercase();
                    if rb == gb {
                        matches += 1;
                    } else {
                        md.extend_from_slice(matches.to_string().as_bytes());
                        matches = 0;
                        md.push(gb);
                        nm += 1;
                    }
                    read_pos += 1;
                    genome_pos += 1;
                }
            }
            b'I' => {
                read_pos += len;
                nm += len as u32;
            }
            b'D' => {
                md.extend_from_slice(matches.to_string().as_bytes());
                matches = 0;
                md.push(b'^');
                for _ in 0..len {
                    md.push(genome.base(chr, genome_pos).to_ascii_uppercase());
                    genome_pos += 1;
                    nm += 1;
                }
            }
            b'N' => {
                genome_pos += len;
            }
            b'S' => {
                read_pos += len;
            }
            b'H' | b'P' => {}
            _ => {}
        }
    }
    md.extend_from_slice(matches.to_string().as_bytes());
    (nm, md)
}

/// Build the XM methylation string for a single aligned read.
///
/// The SAM record from Bowtie2 contains the read sequence as aligned (RC for
/// minus-strand reads) and the 1-based leftmost position `pos`.  We walk the
/// CIGAR to pair each read base with its genomic position and call
/// `methylation_char` with the correct 3-byte context window.
pub fn build_xm(hit: &RawHit, conv: ReadConversion, genome: &Genome) -> Vec<u8> {
    if hit.is_unmapped() || hit.rname == b"*" || hit.cigar == b"*" {
        return vec![b'.'; hit.seq.len()];
    }

    let chr = match std::str::from_utf8(&hit.rname) {
        Ok(s) => s,
        Err(_) => return vec![b'.'; hit.seq.len()],
    };

    let ops = parse_cigar(&hit.cigar);

    let mut xm = Vec::with_capacity(hit.seq.len());
    let mut read_pos: usize = 0;
    // genome_pos is 0-based
    let mut genome_pos: usize = hit.pos.saturating_sub(1) as usize;

    for (len, op) in &ops {
        let len = *len as usize;
        match op {
            b'M' | b'X' | b'=' => {
                for _ in 0..len {
                    let ctx = context3(genome, chr, genome_pos, conv);
                    let base = hit.seq[read_pos];
                    xm.push(methylation_char(base, &ctx, conv));
                    read_pos += 1;
                    genome_pos += 1;
                }
            }
            b'I' => {
                // Inserted bases have no genomic position → dot.
                for _ in 0..len {
                    xm.push(b'.');
                    read_pos += 1;
                }
            }
            b'D' | b'N' => {
                // Deleted / skipped bases: reference advances, read does not.
                genome_pos += len;
            }
            b'S' => {
                // Soft-clipped read bases → dot.
                for _ in 0..len {
                    xm.push(b'.');
                    read_pos += 1;
                }
            }
            b'H' | b'P' => {
                // Hard clip / padding: consume nothing.
            }
            _ => {
                // Unknown op: treat as soft clip.
                for _ in 0..len {
                    xm.push(b'.');
                    read_pos += 1;
                }
            }
        }
    }

    // Pad to read length if CIGAR was shorter (shouldn't happen with valid input).
    while xm.len() < hit.seq.len() {
        xm.push(b'.');
    }
    xm.truncate(hit.seq.len());
    xm
}

/// Return the 3-byte context slice for `methylation_char` at the given
/// 0-based genome position.
///
/// Strand convention:
/// - CT / forward (OT): context = [P, P+1, P+2] on plus strand
/// - CT / reverse (OB): context = [comp(P), comp(P-1), comp(P-2)] (minus strand 5'→3')
/// - GA / either  (CTOT or CTOB): context = [P-2, P-1, P] on plus strand
fn context3(
    genome: &Genome,
    chr: &str,
    pos: usize, // 0-based genome position
    conv: ReadConversion,
) -> [u8; 3] {
    match conv {
        ReadConversion::CT => {
            // OT: downstream = increasing genomic coordinate
            [
                genome.base(chr, pos),
                genome.base(chr, pos + 1),
                genome.base(chr, pos + 2),
            ]
        }
        ReadConversion::GA => {
            // CTOT / CTOB: upstream on plus strand (P-2, P-1, P)
            [
                genome.base(chr, pos.saturating_sub(2)),
                genome.base(chr, pos.saturating_sub(1)),
                genome.base(chr, pos),
            ]
        }
    }
}

// ── CIGAR parsing ─────────────────────────────────────────────────────────────

pub fn parse_cigar(cigar: &[u8]) -> Vec<(u32, u8)> {
    let mut ops = Vec::new();
    let mut num: u32 = 0;
    for &b in cigar {
        if b.is_ascii_digit() {
            num = num * 10 + (b - b'0') as u32;
        } else {
            ops.push((num, b));
            num = 0;
        }
    }
    ops
}

// ── helpers ───────────────────────────────────────────────────────────────────

fn unmapped_read(seq: &[u8]) -> AlignedRead {
    AlignedRead {
        qname: b"*".to_vec(),
        flag: 4,
        rname: b"*".to_vec(),
        pos: 0,
        mapq: 0,
        cigar: b"*".to_vec(),
        rnext: b"*".to_vec(),
        pnext: 0,
        tlen: 0,
        seq: seq.to_vec(),
        qual: b"*".to_vec(),
        xm: vec![b'.'; seq.len()],
        xr: ReadConversion::CT,
        xg: GenomeConversion::CT,
        nm: 0,
        as_score: 0,
        xs_score: None,
        extra_tags: Vec::new(),
        outcome: Outcome::Unmapped,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_cigar_basic() {
        assert_eq!(parse_cigar(b"50M"), vec![(50, b'M')]);
        assert_eq!(
            parse_cigar(b"10M2I38M"),
            vec![(10, b'M'), (2, b'I'), (38, b'M')]
        );
        assert_eq!(parse_cigar(b"5S45M"), vec![(5, b'S'), (45, b'M')]);
    }
}
