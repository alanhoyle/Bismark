use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use anyhow::{bail, Context, Result};
use bismark_lib::sam::{GenomeConversion, ReadConversion};

use crate::convert::{ct_convert, ga_convert};
use crate::fastq::{write_record, FastqRecord};

/// One of the four bisulfite strand aligners.
pub struct StrandConfig {
    pub name: &'static str,
    /// Conversion applied to R1 (and to SE reads).
    pub read_conv: ReadConversion,
    /// Conversion applied to R2 in PE mode (opposite of R1 for standard bisulfite PE).
    pub read_conv_r2: ReadConversion,
    pub genome_conv: GenomeConversion,
    pub index: PathBuf,
    /// Pass --norc to Bowtie2 (only forward-strand alignments).
    pub norc: bool,
    /// Pass --nofw to Bowtie2 (only reverse-complement alignments).
    pub nofw: bool,
}

/// Build strand configs for directional mode (OT + OB).
pub fn strand_configs_directional(genome_folder: &Path) -> [StrandConfig; 2] {
    let ct_idx = genome_folder.join("Bisulfite_Genome/CT_conversion/BS_CT");
    let ga_idx = genome_folder.join("Bisulfite_Genome/GA_conversion/BS_GA");
    [
        StrandConfig {
            name: "OT",
            read_conv: ReadConversion::CT,
            read_conv_r2: ReadConversion::GA,
            genome_conv: GenomeConversion::CT,
            index: ct_idx,
            norc: true,
            nofw: false,
        },
        StrandConfig {
            name: "OB",
            read_conv: ReadConversion::CT,
            read_conv_r2: ReadConversion::GA,
            genome_conv: GenomeConversion::GA,
            index: ga_idx,
            norc: false,
            nofw: true,
        },
    ]
}

/// Build strand configs for non-directional mode (OT + OB + CTOT + CTOB).
pub fn strand_configs_nondirectional(genome_folder: &Path) -> Vec<StrandConfig> {
    let ct_idx = genome_folder.join("Bisulfite_Genome/CT_conversion/BS_CT");
    let ga_idx = genome_folder.join("Bisulfite_Genome/GA_conversion/BS_GA");
    vec![
        StrandConfig {
            name: "OT",
            read_conv: ReadConversion::CT,
            read_conv_r2: ReadConversion::GA,
            genome_conv: GenomeConversion::CT,
            index: ct_idx.clone(),
            norc: false,
            nofw: false,
        },
        StrandConfig {
            name: "OB",
            read_conv: ReadConversion::CT,
            read_conv_r2: ReadConversion::GA,
            genome_conv: GenomeConversion::GA,
            index: ga_idx.clone(),
            norc: false,
            nofw: false,
        },
        StrandConfig {
            name: "CTOT",
            read_conv: ReadConversion::GA,
            read_conv_r2: ReadConversion::CT,
            genome_conv: GenomeConversion::CT,
            index: ct_idx,
            norc: false,
            nofw: false,
        },
        StrandConfig {
            name: "CTOB",
            read_conv: ReadConversion::GA,
            read_conv_r2: ReadConversion::CT,
            genome_conv: GenomeConversion::GA,
            index: ga_idx,
            norc: false,
            nofw: false,
        },
    ]
}

/// Build strand configs for PBAT mode (CTOT + CTOB, i.e. GA reads to both genomes).
pub fn strand_configs_pbat(genome_folder: &Path) -> [StrandConfig; 2] {
    let ct_idx = genome_folder.join("Bisulfite_Genome/CT_conversion/BS_CT");
    let ga_idx = genome_folder.join("Bisulfite_Genome/GA_conversion/BS_GA");
    [
        StrandConfig {
            name: "CTOT",
            read_conv: ReadConversion::GA,
            read_conv_r2: ReadConversion::CT,
            genome_conv: GenomeConversion::CT,
            index: ct_idx,
            norc: true,
            nofw: false,
        },
        StrandConfig {
            name: "CTOB",
            read_conv: ReadConversion::GA,
            read_conv_r2: ReadConversion::CT,
            genome_conv: GenomeConversion::GA,
            index: ga_idx,
            norc: false,
            nofw: true,
        },
    ]
}

/// A raw alignment hit returned from Bowtie2 (one read/mate).
#[derive(Debug, Clone)]
pub struct RawHit {
    pub qname: Vec<u8>,
    pub flag: u16,
    pub rname: Vec<u8>,
    /// 1-based alignment position.
    pub pos: u32,
    pub mapq: u8,
    pub cigar: Vec<u8>,
    pub rnext: Vec<u8>,
    pub pnext: u32,
    pub tlen: i32,
    /// Sequence as it appears in SAM (may be RC of original).
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    /// AS:i alignment score (absent means unmapped).
    pub as_score: Option<i32>,
    /// XS:i secondary alignment score (present → possible multi-mapper).
    pub xs_score: Option<i32>,
}

impl RawHit {
    /// True if this record is unmapped (FLAG bit 4).
    pub fn is_unmapped(&self) -> bool {
        self.flag & 0x4 != 0
    }
}

/// Run Bowtie2 on SE reads for one strand and return one `Option<RawHit>` per
/// input read (None = unmapped or filtered out).
pub fn align_strand_se(
    bowtie2: &str,
    cfg: &StrandConfig,
    reads: &[FastqRecord],
    extra_args: &[String],
    threads: usize,
) -> Result<Vec<Option<RawHit>>> {
    let mut cmd = bowtie2_base_cmd(bowtie2, cfg, threads);
    cmd.arg("-U").arg("-");
    for a in extra_args {
        cmd.arg(a);
    }
    cmd.stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit());

    let mut child = cmd
        .spawn()
        .with_context(|| format!("failed to spawn bowtie2 for strand {}", cfg.name))?;

    let stdin = child.stdin.take().unwrap();
    let conv = cfg.read_conv;
    let reads_clone: Vec<FastqRecord> = reads.to_vec();

    let writer = std::thread::spawn(move || -> Result<()> {
        let mut w = BufWriter::new(stdin);
        for rec in &reads_clone {
            let converted_seq = convert_seq(&rec.seq, conv);
            let mut converted = rec.clone();
            converted.seq = converted_seq;
            write_record(&mut w, &converted)?;
        }
        Ok(())
    });

    let stdout = child.stdout.take().unwrap();
    let sam_hits = collect_se_hits(BufReader::new(stdout), reads.len())?;

    writer
        .join()
        .expect("writer thread panicked")
        .context("write to bowtie2 stdin")?;
    let status = child.wait()?;
    if !status.success() {
        bail!(
            "bowtie2 exited with status {} for strand {}",
            status,
            cfg.name
        );
    }

    Ok(sam_hits)
}

/// Run Bowtie2 on PE reads for one strand using `--interleaved` stdin.
/// Returns one `Option<(RawHit, RawHit)>` per pair (None = concordant pair not
/// found / either mate unmapped).
pub fn align_strand_pe(
    bowtie2: &str,
    cfg: &StrandConfig,
    pairs: &[(FastqRecord, FastqRecord)],
    extra_args: &[String],
    threads: usize,
) -> Result<Vec<Option<(RawHit, RawHit)>>> {
    let mut cmd = bowtie2_base_cmd(bowtie2, cfg, threads);
    cmd.arg("--interleaved").arg("-");
    for a in extra_args {
        cmd.arg(a);
    }
    cmd.stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit());

    let mut child = cmd
        .spawn()
        .with_context(|| format!("failed to spawn bowtie2 for strand {}", cfg.name))?;

    let stdin = child.stdin.take().unwrap();
    let conv_r1 = cfg.read_conv;
    let conv_r2 = cfg.read_conv_r2;
    let pairs_clone: Vec<(FastqRecord, FastqRecord)> = pairs.to_vec();

    let writer = std::thread::spawn(move || -> Result<()> {
        let mut w = BufWriter::new(stdin);
        for (r1, r2) in &pairs_clone {
            let mut r1c = r1.clone();
            let mut r2c = r2.clone();
            r1c.seq = convert_seq(&r1.seq, conv_r1);
            r2c.seq = convert_seq(&r2.seq, conv_r2);
            write_record(&mut w, &r1c)?;
            write_record(&mut w, &r2c)?;
        }
        Ok(())
    });

    let stdout = child.stdout.take().unwrap();
    let sam_hits = collect_pe_hits(BufReader::new(stdout), pairs.len())?;

    writer
        .join()
        .expect("writer thread panicked")
        .context("write to bowtie2 stdin")?;
    let status = child.wait()?;
    if !status.success() {
        bail!(
            "bowtie2 exited with status {} for strand {}",
            status,
            cfg.name
        );
    }

    Ok(sam_hits)
}

// ── internal helpers ──────────────────────────────────────────────────────────

fn bowtie2_base_cmd(bowtie2: &str, cfg: &StrandConfig, threads: usize) -> Command {
    let mut cmd = Command::new(bowtie2);
    cmd.arg("-p")
        .arg(threads.to_string())
        .arg("-x")
        .arg(&cfg.index);
    if cfg.norc {
        cmd.arg("--norc");
    }
    if cfg.nofw {
        cmd.arg("--nofw");
    }
    cmd
}

fn convert_seq(seq: &[u8], conv: ReadConversion) -> Vec<u8> {
    match conv {
        ReadConversion::CT => ct_convert(seq),
        ReadConversion::GA => ga_convert(seq),
    }
}

/// Collect SAM records for SE mode.  Bowtie2 emits one record per input read
/// (when `--no-unal` is NOT set — we want unmapped so we can count them).
/// We key them by QNAME order in the SAM output, which matches input order.
fn collect_se_hits<R: BufRead>(reader: R, n_reads: usize) -> Result<Vec<Option<RawHit>>> {
    let mut out: Vec<Option<RawHit>> = vec![None; n_reads];
    let mut idx = 0usize;
    for line in reader.split(b'\n') {
        let line = line?;
        if line.starts_with(b"@") {
            continue;
        }
        if line.is_empty() {
            continue;
        }
        if idx >= n_reads {
            break;
        }
        let hit = parse_sam_line(&line)?;
        if !hit.is_unmapped() {
            out[idx] = Some(hit);
        }
        idx += 1;
    }
    Ok(out)
}

/// Collect SAM records for PE mode.  Bowtie2 --interleaved emits records in
/// pairs; we group them back into (R1, R2) tuples, keeping only concordant
/// aligned pairs where both mates mapped.
fn collect_pe_hits<R: BufRead>(reader: R, n_pairs: usize) -> Result<Vec<Option<(RawHit, RawHit)>>> {
    let mut out: Vec<Option<(RawHit, RawHit)>> = vec![None; n_pairs];
    let mut pending: Option<RawHit> = None;
    let mut pair_idx = 0usize;

    for line in reader.split(b'\n') {
        let line = line?;
        if line.starts_with(b"@") {
            continue;
        }
        if line.is_empty() {
            continue;
        }
        let hit = parse_sam_line(&line)?;
        match pending.take() {
            None => {
                pending = Some(hit);
            }
            Some(r1) => {
                if pair_idx >= n_pairs {
                    break;
                }
                // Both mates must be mapped and concordant (FLAG & 0x2 = proper pair).
                let both_mapped = !r1.is_unmapped() && !hit.is_unmapped();
                let concordant = r1.flag & 0x2 != 0 && hit.flag & 0x2 != 0;
                if both_mapped && concordant {
                    out[pair_idx] = Some((r1, hit));
                }
                pair_idx += 1;
            }
        }
    }
    Ok(out)
}

/// Parse a single SAM data line.
fn parse_sam_line(line: &[u8]) -> Result<RawHit> {
    let fields: Vec<&[u8]> = line.splitn(12, |&b| b == b'\t').collect();
    if fields.len() < 11 {
        bail!("malformed SAM line: {}", String::from_utf8_lossy(line));
    }

    let flag = parse_int::<u16>(fields[1])?;
    let pos = parse_int::<u32>(fields[3])?;
    let mapq = parse_int::<u8>(fields[4])?;
    let pnext = parse_int::<u32>(fields[7])?;
    let tlen = parse_signed::<i32>(fields[8])?;

    let mut as_score: Option<i32> = None;
    let mut xs_score: Option<i32> = None;
    if fields.len() == 12 {
        for tag in fields[11].split(|&b| b == b'\t') {
            if let Some(v) = parse_i_tag(b"AS", tag) {
                as_score = Some(v);
            } else if let Some(v) = parse_i_tag(b"XS", tag) {
                xs_score = Some(v);
            }
        }
    }

    // Treat unmapped reads (flag & 4) as having no alignment score.
    if flag & 0x4 != 0 {
        as_score = None;
    }

    Ok(RawHit {
        qname: fields[0].to_vec(),
        flag,
        rname: normalise_bismark_rname(fields[2]),
        pos,
        mapq,
        cigar: fields[5].to_vec(),
        rnext: normalise_bismark_rname(fields[6]),
        pnext,
        tlen,
        seq: fields[9].to_vec(),
        qual: fields[10].to_vec(),
        as_score,
        xs_score,
    })
}

fn normalise_bismark_rname(name: &[u8]) -> Vec<u8> {
    for suffix in [b"_CT_converted".as_slice(), b"_GA_converted".as_slice()] {
        if let Some(prefix) = name.strip_suffix(suffix) {
            return prefix.to_vec();
        }
    }
    name.to_vec()
}

fn parse_i_tag(name: &[u8], tag: &[u8]) -> Option<i32> {
    if tag.len() > 5 && &tag[0..2] == name && &tag[2..5] == b":i:" {
        std::str::from_utf8(&tag[5..])
            .ok()
            .and_then(|s| s.parse().ok())
    } else {
        None
    }
}

fn parse_int<T: std::str::FromStr>(b: &[u8]) -> Result<T>
where
    T::Err: std::fmt::Display,
{
    std::str::from_utf8(b)
        .context("non-UTF8 SAM field")?
        .parse::<T>()
        .map_err(|e| anyhow::anyhow!("parse error: {}", e))
}

fn parse_signed<T: std::str::FromStr>(b: &[u8]) -> Result<T>
where
    T::Err: std::fmt::Display,
{
    parse_int(b)
}

#[cfg(test)]
mod tests {
    use super::normalise_bismark_rname;

    #[test]
    fn normalises_converted_reference_names() {
        assert_eq!(
            normalise_bismark_rname(b"Ecoli_K12_CT_converted"),
            b"Ecoli_K12"
        );
        assert_eq!(
            normalise_bismark_rname(b"Ecoli_K12_GA_converted"),
            b"Ecoli_K12"
        );
        assert_eq!(normalise_bismark_rname(b"Ecoli_K12"), b"Ecoli_K12");
        assert_eq!(normalise_bismark_rname(b"*"), b"*");
        assert_eq!(normalise_bismark_rname(b"="), b"=");
    }
}
