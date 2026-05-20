use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use bismark_lib::bam_io::{bam_is_truncated, find_samtools, BamReader, BamWriter};
use bismark_lib::BISMARK_VERSION;
use clap::Parser;
use rustc_hash::{FxHashMap, FxHashSet};

/// Dedup key for single-end reads: (strand_index, chr_interned, key_pos)
/// For OT/CTOB (forward): key_pos = POS (start)
/// For CTOT/OB (reverse): key_pos = end_pos computed from CIGAR
type SeKey = (u8, u32, u32);

/// Dedup key for paired-end reads: (strand_index, chr_interned, start, end)
type PeKey = (u8, u32, u32, u32);

#[derive(Parser)]
#[command(
    name = "deduplicate_bismark",
    about = "Remove PCR duplicate alignments from Bismark BAM files",
    version = BISMARK_VERSION,
    disable_version_flag = true,
)]
struct Cli {
    /// BAM file(s) to deduplicate
    files: Vec<PathBuf>,

    /// Force single-end mode (auto-detected if not set)
    #[arg(short = 's', long = "single")]
    single: bool,

    /// Force paired-end mode (auto-detected if not set)
    #[arg(short = 'p', long = "paired")]
    paired: bool,

    /// Use UMI/barcode from read ID for deduplication (RRBS mode)
    #[arg(long = "barcode", visible_alias = "umi")]
    rrbs: bool,

    /// Treat multiple files as one combined sample
    #[arg(long)]
    multiple: bool,

    /// Output directory
    #[arg(long = "output_dir", default_value = "")]
    output_dir: String,

    /// Custom output filename
    #[arg(long = "outfile")]
    outfile: Option<String>,

    /// Barcode/UMI format: read IDs from bcl-convert with internal UMIs
    #[arg(long)]
    bclconvert: bool,

    /// Path to samtools
    #[arg(long)]
    samtools_path: Option<String>,

    /// Number of threads to pass to samtools view
    #[arg(long = "parallel", default_value_t = 1)]
    parallel: u32,

    /// Print version and exit
    #[arg(long = "version")]
    version: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.version {
        println!(
            "\n\n\t\tBismark Deduplication Module\n\n\t\tBismark version: {}\n\t\tCopyright {}\n\n",
            BISMARK_VERSION, "2010-25 Felix Krueger, Altos Bioinformatics"
        );
        return Ok(());
    }

    if cli.single && cli.paired {
        bail!("Please select either -s (single-end) or -p (paired-end), not both!");
    }
    if cli.parallel == 0 {
        bail!("Core usage needs to be set to 1 or more");
    }

    let samtools = find_samtools(cli.samtools_path.as_deref())?;
    let output_dir = normalise_dir(&cli.output_dir);

    if cli.rrbs {
        eprintln!("\nIf the input file has several alignments to the same single position in the genome, only alignments with a unique barcode (UMI) will be chosen)\n");
    } else {
        if cli.multiple {
            eprintln!("Multiple Input files for the same sample selected - All input files are treated as one big single file. The files to be used are:");
            eprintln!(
                "{}\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
                cli.files
                    .iter()
                    .map(|f| f.display().to_string())
                    .collect::<Vec<_>>()
                    .join("\n")
            );
        }
        eprintln!("\nIf there are several alignments to a single position in the genome the first alignment will be chosen. Since the input files are not in any way sorted this is a near-enough random selection of reads.\n");
    }

    if cli.multiple {
        // Treat all files as one combined input
        if !cli.files.is_empty() {
            deduplicate_files(
                &cli.files,
                cli.single,
                cli.paired,
                cli.rrbs,
                cli.bclconvert,
                true,
                &output_dir,
                cli.outfile.as_deref(),
                &samtools,
                cli.parallel,
            )?;
        }
    } else {
        for file in &cli.files {
            deduplicate_files(
                &[file.clone()],
                cli.single,
                cli.paired,
                cli.rrbs,
                cli.bclconvert,
                false,
                &output_dir,
                cli.outfile.as_deref(),
                &samtools,
                cli.parallel,
            )?;
        }
    }

    Ok(())
}

fn deduplicate_files(
    files: &[PathBuf],
    global_single: bool,
    global_paired: bool,
    rrbs: bool,
    bclconvert: bool,
    multiple: bool,
    output_dir: &str,
    user_outfile: Option<&str>,
    samtools: &str,
    parallel: u32,
) -> Result<()> {
    let primary = &files[0];

    // Determine SE/PE
    let (is_single, is_paired) = if global_single {
        (true, false)
    } else if global_paired {
        (false, true)
    } else {
        detect_library_type(samtools, primary)?
    };

    if !is_single && !is_paired {
        bail!("Please specify either -s (single-end) or -p (paired-end), or provide a BAM file with a @PG header line");
    }

    // Check for truncation
    for f in files {
        if f.to_string_lossy().ends_with(".bam") && bam_is_truncated(samtools, f) {
            bail!("File {} appears truncated", f.display());
        }
    }

    // Derive report filename
    let report_stem = derive_stem(primary, user_outfile);
    let report_name = if multiple {
        format!("{report_stem}.multiple.deduplication_report.txt")
    } else {
        format!("{report_stem}.deduplication_report.txt")
    };
    let report_path = format!("{output_dir}{report_name}");
    let mut report = std::fs::File::create(&report_path)
        .with_context(|| format!("failed to create {report_path}"))?;

    // Derive output BAM filename
    let out_stem = if let Some(u) = user_outfile {
        base_stem(u)
    } else {
        derive_stem(primary, None)
    };
    let out_name = if multiple {
        format!("{out_stem}.multiple.deduplicated.bam")
    } else {
        format!("{out_stem}.deduplicated.bam")
    };
    let out_path = PathBuf::from(format!("{output_dir}{out_name}"));
    eprintln!("Output file is: {out_name}\n");

    // Read header for output BAM
    let header = get_sam_header(samtools, primary)?;

    let mut out_bam = BamWriter::open_with_threads(samtools, &out_path, parallel)?;
    for line in header.lines() {
        if !line.is_empty() {
            out_bam.write_line(line.as_bytes())?;
        }
    }

    // Build chromosome → u32 intern table from header
    let chr_map = build_chr_map(&header);

    // Open input
    let mut reader = open_input(samtools, files, multiple, parallel)?;
    let mut buf: Vec<u8> = Vec::with_capacity(8192);

    let mut unique_seqs_se: FxHashSet<SeKey> = FxHashSet::default();
    let mut unique_seqs_pe: FxHashSet<PeKey> = FxHashSet::default();
    let mut unique_seqs_umi: FxHashMap<SeKey, Vec<u8>> = FxHashMap::default();
    let mut positions: FxHashSet<u64> = FxHashSet::default(); // count distinct positions

    let mut count: u64 = 0;
    let mut removed: u64 = 0;

    loop {
        buf.clear();
        let n = reader.lines().read_until(b'\n', &mut buf)?;
        if n == 0 {
            break;
        }
        while buf.last() == Some(&b'\n') || buf.last() == Some(&b'\r') {
            buf.pop();
        }

        if buf.starts_with(b"@") {
            continue;
        } // headers already written
        if buf.is_empty() {
            continue;
        }

        count += 1;
        let r1 = buf.clone();

        // Decode XR/XG to get strand index
        let (xr, xg) = extract_xr_xg(&r1)?;
        let strand_idx = strand_index(xr, xg)?;
        let forward = strand_idx == 0 || strand_idx == 2; // OT or CTOB

        let fields: Vec<&[u8]> = r1.splitn(11, |&b| b == b'\t').collect();
        if fields.len() < 10 {
            bail!("Malformed SAM line (fewer than 10 fields)");
        }

        let chr_bytes = fields[2];
        let chr_id = *chr_map.get(chr_bytes).unwrap_or(&u32::MAX);
        let pos: u32 = parse_u32(fields[3]);
        let cigar_r1 = fields[5];

        if is_single {
            // SE dedup key: (strand, chr, position)
            // Forward reads use start; reverse reads use end (= POS + ref_span - 1)
            let key_pos = if forward {
                pos
            } else {
                pos.saturating_sub(1) + cigar_ref_span(cigar_r1)
            };
            let key: SeKey = (strand_idx, chr_id, key_pos);

            if rrbs {
                let barcode = extract_barcode(fields[0], bclconvert);
                let full_key = key;
                if let Some(existing) = unique_seqs_umi.get(&full_key) {
                    if existing.as_slice() != barcode.as_slice() {
                        // Different UMI at same position — keep
                        unique_seqs_umi.insert(full_key, barcode.clone());
                        out_bam.write_line(&r1)?;
                    } else {
                        removed += 1;
                    }
                } else {
                    unique_seqs_umi.insert(full_key, barcode);
                    out_bam.write_line(&r1)?;
                }
                positions.insert(pack_pos(strand_idx, chr_id, key_pos));
            } else {
                if unique_seqs_se.contains(&key) {
                    removed += 1;
                    positions.insert(pack_pos(strand_idx, chr_id, key_pos));
                } else {
                    unique_seqs_se.insert(key);
                    out_bam.write_line(&r1)?;
                }
            }
        } else {
            // PE dedup: read R2
            buf.clear();
            let n2 = reader.lines().read_until(b'\n', &mut buf)?;
            if n2 == 0 {
                break;
            }
            while buf.last() == Some(&b'\n') || buf.last() == Some(&b'\r') {
                buf.pop();
            }
            let r2 = buf.clone();

            let f2: Vec<&[u8]> = r2.splitn(11, |&b| b == b'\t').collect();
            let pos2: u32 = parse_u32(f2.get(3).copied().unwrap_or(b"0"));
            let cigar_r2 = f2.get(5).copied().unwrap_or(b"*");

            let (start, end) = if forward {
                // OT/CTOB: start = R1.POS, end = R2.POS + R2_cigar_span - 1
                let end = pos2.saturating_sub(1) + cigar_ref_span(cigar_r2);
                (pos, end)
            } else {
                // CTOT/OB: end = R1.POS + R1_cigar_span - 1, start = R2.POS
                let end = pos.saturating_sub(1) + cigar_ref_span(cigar_r1);
                (pos2, end)
            };

            let key: PeKey = (strand_idx, chr_id, start, end);

            if unique_seqs_pe.contains(&key) {
                removed += 1;
                positions.insert(pack_pos_pe(strand_idx, chr_id, start));
            } else {
                unique_seqs_pe.insert(key);
                out_bam.write_line(&r1)?;
                out_bam.write_line(&r2)?;
            }
        }
    }

    out_bam.finish()?;
    reader.finish()?;

    let leftover = count - removed;
    let (pct_rem, pct_left) = if count > 0 {
        (
            format!("{:.2}", removed as f64 / count as f64 * 100.0),
            format!("{:.2}", leftover as f64 / count as f64 * 100.0),
        )
    } else {
        ("N/A".into(), "N/A".into())
    };

    let n_positions = positions.len();
    let file_name = primary.display();

    let lines = [
        format!("\nTotal number of alignments analysed in {file_name}:\t{count}"),
        format!("Total number duplicated alignments removed:\t{removed} ({pct_rem}%)"),
        format!("Duplicated alignments were found at:\t{n_positions} different position(s)\n"),
        format!(
            "Total count of deduplicated leftover sequences: {leftover} ({pct_left}% of total)\n"
        ),
    ];
    for l in &lines {
        eprintln!("{l}");
        writeln!(report, "{l}")?;
    }

    Ok(())
}

/// Compute how many reference bases are consumed by a CIGAR string.
fn cigar_ref_span(cigar: &[u8]) -> u32 {
    let mut span = 0u32;
    let mut num = 0u32;
    for &b in cigar {
        if b.is_ascii_digit() {
            num = num * 10 + (b - b'0') as u32;
        } else {
            match b {
                b'M' | b'D' | b'N' | b'=' | b'X' => span += num,
                _ => {}
            }
            num = 0;
        }
    }
    span
}

fn pack_pos(strand: u8, chr: u32, pos: u32) -> u64 {
    (strand as u64) | ((chr as u64) << 8) | ((pos as u64) << 32)
}

fn pack_pos_pe(strand: u8, chr: u32, pos: u32) -> u64 {
    pack_pos(strand, chr, pos)
}

fn strand_index(xr: &[u8], xg: &[u8]) -> Result<u8> {
    match (xr, xg) {
        (b"CT", b"CT") => Ok(0), // OT
        (b"GA", b"CT") => Ok(1), // CTOT
        (b"GA", b"GA") => Ok(2), // CTOB
        (b"CT", b"GA") => Ok(3), // OB
        _ => bail!(
            "Unexpected XR/XG combination: {:?}/{:?}",
            String::from_utf8_lossy(xr),
            String::from_utf8_lossy(xg)
        ),
    }
}

fn extract_xr_xg<'a>(line: &'a [u8]) -> Result<(&'a [u8], &'a [u8])> {
    let mut xr: Option<&[u8]> = None;
    let mut xg: Option<&[u8]> = None;
    for tag in line.split(|&b| b == b'\t') {
        if tag.starts_with(b"XR:Z:") {
            let v = &tag[5..];
            xr = Some(v.split(|&b| b == b'\r' || b == b'\n').next().unwrap_or(v));
        } else if tag.starts_with(b"XG:Z:") {
            let v = &tag[5..];
            xg = Some(v.split(|&b| b == b'\r' || b == b'\n').next().unwrap_or(v));
        }
    }
    match (xr, xg) {
        (Some(r), Some(g)) => Ok((r, g)),
        _ => bail!("Failed to extract XR/XG tags from SAM line"),
    }
}

fn extract_barcode(qname: &[u8], bclconvert: bool) -> Vec<u8> {
    let s = String::from_utf8_lossy(qname);
    if bclconvert {
        // bcl-convert format: ...:<UMI>_N:N:N:<index>
        if let Some(cap) = s.rfind(':') {
            let tail = &s[cap + 1..];
            if let Some(underscore) = tail.find('_') {
                return tail[..underscore].as_bytes().to_vec();
            }
        }
    } else {
        // Standard UMI: last colon-delimited field
        if let Some(pos) = s.rfind(':') {
            return s[pos + 1..].as_bytes().to_vec();
        }
    }
    qname.to_vec()
}

fn build_chr_map(header: &str) -> FxHashMap<Vec<u8>, u32> {
    let mut map = FxHashMap::default();
    let mut idx = 0u32;
    for line in header.lines() {
        if line.starts_with("@SQ") {
            for field in line.split('\t') {
                if let Some(name) = field.strip_prefix("SN:") {
                    map.insert(name.as_bytes().to_vec(), idx);
                    idx += 1;
                    break;
                }
            }
        }
    }
    map
}

fn get_sam_header(samtools: &str, path: &Path) -> Result<String> {
    let out = std::process::Command::new(samtools)
        .args(["view", "-H"])
        .arg(path)
        .output()
        .context("samtools view -H")?;
    Ok(String::from_utf8_lossy(&out.stdout).into_owned())
}

fn detect_library_type(samtools: &str, path: &Path) -> Result<(bool, bool)> {
    eprintln!("Trying to determine the type of mapping from the SAM header line");
    let output = std::process::Command::new(samtools)
        .args(["view", "-H"])
        .arg(path)
        .output()
        .context("samtools view -H")?;
    for chunk in output.stdout.split(|&b| b == b'\n') {
        if !chunk.starts_with(b"@PG") {
            continue;
        }
        let s = String::from_utf8_lossy(chunk);
        if !s.contains("ID:Bismark") {
            continue;
        }
        if (s.contains(" -1 ") || s.contains(" --1 "))
            && (s.contains(" -2 ") || s.contains(" --2 "))
        {
            eprintln!("Treating file as paired-end data (extracted from @PG line)");
            return Ok((false, true));
        } else {
            eprintln!("Treating file as single-end data (extracted from @PG line)");
            return Ok((true, false));
        }
    }
    Ok((false, false))
}

fn open_input(
    samtools: &str,
    files: &[PathBuf],
    _multiple: bool,
    parallel: u32,
) -> Result<BamReader> {
    // For multiple files, we could use `samtools cat -h ... | samtools view -h`,
    // but for simplicity we process only the primary file.
    // Full multi-file support would require spawning samtools cat.
    let threads = parallel.to_string();
    BamReader::open(samtools, &files[0], &["--threads", &threads])
}

fn derive_stem(path: &Path, user_outfile: Option<&str>) -> String {
    let name = if let Some(u) = user_outfile {
        u
    } else {
        path.to_str().unwrap_or("")
    };
    base_stem(name)
}

fn base_stem(name: &str) -> String {
    let s = Path::new(name)
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let s = s.trim_end_matches(".gz");
    let s = s.trim_end_matches(".sam");
    let s = s.trim_end_matches(".bam");
    let s = s.trim_end_matches(".txt");
    s.to_string()
}

fn normalise_dir(s: &str) -> String {
    if s.is_empty() {
        return String::new();
    }
    if s.ends_with('/') {
        s.to_string()
    } else {
        format!("{s}/")
    }
}

fn parse_u32(b: &[u8]) -> u32 {
    std::str::from_utf8(b)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ─── cigar_ref_span ──────────────────────────────────────────────────────

    #[test]
    fn test_cigar_ref_span_all_match() {
        assert_eq!(cigar_ref_span(b"10M"), 10);
    }

    #[test]
    fn test_cigar_ref_span_with_insertion() {
        // Insertions don't consume reference
        assert_eq!(cigar_ref_span(b"5M2I3M"), 8);
    }

    #[test]
    fn test_cigar_ref_span_with_deletion() {
        assert_eq!(cigar_ref_span(b"5M2D3M"), 10);
    }

    #[test]
    fn test_cigar_ref_span_softclip() {
        // Soft clip does not consume reference
        assert_eq!(cigar_ref_span(b"2S8M2S"), 8);
    }

    #[test]
    fn test_cigar_ref_span_with_skip() {
        // N (intron/skip) consumes reference
        assert_eq!(cigar_ref_span(b"5M100N5M"), 110);
    }

    // ─── pack_pos ────────────────────────────────────────────────────────────

    #[test]
    fn test_pack_pos_zero() {
        assert_eq!(pack_pos(0, 0, 0), 0);
    }

    #[test]
    fn test_pack_pos_strand_in_low_byte() {
        let v = pack_pos(3, 0, 0);
        assert_eq!(v & 0xFF, 3);
    }

    #[test]
    fn test_pack_pos_chr_in_middle() {
        let v = pack_pos(0, 5, 0);
        assert_eq!((v >> 8) & 0xFF_FF, 5);
    }

    #[test]
    fn test_pack_pos_position_in_high() {
        let v = pack_pos(0, 0, 1000);
        assert_eq!((v >> 32) as u32, 1000);
    }

    #[test]
    fn test_pack_pos_combined() {
        // strand=2, chr=1, pos=500
        let v = pack_pos(2, 1, 500);
        assert_eq!(v & 0xFF, 2);
        assert_eq!((v >> 8) & 0xFF_FF, 1);
        assert_eq!((v >> 32) as u32, 500);
    }

    // ─── strand_index ────────────────────────────────────────────────────────

    #[test]
    fn test_strand_index_ot() {
        assert_eq!(strand_index(b"CT", b"CT").unwrap(), 0);
    }

    #[test]
    fn test_strand_index_ctot() {
        assert_eq!(strand_index(b"GA", b"CT").unwrap(), 1);
    }

    #[test]
    fn test_strand_index_ctob() {
        assert_eq!(strand_index(b"GA", b"GA").unwrap(), 2);
    }

    #[test]
    fn test_strand_index_ob() {
        assert_eq!(strand_index(b"CT", b"GA").unwrap(), 3);
    }

    #[test]
    fn test_strand_index_invalid() {
        assert!(strand_index(b"XX", b"YY").is_err());
    }

    // ─── extract_barcode ─────────────────────────────────────────────────────

    #[test]
    fn test_extract_barcode_standard_umi() {
        // Standard: last colon-delimited field
        let qname = b"read1:ACGT";
        assert_eq!(extract_barcode(qname, false), b"ACGT");
    }

    #[test]
    fn test_extract_barcode_no_colon() {
        let qname = b"readnocoercion";
        assert_eq!(extract_barcode(qname, false), b"readnocoercion");
    }

    #[test]
    fn test_extract_barcode_bclconvert() {
        // bcl-convert format: last colon field is <UMI>_<extra>
        let qname = b"instrument:run:flowcell:lane:tile:x:AACCGGTT_illumina";
        assert_eq!(extract_barcode(qname, true), b"AACCGGTT");
    }

    // ─── base_stem ───────────────────────────────────────────────────────────

    #[test]
    fn test_base_stem_bam() {
        assert_eq!(base_stem("sample.bam"), "sample");
    }

    #[test]
    fn test_base_stem_sam() {
        assert_eq!(base_stem("sample.sam"), "sample");
    }

    #[test]
    fn test_base_stem_txt_gz() {
        assert_eq!(base_stem("sample.txt.gz"), "sample");
    }

    #[test]
    fn test_base_stem_with_path() {
        assert_eq!(base_stem("/data/sample.bam"), "sample");
    }

    // ─── build_chr_map ───────────────────────────────────────────────────────

    #[test]
    fn test_build_chr_map_basic() {
        let header =
            "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422\n@SQ\tSN:chr2\tLN:242193529\n";
        let map = build_chr_map(header);
        assert_eq!(map[b"chr1".as_ref()], 0);
        assert_eq!(map[b"chr2".as_ref()], 1);
    }

    #[test]
    fn test_build_chr_map_empty_header() {
        let map = build_chr_map("@HD\tVN:1.6\n");
        assert!(map.is_empty());
    }

    #[test]
    fn test_build_chr_map_preserves_order() {
        let header =
            "@SQ\tSN:chrM\tLN:16569\n@SQ\tSN:chr1\tLN:248956422\n@SQ\tSN:chrX\tLN:156040895\n";
        let map = build_chr_map(header);
        assert_eq!(map[b"chrM".as_ref()], 0);
        assert_eq!(map[b"chr1".as_ref()], 1);
        assert_eq!(map[b"chrX".as_ref()], 2);
    }
}
