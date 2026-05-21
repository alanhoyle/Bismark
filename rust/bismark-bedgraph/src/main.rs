use std::collections::{BTreeMap, BinaryHeap};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use bismark_lib::BISMARK_VERSION;
use clap::Parser;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

#[derive(Parser)]
#[command(
    name = "bismark2bedGraph",
    about = "Convert Bismark methylation extractor output to bedGraph/coverage format",
    version = BISMARK_VERSION,
    disable_version_flag = true,
)]
struct Cli {
    /// Input methylation call file(s) from bismark_methylation_extractor
    #[arg(required = true)]
    files: Vec<PathBuf>,

    /// Output bedGraph filename (required)
    #[arg(short = 'o', long = "output")]
    output: String,

    /// Output directory
    #[arg(long = "dir", default_value = "")]
    output_dir: String,

    /// Input files have no header line
    #[arg(long = "no_header")]
    no_header: bool,

    /// Minimum read coverage to call methylation [default: 1]
    #[arg(long = "cutoff", default_value_t = 1)]
    cutoff: u32,

    /// Replace whitespace in read IDs with underscore
    #[arg(long = "remove_spaces")]
    remove_spaces: bool,

    /// Include all cytosine contexts (CX), not just CpG
    #[arg(long = "CX", alias = "CX_context")]
    cx_context: bool,

    /// Enable disk-based external merge sort; value sets the in-memory run size
    /// (e.g. 2G, 500M).  Omit to sort entirely in memory (default).
    #[arg(long = "buffer_size")]
    buffer_size: Option<String>,

    /// Input has many scaffolds (accepted for compatibility; external sort handles this natively)
    #[arg(long = "gazillion", alias = "scaffolds")]
    gazillion: bool,

    /// Sort in memory — the default; accepted for Perl compatibility
    #[arg(long = "ample_memory")]
    ample_memory: bool,

    /// Also write 0-based half-open coverage file
    #[arg(long = "zero_based")]
    zero_based: bool,

    /// Write additional UCSC-compatible bedGraph (prefix chr, MT→chrM)
    #[arg(long = "ucsc")]
    ucsc: bool,

    /// Print version and exit
    #[arg(long = "version")]
    version: bool,
}

// ─── External merge sort helpers ─────────────────────────────────────────────

/// Parse a size string like "2G", "500M", "1024K", or a plain byte count.
fn parse_buffer_size(s: &str) -> u64 {
    let s = s.trim();
    let (num, mult) = if let Some(n) = s.strip_suffix(['G', 'g']) {
        (n, 1u64 << 30)
    } else if let Some(n) = s.strip_suffix(['M', 'm']) {
        (n, 1u64 << 20)
    } else if let Some(n) = s.strip_suffix(['K', 'k']) {
        (n, 1u64 << 10)
    } else {
        (s, 1u64)
    };
    num.parse::<u64>().unwrap_or(2 << 30) * mult
}

/// Sort `buf` by (chr, pos) and flush it as a tab-delimited temp file.
/// The buffer is cleared on return.
fn flush_sorted_run(
    buf: &mut Vec<(String, u32, bool)>,
    temp_files: &mut Vec<tempfile::NamedTempFile>,
) -> Result<()> {
    buf.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut tf = tempfile::NamedTempFile::new().context("creating temp file")?;
    for (chr, pos, is_meth) in buf.drain(..) {
        writeln!(tf, "{chr}\t{pos}\t{}", is_meth as u8)?;
    }
    tf.flush()?;
    temp_files.push(tf);
    Ok(())
}

/// Parse one line from a sorted-run temp file: "chr\tpos\t0_or_1".
fn parse_temp_line(line: &str) -> Option<(String, u32, bool)> {
    let mut f = line.splitn(3, '\t');
    let chr = f.next()?.to_string();
    let pos: u32 = f.next()?.parse().ok()?;
    let is_meth = f.next()?.trim() == "1";
    Some((chr, pos, is_meth))
}

/// Write one aggregated cytosine record to all active output writers.
fn emit_record(
    chr: &str,
    pos: u32,
    meth: u32,
    unmeth: u32,
    cutoff: u32,
    bg_out: &mut dyn Write,
    cov_out: &mut dyn Write,
    zero_out: &mut Option<std::fs::File>,
) -> Result<()> {
    let total = meth + unmeth;
    if total < cutoff {
        return Ok(());
    }
    let pct = format_percentage(meth, total);
    let bed_pos = pos - 1;
    writeln!(bg_out, "{chr}\t{bed_pos}\t{pos}\t{pct}")?;
    writeln!(cov_out, "{chr}\t{pos}\t{pos}\t{pct}\t{meth}\t{unmeth}")?;
    if let Some(ref mut zf) = zero_out {
        writeln!(zf, "{chr}\t{bed_pos}\t{pos}\t{pct}\t{meth}\t{unmeth}")?;
    }
    Ok(())
}

/// K-way merge of sorted temp files, aggregating same-position records,
/// writing directly to the output writers.
fn merge_and_write(
    temp_files: Vec<tempfile::NamedTempFile>,
    cutoff: u32,
    bg_out: &mut dyn Write,
    cov_out: &mut dyn Write,
    zero_out: &mut Option<std::fs::File>,
) -> Result<()> {
    // Open a fresh read handle for each temp file.
    let mut readers: Vec<std::io::Lines<BufReader<std::fs::File>>> = temp_files
        .iter()
        .map(|tf| {
            let f = std::fs::File::open(tf.path())
                .with_context(|| format!("reopening temp file {}", tf.path().display()))?;
            Ok(BufReader::new(f).lines())
        })
        .collect::<Result<_>>()?;

    // Min-heap entries: Reverse so BinaryHeap (max) behaves as a min-heap.
    // Tuple: (chr, pos, is_meth_u8, file_idx)
    let mut heap: BinaryHeap<std::cmp::Reverse<(String, u32, u8, usize)>> = BinaryHeap::new();

    // Seed one record from each file.
    for (i, reader) in readers.iter_mut().enumerate() {
        if let Some(Ok(line)) = reader.next() {
            if let Some((chr, pos, is_meth)) = parse_temp_line(&line) {
                heap.push(std::cmp::Reverse((chr, pos, is_meth as u8, i)));
            }
        }
    }

    let mut cur_chr = String::new();
    let mut cur_pos: u32 = 0;
    let mut meth: u32 = 0;
    let mut unmeth: u32 = 0;
    let mut started = false;

    while let Some(std::cmp::Reverse((chr, pos, is_meth_u8, file_idx))) = heap.pop() {
        // New position: emit the accumulated record for the previous position.
        if started && (chr != cur_chr || pos != cur_pos) {
            emit_record(&cur_chr, cur_pos, meth, unmeth, cutoff, bg_out, cov_out, zero_out)?;
            meth = 0;
            unmeth = 0;
        }
        cur_chr = chr;
        cur_pos = pos;
        started = true;
        if is_meth_u8 == 1 { meth += 1; } else { unmeth += 1; }

        // Advance this file and push its next record.
        if let Some(Ok(line)) = readers[file_idx].next() {
            if let Some((chr2, pos2, is_meth2)) = parse_temp_line(&line) {
                heap.push(std::cmp::Reverse((chr2, pos2, is_meth2 as u8, file_idx)));
            }
        }
    }

    // Emit the final accumulated record.
    if started {
        emit_record(&cur_chr, cur_pos, meth, unmeth, cutoff, bg_out, cov_out, zero_out)?;
    }

    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.version {
        println!(
            "\n\n                      Bismark Methylation Extractor Module -\n                                bismark2bedGraph\n\n                      Bismark Extractor Version: {}\n              Copyright 2010-22 Felix Krueger, Altos Bioinformatics\n                     https://github.com/FelixKrueger/Bismark\n\n",
            BISMARK_VERSION
        );
        return Ok(());
    }

    if cli.output.contains('/') {
        bail!("Please specify a file name without any path information (or use --dir if necessary)");
    }

    let output_dir = normalise_dir(&cli.output_dir);

    // Ensure output ends in .gz
    let bedgraph_name = if cli.output.ends_with(".gz") {
        cli.output.clone()
    } else {
        format!("{}.gz", cli.output)
    };

    // Determine which files to use based on CX flag
    let input_files: Vec<PathBuf> = if cli.cx_context {
        cli.files.clone()
    } else {
        cli.files
            .iter()
            .filter(|f| {
                f.file_name()
                    .and_then(|n| n.to_str())
                    .map(|n| n.starts_with("CpG"))
                    .unwrap_or(false)
            })
            .cloned()
            .collect()
    };

    if input_files.is_empty() {
        bail!("It seems that you are trying to generate bedGraph files for files not starting with CpG.... Please specify the option '--CX' and try again");
    }

    eprintln!("Using the following files as Input:");
    for f in &input_files {
        eprint!("\t{}", f.display());
    }
    eprintln!("\n");

    eprintln!("\nSummary of parameters for bismark2bedGraph conversion:");
    eprintln!("{}", "=".repeat(54));
    eprintln!("bedGraph output:\t\t{bedgraph_name}");
    eprintln!("output directory:\t\t>{output_dir}<");
    eprintln!("remove whitespaces:\t\t{}", if cli.remove_spaces { "yes" } else { "no" });
    eprintln!("CX context:\t\t\t{}", if cli.cx_context { "yes" } else { "no (CpG context only, default)" });
    eprintln!("No-header selected:\t\t{}", if cli.no_header { "yes" } else { "no" });
    if let Some(ref bs) = cli.buffer_size {
        eprintln!("Sorting method:\t\t\texternal merge sort (disk-based)");
        eprintln!("Sort buffer size:\t\t{bs}");
    } else {
        eprintln!("Sorting method:\t\t\tin-memory sort (default)");
    }
    eprintln!("Coverage threshold:\t\t{}", cli.cutoff);
    eprintln!("{}", "=".repeat(77));
    eprintln!("Methylation information will now be written into a bedGraph and coverage file");
    eprintln!("{}\n", "=".repeat(77));

    // Coverage output filename
    let coverage_name = if bedgraph_name.ends_with("bedGraph.gz") {
        bedgraph_name.replace("bedGraph.gz", "bismark.cov.gz")
    } else {
        format!("{}.bismark.cov.gz", bedgraph_name)
    };

    // Zero-based coverage output filename
    let zero_name = if cli.zero_based {
        let n = bedgraph_name.trim_end_matches(".gz");
        let n = if n.ends_with("bedGraph") {
            format!("{}.bismark.zero.cov", &n[..n.len() - "bedGraph".len()])
        } else {
            format!("{}.bismark.zero.cov", n)
        };
        Some(n)
    } else {
        None
    };

    eprintln!("Writing bedGraph to file: {bedgraph_name}");
    eprintln!("Also writing out a coverage file including counts methylated and unmethylated residues to file: {coverage_name}");
    if let Some(ref zn) = zero_name {
        eprintln!("Also writing out a 0-based, half-open coverage file including counts methylated and unmethylated residues to file: {zn}");
    }
    eprintln!();

    if cli.ucsc {
        eprintln!("Creating additional bedGraph file that has known Ensembl chromosome names replaced to work with the UCSC genome browser");
        eprintln!("This option:\n- prefixes chromosome names with 'chr'");
        eprintln!("- changes 'MT' to 'chrM'\n");
    }

    let bedgraph_path = format!("{}{}", output_dir, bedgraph_name);
    let coverage_path = format!("{}{}", output_dir, coverage_name);

    let mut bg_out = GzEncoder::new(
        std::fs::File::create(&bedgraph_path)
            .with_context(|| format!("creating {bedgraph_path}"))?,
        Compression::default(),
    );
    let mut cov_out = GzEncoder::new(
        std::fs::File::create(&coverage_path)
            .with_context(|| format!("creating {coverage_path}"))?,
        Compression::default(),
    );
    let mut zero_out: Option<std::fs::File> = if let Some(ref zn) = zero_name {
        let p = format!("{}{}", output_dir, zn);
        Some(std::fs::File::create(&p).with_context(|| format!("creating {p}"))?)
    } else {
        None
    };

    // Write bedGraph track header
    writeln!(bg_out, "track type=bedGraph")?;

    if let Some(ref buf_size_str) = cli.buffer_size {
        // ── External merge sort path ─────────────────────────────────────────
        let buffer_bytes = parse_buffer_size(buf_size_str);
        eprintln!("External sort enabled; in-memory run size: {buf_size_str} ({buffer_bytes} bytes)");

        let mut buf: Vec<(String, u32, bool)> = Vec::new();
        let mut buf_est: u64 = 0;
        let mut temp_files: Vec<tempfile::NamedTempFile> = Vec::new();

        for infile in &input_files {
            // Re-use read_methylation_file but accumulate into a flat buffer
            // by temporarily wrapping with a BTreeMap shim.
            let mut shard: BTreeMap<String, Vec<(u32, bool)>> = BTreeMap::new();
            read_methylation_file(infile, &cli, &mut shard)
                .with_context(|| format!("reading {}", infile.display()))?;
            for (chr, positions) in shard {
                for (pos, is_meth) in positions {
                    buf_est += chr.len() as u64 + 16;
                    buf.push((chr.clone(), pos, is_meth));
                }
                if buf_est >= buffer_bytes {
                    flush_sorted_run(&mut buf, &mut temp_files)?;
                    buf_est = 0;
                }
            }
        }
        if !buf.is_empty() {
            flush_sorted_run(&mut buf, &mut temp_files)?;
        }

        eprintln!("Merging {} sorted run(s)…", temp_files.len());
        merge_and_write(temp_files, cli.cutoff, &mut bg_out, &mut cov_out, &mut zero_out)?;
    } else {
        // ── In-memory sort path (default) ────────────────────────────────────
        let mut data: BTreeMap<String, Vec<(u32, bool)>> = BTreeMap::new();
        for infile in &input_files {
            read_methylation_file(infile, &cli, &mut data)
                .with_context(|| format!("reading {}", infile.display()))?;
        }

        for (chr, mut positions) in data {
            positions.sort_unstable_by_key(|&(pos, _)| pos);

            let mut i = 0;
            while i < positions.len() {
                let pos = positions[i].0;
                let mut meth: u32 = 0;
                let mut unmeth: u32 = 0;
                while i < positions.len() && positions[i].0 == pos {
                    if positions[i].1 { meth += 1; } else { unmeth += 1; }
                    i += 1;
                }
                emit_record(&chr, pos, meth, unmeth, cli.cutoff, &mut bg_out, &mut cov_out, &mut zero_out)?;
            }
        }
    }

    bg_out.finish().context("finalizing bedGraph gzip")?;
    cov_out.finish().context("finalizing coverage gzip")?;

    eprintln!("Finished writing bedGraph and coverage files.");

    if cli.ucsc {
        write_ucsc_bedgraph(&bedgraph_path, &output_dir, &bedgraph_name)?;
    }

    Ok(())
}

fn format_percentage(meth: u32, total: u32) -> String {
    let pct = meth as f64 / total as f64 * 100.0;
    if pct == 0.0 {
        return "0".to_string();
    }

    // Perl's default numeric stringification uses roughly 15 significant
    // digits. Match that so differential output does not diverge on values
    // like 1/3 and 2/3.
    let digits_before_decimal = if pct >= 1.0 {
        pct.log10().floor() as i32 + 1
    } else {
        0
    };
    let decimals = (15 - digits_before_decimal).max(0) as usize;
    let mut s = format!("{pct:.decimals$}");
    if s.contains('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
    }
    s
}

fn read_methylation_file(
    path: &Path,
    cli: &Cli,
    data: &mut BTreeMap<String, Vec<(u32, bool)>>,
) -> Result<()> {
    let reader: Box<dyn BufRead> = if path.to_string_lossy().ends_with(".gz") {
        let f = std::fs::File::open(path)?;
        Box::new(BufReader::new(MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::new(std::fs::File::open(path)?))
    };

    let fname = path.file_name().unwrap().to_string_lossy();
    eprintln!("Now writing methylation information for file >>{fname}<<");

    let mut lines = reader.lines();

    // Skip header line unless --no_header
    if !cli.no_header {
        if let Some(first) = lines.next() {
            let first = first?;
            // If it doesn't look like a Bismark header, treat it as data
            if !first.starts_with("Bismark") && !first.starts_with("ReadID") {
                process_meth_line(&first, cli, data);
            }
        }
    }

    for line in lines {
        let line = line?;
        if line.starts_with("Bismark") || line.is_empty() {
            continue;
        }
        process_meth_line(&line, cli, data);
    }

    eprintln!("Finished writing out individual chromosome files for {fname}");
    Ok(())
}

fn process_meth_line(line: &str, cli: &Cli, data: &mut BTreeMap<String, Vec<(u32, bool)>>) {
    let fields: Vec<&str> = line.splitn(5, '\t').collect();
    if fields.len() < 5 {
        return;
    }

    let mut _name = fields[0];
    let mut _name_owned;
    if cli.remove_spaces && _name.contains(' ') {
        _name_owned = _name.replace(' ', "_");
        _name = &_name_owned;
    }

    let meth_state = fields[1]; // '+' or '-'
    let chr = fields[2];
    let pos: u32 = match fields[3].parse() {
        Ok(p) => p,
        Err(_) => return,
    };
    let ctx = fields[4].trim();

    if !validate_methylation_call(meth_state, ctx) {
        return;
    }

    let is_meth = meth_state == "+";
    data.entry(chr.to_string()).or_default().push((pos, is_meth));
}

fn validate_methylation_call(meth_state: &str, ctx: &str) -> bool {
    let ctx_char = ctx.chars().next().unwrap_or('?');
    match meth_state {
        "+" => matches!(ctx_char, 'Z' | 'X' | 'H'),
        "-" => matches!(ctx_char, 'z' | 'x' | 'h'),
        _ => false,
    }
}

fn write_ucsc_bedgraph(bedgraph_gz: &str, output_dir: &str, bedgraph_name: &str) -> Result<()> {
    eprintln!("Finally, creating a UCSC compatible version of the bedGraph file. This shouldn't take long...");
    eprintln!("     ==============     ==============     ==============     ==============     ==============");

    let ucsc_stem = bedgraph_name.trim_end_matches(".gz");
    let ucsc_name = format!("{}_UCSC.bedGraph.gz", ucsc_stem);
    let ucsc_path = format!("{}{}", output_dir, ucsc_name);

    eprintln!("Writing a new version of the bedGraph file compatible to UCSC  genomes to file >{ucsc_name}<");
    eprintln!("Chromosomes names will start with 'chr'");
    eprintln!("Chromosome MT (Ensembl mitochondrial DNA http://www.ensembl.org/Homo_sapiens/Location/Chromosome?chr=MT;r=MT:1-16569) will be renamed to 'chrM'");

    let reader: Box<dyn BufRead> = if bedgraph_gz.ends_with(".gz") {
        let f = std::fs::File::open(bedgraph_gz)
            .with_context(|| format!("opening {bedgraph_gz}"))?;
        Box::new(BufReader::new(MultiGzDecoder::new(f)))
    } else {
        let f = std::fs::File::open(bedgraph_gz)
            .with_context(|| format!("opening {bedgraph_gz}"))?;
        Box::new(BufReader::new(f))
    };

    let out = std::fs::File::create(&ucsc_path)
        .with_context(|| format!("creating {ucsc_path}"))?;
    let mut gz_out = GzEncoder::new(out, Compression::default());

    let mut lines = reader.lines();
    // Write first line (track header) unchanged
    if let Some(first) = lines.next() {
        writeln!(gz_out, "{}", first?)?;
    }

    for line in lines {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let mut parts = line.splitn(4, '\t');
        let chr = parts.next().unwrap_or("");
        let rest = &line[chr.len()..];

        let ucsc_chr = if chr == "MT" {
            "chrM".to_string()
        } else if chr.starts_with("chr") {
            chr.to_string()
        } else {
            format!("chr{chr}")
        };

        writeln!(gz_out, "{ucsc_chr}{rest}")?;
    }

    gz_out.finish().context("finalizing UCSC bedGraph gzip")?;
    eprintln!("\nAll done. Finished UCSC conversion.\n");
    Ok(())
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

#[cfg(test)]
mod tests {
    use super::*;

    // ─── validate_methylation_call ───────────────────────────────────────────

    #[test]
    fn test_validate_meth_plus_valid() {
        assert!(validate_methylation_call("+", "Z"));
        assert!(validate_methylation_call("+", "X"));
        assert!(validate_methylation_call("+", "H"));
    }

    #[test]
    fn test_validate_unmeth_minus_valid() {
        assert!(validate_methylation_call("-", "z"));
        assert!(validate_methylation_call("-", "x"));
        assert!(validate_methylation_call("-", "h"));
    }

    #[test]
    fn test_validate_wrong_case() {
        // lowercase z with + is invalid
        assert!(!validate_methylation_call("+", "z"));
        // uppercase Z with - is invalid
        assert!(!validate_methylation_call("-", "Z"));
    }

    #[test]
    fn test_validate_unknown_state() {
        assert!(!validate_methylation_call("?", "Z"));
        assert!(!validate_methylation_call("", "Z"));
    }

    #[test]
    fn test_validate_with_longer_context_string() {
        // context string may have extra chars; only first char matters
        assert!(validate_methylation_call("+", "ZCG"));
        assert!(validate_methylation_call("-", "zCG"));
    }

    // ─── normalise_dir ───────────────────────────────────────────────────────

    #[test]
    fn test_normalise_dir_empty() {
        assert_eq!(normalise_dir(""), "");
    }

    #[test]
    fn test_normalise_dir_adds_slash() {
        assert_eq!(normalise_dir("output"), "output/");
    }

    #[test]
    fn test_normalise_dir_keeps_slash() {
        assert_eq!(normalise_dir("output/"), "output/");
    }
}
