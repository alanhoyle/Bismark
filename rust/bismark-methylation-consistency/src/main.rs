use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use bismark_lib::bam_io::{bam_is_truncated, find_samtools, BamReader, BamWriter};
use bismark_lib::BISMARK_VERSION;
use clap::Parser;

#[derive(Parser)]
#[command(
    name = "methylation_consistency",
    about = "Split Bismark BAM files by CpG methylation consistency",
    version = BISMARK_VERSION,
    disable_version_flag = true,
)]
struct Cli {
    /// BAM file(s) to process
    #[arg(required = true)]
    files: Vec<PathBuf>,

    /// Minimum number of CpGs (or CHHs) required per read [default: 5]
    #[arg(long = "min-count", short = 'm', default_value_t = 5)]
    min_count: u32,

    /// Percentage below which a read is considered fully unmethylated [default: 10]
    #[arg(long, default_value_t = 10)]
    lower_threshold: u32,

    /// Percentage above which a read is considered fully methylated [default: 90]
    #[arg(long, default_value_t = 90)]
    upper_threshold: u32,

    /// Treat files as single-end (auto-detected if not set)
    #[arg(short = 's', long = "single_end")]
    single: bool,

    /// Treat files as paired-end (auto-detected if not set)
    #[arg(short = 'p', long = "paired_end")]
    paired: bool,

    /// Use CHH context instead of CpG (experimental)
    #[arg(long)]
    chh: bool,

    /// Path to samtools
    #[arg(long)]
    samtools_path: Option<String>,

    /// Print version and exit
    #[arg(long = "version")]
    version: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.version {
        println!(
            "\n             Bismark Read Methylation Constistency\n\t         v{}\n       Copyright 2019-22 Felix Krueger, Altos Bioinformatics\n              https://github.com/FelixKrueger/Bismark\n",
            BISMARK_VERSION
        );
        return Ok(());
    }

    if cli.single && cli.paired {
        bail!("You cannot select both single-end (SE) as well as paired-end (PE). Please settle for just one type...");
    }
    if cli.upper_threshold < 51 || cli.upper_threshold > 100 {
        bail!("The upper methylation threshold needs to be a number between 51 and 100% [default is 90%].");
    }
    if cli.lower_threshold > 49 {
        bail!("The lower methylation threshold needs to be a number between 0 and 49% [default is 10%].");
    }

    if cli.chh {
        eprintln!("     ~~~~~     \nTHIS IS AN EXPERIMENTAL VERSION that works on CHH context and **NOT** on the usual CpG context. You have been warned!\n     ~~~~\n");
    }

    eprintln!("Upper and lower methylation thresholds given as:\nUpper: {}\nLower: {}\n", cli.upper_threshold, cli.lower_threshold);

    let samtools = find_samtools(cli.samtools_path.as_deref())?;

    for file in &cli.files {
        process_file(file, &cli, &samtools)
            .with_context(|| format!("processing {}", file.display()))?;
    }
    Ok(())
}

fn process_file(path: &Path, cli: &Cli, samtools: &str) -> Result<()> {
    eprintln!("Now processing file: {}\n     ~~~~~", path.display());

    if path.to_string_lossy().ends_with(".bam") && bam_is_truncated(samtools, path) {
        bail!("File {} appears truncated", path.display());
    }

    let (_is_single, is_paired) = if cli.single {
        eprintln!("Single-end (SE) mode selected manually");
        (true, false)
    } else if cli.paired {
        eprintln!("Paired-end (PE) mode selected manually - methylation information from both reads are simply added together");
        (false, true)
    } else {
        detect_library_type(samtools, path)?
    };

    // Fetch SAM header to write into all three output BAMs
    let header = get_sam_header(samtools, path)?;

    let stem = path.to_str().unwrap().trim_end_matches(".bam");
    let chh_tag = if cli.chh { "_CHH" } else { "" };

    let meth_path   = PathBuf::from(format!("{stem}{chh_tag}_all_meth.bam"));
    let unmeth_path = PathBuf::from(format!("{stem}{chh_tag}_all_unmeth.bam"));
    let mixed_path  = PathBuf::from(format!("{stem}{chh_tag}_mixed_meth.bam"));
    let report_path = PathBuf::from(format!("{stem}{chh_tag}_consistency_report.txt"));

    let mut out_meth   = BamWriter::open(samtools, &meth_path)?;
    let mut out_unmeth = BamWriter::open(samtools, &unmeth_path)?;
    let mut out_mixed  = BamWriter::open(samtools, &mixed_path)?;
    let mut report     = std::fs::File::create(&report_path)?;

    // Write SAM headers to all output BAMs upfront
    for line in header.lines() {
        if line.is_empty() { continue; }
        out_meth.write_line(line.as_bytes())?;
        out_unmeth.write_line(line.as_bytes())?;
        out_mixed.write_line(line.as_bytes())?;
    }

    let mut reader = BamReader::open(samtools, path, &[])?;
    let mut buf: Vec<u8> = Vec::with_capacity(4096);

    let mut all_meth_count   = 0u64;
    let mut all_unmeth_count = 0u64;
    let mut mixed_meth_count = 0u64;
    let mut discarded_count  = 0u64;

    loop {
        buf.clear();
        let n = reader.lines().read_until(b'\n', &mut buf)?;
        if n == 0 { break; }
        while buf.last() == Some(&b'\n') || buf.last() == Some(&b'\r') { buf.pop(); }
        if buf.starts_with(b"@") || buf.is_empty() { continue; }

        let r1 = buf.clone();
        let (mut meth, mut unmeth) = count_methylation(&r1, cli.chh);

        if is_paired {
            buf.clear();
            let n2 = reader.lines().read_until(b'\n', &mut buf)?;
            if n2 == 0 { break; }
            while buf.last() == Some(&b'\n') || buf.last() == Some(&b'\r') { buf.pop(); }
            let r2 = buf.clone();

            let id1 = first_field(&r1);
            let id2 = first_field(&r2);
            if id1 != id2 {
                bail!(
                    "READ IDs of R1 and R2 did not match ({} vs {}). This doesn't look like paired-end data.",
                    String::from_utf8_lossy(id1), String::from_utf8_lossy(id2)
                );
            }
            let (m2, u2) = count_methylation(&r2, cli.chh);
            meth += m2;
            unmeth += u2;

            if meth + unmeth < cli.min_count {
                discarded_count += 1;
                continue;
            }
            let pct = format!("{:.1}", meth as f64 / (meth + unmeth) as f64 * 100.0)
                .parse::<f64>().unwrap_or(0.0);

            if pct <= cli.lower_threshold as f64 {
                all_unmeth_count += 1;
                out_unmeth.write_line(&r1)?;
                out_unmeth.write_line(&r2)?;
            } else if pct >= cli.upper_threshold as f64 {
                all_meth_count += 1;
                out_meth.write_line(&r1)?;
                out_meth.write_line(&r2)?;
            } else {
                mixed_meth_count += 1;
                out_mixed.write_line(&r1)?;
                out_mixed.write_line(&r2)?;
            }
        } else {
            if meth + unmeth < cli.min_count {
                discarded_count += 1;
                continue;
            }
            let pct = format!("{:.1}", meth as f64 / (meth + unmeth) as f64 * 100.0)
                .parse::<f64>().unwrap_or(0.0);

            if pct <= cli.lower_threshold as f64 {
                all_unmeth_count += 1;
                out_unmeth.write_line(&r1)?;
            } else if pct >= cli.upper_threshold as f64 {
                all_meth_count += 1;
                out_meth.write_line(&r1)?;
            } else {
                mixed_meth_count += 1;
                out_mixed.write_line(&r1)?;
            }
        }
    }

    out_meth.finish()?;
    out_unmeth.finish()?;
    out_mixed.finish()?;
    reader.finish()?;

    let total = all_meth_count + all_unmeth_count + mixed_meth_count + discarded_count;
    let (pm, pu, pmix, pd) = if total > 0 {
        (
            format!("{:.2}", all_meth_count   as f64 / total as f64 * 100.0),
            format!("{:.2}", all_unmeth_count  as f64 / total as f64 * 100.0),
            format!("{:.2}", mixed_meth_count  as f64 / total as f64 * 100.0),
            format!("{:.2}", discarded_count   as f64 / total as f64 * 100.0),
        )
    } else {
        ("N/A".into(), "N/A".into(), "N/A".into(), "N/A".into())
    };

    let type_str = if is_paired { "paired-end" } else { "single-end" };
    let lo = cli.lower_threshold;
    let hi = cli.upper_threshold;
    let min = cli.min_count;
    let ctx = if cli.chh { "CHH" } else { "CpG" };

    let lines = [
        format!("Total {type_str} records     -\t{total}"),
        "-------------------------------------------------".to_string(),
        format!("All methylated    [ >= {hi}% ] -\t{all_meth_count} ({pm}%)"),
        format!("All unmethylated  [ <= {lo}% ] -\t{all_unmeth_count} ({pu}%)"),
        format!("Mixed methylation [ {lo}-{hi}% ] -\t{mixed_meth_count} ({pmix}%)"),
        format!("Too few {ctx}s   [min-count {min}] -\t{discarded_count} ({pd}%)"),
    ];

    eprintln!("Summary for {}:\n", path.display());
    for l in &lines {
        eprintln!("{l}");
        writeln!(report, "{l}")?;
    }

    Ok(())
}

fn count_methylation(line: &[u8], use_chh: bool) -> (u32, u32) {
    let xm = match extract_xm(line) { Some(x) => x, None => return (0, 0) };
    let (mc, uc) = if use_chh { (b'H', b'h') } else { (b'Z', b'z') };
    let (mut m, mut u) = (0u32, 0u32);
    for &b in xm { if b == mc { m += 1; } else if b == uc { u += 1; } }
    (m, u)
}

fn extract_xm(line: &[u8]) -> Option<&[u8]> {
    const M: &[u8] = b"XM:Z:";
    let p = line.windows(M.len()).position(|w| w == M)?;
    let s = p + M.len();
    let e = line[s..].iter().position(|&b| b == b'\t').map(|x| s + x).unwrap_or(line.len());
    Some(&line[s..e])
}

fn first_field(line: &[u8]) -> &[u8] {
    line.iter().position(|&b| b == b'\t').map(|p| &line[..p]).unwrap_or(line)
}

fn get_sam_header(samtools: &str, path: &Path) -> Result<String> {
    let out = std::process::Command::new(samtools).args(["view", "-H"]).arg(path)
        .output().context("samtools view -H failed")?;
    Ok(String::from_utf8_lossy(&out.stdout).into_owned())
}

fn detect_library_type(samtools: &str, path: &Path) -> Result<(bool, bool)> {
    eprintln!("Trying to determine the type of mapping from the SAM header line");
    let output = std::process::Command::new(samtools)
        .args(["view", "-H"]).arg(path).output().context("samtools view -H")?;
    for chunk in output.stdout.split(|&b| b == b'\n') {
        if !chunk.starts_with(b"@PG") { continue; }
        let s = String::from_utf8_lossy(chunk);
        if !s.contains("ID:Bismark") { continue; }
        if (s.contains(" -1 ") || s.contains(" --1 ")) && (s.contains(" -2 ") || s.contains(" --2 ")) {
            eprintln!("Paired-end (PE) mode selected (auto-detected) - methylation information from both reads are simply added together");
            return Ok((false, true));
        } else {
            eprintln!("Single-end (SE) mode selected (auto-detected)");
            return Ok((true, false));
        }
    }
    Ok((true, false))
}
