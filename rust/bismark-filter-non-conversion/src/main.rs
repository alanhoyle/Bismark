use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use bismark_lib::bam_io::{bam_is_truncated, find_samtools, BamReader, BamWriter};
use bismark_lib::BISMARK_VERSION;
use clap::Parser;

#[derive(Parser)]
#[command(
    name = "filter_non_conversion",
    about = "Filtering incomplete bisulfite conversion from Bismark BAM files",
    version = BISMARK_VERSION,
    disable_version_flag = true,
)]
struct Cli {
    /// BAM file(s) to process
    #[arg(required = true)]
    files: Vec<PathBuf>,

    /// Filter single-end files (auto-detected if not set)
    #[arg(short = 's', long = "single")]
    single: bool,

    /// Filter paired-end files (auto-detected if not set)
    #[arg(short = 'p', long = "paired")]
    paired: bool,

    /// Number of methylated non-CG calls at which a read is removed [default: 3]
    #[arg(long, default_value_t = 3)]
    threshold: u32,

    /// Filter by percentage of non-CG methylation instead of absolute count
    #[arg(long)]
    percentage_cutoff: Option<u32>,

    /// Minimum non-CG cytosine count required before --percentage_cutoff applies [default: 5]
    #[arg(long, default_value_t = 5)]
    minimum_count: u32,

    /// Non-CG methylation must be consecutive; any unmethylated C resets the counter
    #[arg(long)]
    consecutive: bool,

    /// Path to samtools installation
    #[arg(long)]
    samtools_path: Option<String>,

    /// Print version and exit
    #[arg(long = "version", short = 'V')]
    version: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.version {
        println!(
            "\n\t\tBismark non-conversion filtering\n\t    \n\t\t   Bismark non-conversion version: {}\n\t    Copyright 2010-22 Felix Krueger, Altos Bioinformatics\n\t            https://github.com/FelixKrueger/Bismark\n\t\t\n",
            BISMARK_VERSION
        );
        return Ok(());
    }

    if cli.single && cli.paired {
        bail!("Please select either -s for single-end files or -p for paired-end files, but not both at the same time!");
    }
    if cli.percentage_cutoff.is_some() && cli.consecutive {
        bail!("The options --percentage_cutoff and --consecutive are mutually exclusive. Please respecify!");
    }
    if let Some(pct) = cli.percentage_cutoff {
        if pct > 100 {
            bail!("The percentage cutoff value has to be within the range of 0-100 [%]. Please respecify!");
        }
    }

    let samtools = find_samtools(cli.samtools_path.as_deref())?;

    let start = Instant::now();

    if !cli.single && !cli.paired {
        eprintln!("\nNeither -s (single-end) nor -p (paired-end) selected for non-bisulfite conversion filtering. Trying to extract this information for each file separately from the @PG line of the SAM/BAM file");
    }

    for file in &cli.files {
        process_file(
            file,
            cli.single,
            cli.paired,
            cli.threshold,
            cli.percentage_cutoff,
            cli.minimum_count,
            cli.consecutive,
            &samtools,
        )
        .with_context(|| format!("processing {}", file.display()))?;
    }

    eprintln!("Please continue with deduplication or methylation extraction now (depending on your application)\n");

    let elapsed = start.elapsed().as_secs();
    let days = elapsed / 86400;
    let hours = (elapsed % 86400) / 3600;
    let mins = (elapsed % 3600) / 60;
    let secs = elapsed % 60;
    eprintln!("filter_non_conversion completed in {days}d {hours}h {mins}m {secs}s\n");

    Ok(())
}

fn process_file(
    path: &Path,
    global_single: bool,
    global_paired: bool,
    threshold: u32,
    percentage_cutoff: Option<u32>,
    minimum_count: u32,
    consecutive: bool,
    samtools: &str,
) -> Result<()> {
    if !path.to_string_lossy().ends_with(".bam") {
        bail!("Please provide a BAM file to continue!");
    }

    if bam_is_truncated(samtools, path) {
        bail!("File {} appears truncated — please re-check", path.display());
    }

    let (is_single, is_paired) = if global_single {
        (true, false)
    } else if global_paired {
        (false, true)
    } else {
        detect_library_type(samtools, path)?
    };

    if !is_single && !is_paired {
        bail!("Please specify either -s (single-end) or -p (paired-end), or provide a SAM/BAM file with an @PG header line");
    }

    let name = path.display().to_string();

    if is_single {
        if let Some(pct) = percentage_cutoff {
            eprintln!("Using an overall percentage of >> {pct}% << and a minimum count of >> {minimum_count} << cytosines in non-CG context as filtering criteria before a read gets removed\n");
        } else {
            eprintln!("Using a threshold of >> {threshold} << methylation calls in non-CG context before a read gets removed\n");
        }
    } else if let Some(pct) = percentage_cutoff {
        eprintln!("Using an overall percentage of >> {pct}% << and a minimum count of >> {minimum_count} << cytosines in non-CG context as filtering criteria before a read pair gets removed (either read can fail the entire read pair)\n");
    } else {
        eprintln!("Using a threshold of >> {threshold} << methylation calls in non-CG context before a read pair gets removed (either read can fail the entire read pair)\n");
    }

    let stem = path.to_str().unwrap().trim_end_matches(".bam");
    let filtered_path = PathBuf::from(format!("{stem}.nonCG_filtered.bam"));
    let removed_path  = PathBuf::from(format!("{stem}.nonCG_removed_seqs.bam"));
    let report_path   = PathBuf::from(format!("{stem}.non-conversion_filtering.txt"));

    let mut reader  = BamReader::open(samtools, path, &[])?;
    let mut out     = BamWriter::open(samtools, &filtered_path)?;
    let mut removed = BamWriter::open(samtools, &removed_path)?;
    let mut report  = std::fs::File::create(&report_path)
        .with_context(|| format!("failed to open report {}", report_path.display()))?;

    let mut count: u64 = 0;
    let mut kicked: u64 = 0;
    let mut line_buf: Vec<u8> = Vec::with_capacity(4096);

    loop {
        line_buf.clear();
        let n = reader.lines().read_until(b'\n', &mut line_buf)?;
        if n == 0 {
            break;
        }
        while line_buf.last() == Some(&b'\n') || line_buf.last() == Some(&b'\r') {
            line_buf.pop();
        }

        if line_buf.starts_with(b"@") {
            out.write_line(&line_buf)?;
            removed.write_line(&line_buf)?;
            continue;
        }

        if is_single {
            count += 1;
            if evaluate_line(&line_buf, threshold, percentage_cutoff, minimum_count, consecutive) {
                kicked += 1;
                removed.write_line(&line_buf)?;
            } else {
                out.write_line(&line_buf)?;
            }
        } else {
            let r1 = line_buf.clone();

            // Read R2
            line_buf.clear();
            let n2 = reader.lines().read_until(b'\n', &mut line_buf)?;
            if n2 == 0 {
                out.write_line(&r1)?;
                break;
            }
            while line_buf.last() == Some(&b'\n') || line_buf.last() == Some(&b'\r') {
                line_buf.pop();
            }
            let r2 = line_buf.clone();
            count += 1;

            let fails_r1 = evaluate_line(&r1, threshold, percentage_cutoff, minimum_count, consecutive);
            let fails_r2 = if !fails_r1 {
                evaluate_line(&r2, threshold, percentage_cutoff, minimum_count, consecutive)
            } else {
                false
            };

            if fails_r1 || fails_r2 {
                kicked += 1;
                removed.write_line(&r1)?;
                removed.write_line(&r2)?;
            } else {
                out.write_line(&r1)?;
                out.write_line(&r2)?;
            }
        }
    }

    out.finish()?;
    removed.finish()?;
    reader.finish()?;

    let percent = if count == 0 {
        "N/A".to_string()
    } else {
        format!("{:.1}", kicked as f64 / count as f64 * 100.0)
    };
    let insert = if consecutive { "consecutive " } else { "" };

    eprintln!("NON-CONVERSION SUMMARY\n======================");

    let (summary_line, filter_line) = if is_paired {
        let s = format!("Analysed read pairs (paired-end) in file >> {name} <<  in total:\t{count}");
        let f = if let Some(pct) = percentage_cutoff {
            format!("Sequences removed because of apparent non-bisulfite conversion (at least {pct}% methylation and {minimum_count} non-CG calls in total in at least one of the reads):\t{kicked} ({percent}%)\n")
        } else {
            format!("Sequences removed because of apparent non-bisulfite conversion (at least {threshold} {insert}non-CG calls in one of the reads):\t{kicked} ({percent}%)\n")
        };
        (s, f)
    } else {
        let s = format!("Analysed sequences (single-end) in file >> {name} << in total:\t{count}");
        let f = if let Some(pct) = percentage_cutoff {
            format!("Sequences removed because of apparent non-bisulfite conversion (at least {pct}% methylation and {minimum_count} non-CG calls in total per read):\t{kicked} ({percent}%)\n")
        } else {
            format!("Sequences removed because of apparent non-bisulfite conversion (at least {threshold} {insert}non-CG calls per read):\t{kicked} ({percent}%)\n")
        };
        (s, f)
    };

    eprintln!("{summary_line}");
    eprintln!("{filter_line}");
    writeln!(report, "{summary_line}")?;
    writeln!(report, "{filter_line}")?;

    let elapsed = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default();
    let _ = elapsed; // runtime is reported globally in main()

    Ok(())
}

/// Returns true if the read should be filtered out.
fn evaluate_line(
    line: &[u8],
    threshold: u32,
    percentage_cutoff: Option<u32>,
    minimum_count: u32,
    consecutive: bool,
) -> bool {
    let xm = match extract_xm(line) {
        Some(x) => x,
        None => return false,
    };

    let mut noncpg_meth: u32 = 0; // H or X (methylated non-CG)
    let mut total_noncg: u32 = 0; // H, X, h, x

    for &b in xm {
        match b {
            b'H' | b'X' => {
                noncpg_meth += 1;
                total_noncg += 1;
            }
            b'h' | b'x' => {
                total_noncg += 1;
            }
            _ => {}
        }

        if consecutive && matches!(b, b'z' | b'h' | b'x') {
            noncpg_meth = 0;
        }

        if percentage_cutoff.is_none() && noncpg_meth >= threshold {
            return true;
        }
    }

    if let Some(pct_cutoff) = percentage_cutoff {
        if total_noncg >= minimum_count {
            let pct: f64 = format!("{:.1}", noncpg_meth as f64 / total_noncg as f64 * 100.0)
                .parse()
                .unwrap_or(0.0);
            if pct as u32 >= pct_cutoff {
                return true;
            }
        }
    }

    false
}

/// Extract the bytes of the XM:Z: tag from a raw SAM line.
fn extract_xm(line: &[u8]) -> Option<&[u8]> {
    const MARKER: &[u8] = b"XM:Z:";
    let pos = line.windows(MARKER.len()).position(|w| w == MARKER)?;
    let start = pos + MARKER.len();
    let end = line[start..]
        .iter()
        .position(|&b| b == b'\t')
        .map(|p| start + p)
        .unwrap_or(line.len());
    Some(&line[start..end])
}

fn detect_library_type(samtools: &str, path: &Path) -> Result<(bool, bool)> {
    eprintln!("Trying to determine the type of mapping from the SAM header line");
    let output = std::process::Command::new(samtools)
        .args(["view", "-H"])
        .arg(path)
        .output()
        .context("samtools view -H failed")?;

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
