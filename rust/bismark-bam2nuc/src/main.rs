use std::collections::HashMap;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use bismark_lib::bam_io::{find_samtools, BamReader};
use bismark_lib::fasta::Genome;
use bismark_lib::BISMARK_VERSION;
use clap::Parser;

const MONO: [&str; 4] = ["A", "C", "G", "T"];
const DI: [&str; 16] = [
    "AA","AC","AG","AT","CA","CC","CG","CT",
    "GA","GC","GG","GT","TA","TC","TG","TT",
];

#[derive(Parser)]
#[command(
    name = "bam2nuc",
    about = "Calculates nucleotide coverage from Bismark BAM/CRAM files",
    version = BISMARK_VERSION,
    disable_version_flag = true,
)]
struct Cli {
    /// BAM/CRAM file(s) to process
    files: Vec<PathBuf>,

    /// Directory to write output files into
    #[arg(long = "dir", default_value = "")]
    output_dir: String,

    /// Folder containing the reference genome FASTA files
    #[arg(short = 'g', long = "genome_folder", required = true)]
    genome_folder: PathBuf,

    /// Only calculate genomic nucleotide frequencies, skip BAM processing
    #[arg(long)]
    genomic_composition_only: bool,

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
            "\n\n                        Bismark Nucleotide Coverage Module -\n                                     bam2nuc\n\n                           Bismark Version: {}\n              Copyright 2010-22 Felix Krueger, Altos Bioinformatics\n                     https://github.com/FelixKrueger/Bismark\n\t       \n\n",
            BISMARK_VERSION
        );
        return Ok(());
    }

    if cli.files.is_empty() && !cli.genomic_composition_only {
        bail!("You need to provide one or more BAM files to continue. Please respecify!");
    }

    let samtools = find_samtools(cli.samtools_path.as_deref())?;
    let output_dir = normalise_dir(&cli.output_dir);

    eprintln!("Summary of parameters for nucleotide coverage report:");
    eprintln!("{}", "=".repeat(53));
    eprintln!("Output directory:\t\t\t>{output_dir}<");
    eprintln!("Genome directory:\t\t\t>{}<", cli.genome_folder.display());
    eprintln!("Samtools installation:\t\t\t>{samtools}<\n");

    let genome = Genome::load(&cli.genome_folder)?;
    eprintln!("Stored sequence information of {} chromosomes/scaffolds in total\n", genome.sequences.len());

    if cli.genomic_composition_only {
        get_genomic_frequencies(&genome, &cli.genome_folder, &output_dir)?;
        eprintln!("Finished processing genomic nucleotide frequencies\n");
        return Ok(());
    }

    for file in &cli.files {
        generate_nucleotide_report(file, &genome, &cli.genome_folder, &output_dir, &samtools)?;
    }
    Ok(())
}

fn generate_nucleotide_report(
    infile: &Path,
    genome: &Genome,
    genome_folder: &Path,
    output_dir: &str,
    samtools: &str,
) -> Result<()> {
    eprintln!("{}", "=".repeat(66));
    eprintln!("Mono- and di-nucleotide coverage will now be written into a report");
    eprintln!("{}\n", "=".repeat(66));

    let genomic_freqs = get_genomic_frequencies(genome, genome_folder, output_dir)?;

    eprintln!("\n\nCalculating read frequencies from file '{}'", infile.display());
    eprintln!("{}", "=".repeat(90));

    let is_paired = detect_is_paired(samtools, infile)?;
    if is_paired {
        eprintln!("Determined the file to be paired-end");
    } else {
        eprintln!("Determined the file to be single-end");
    }

    let mut read_freqs: HashMap<String, u64> = HashMap::new();
    let mut reader = BamReader::open(samtools, infile, &[])?;
    let mut buf: Vec<u8> = Vec::with_capacity(4096);
    let mut count: u64 = 0;

    loop {
        buf.clear();
        let n = reader.lines().read_until(b'\n', &mut buf)?;
        if n == 0 { break; }
        while buf.last() == Some(&b'\n') || buf.last() == Some(&b'\r') { buf.pop(); }
        if buf.starts_with(b"@") || buf.is_empty() { continue; }

        count += 1;
        if count % 500_000 == 0 { eprintln!("Processed {count} lines"); }

        let fields: Vec<&[u8]> = buf.splitn(12, |&b| b == b'\t').collect();
        if fields.len() < 10 { continue; }

        let flag: u16 = parse_u16(fields[1]);
        let chr = std::str::from_utf8(fields[2]).unwrap_or("*");
        let start: usize = parse_usize(fields[3]);
        let cigar = fields[5];
        let seq_len = fields[9].len();

        // Skip reads with indels/splicing/soft-clip in CIGAR
        if cigar.iter().any(|&b| matches!(b, b'I' | b'D' | b'S' | b'N')) {
            continue;
        }

        let genomic_seq = match genome.slice(chr, start.saturating_sub(1), seq_len) {
            Some(s) => s,
            None => { continue; }
        };

        let seq: Vec<u8> = if needs_revcomp(flag, is_paired) {
            rev_comp(genomic_seq)
        } else {
            genomic_seq.to_vec()
        };

        count_nucleotides(&seq, &mut read_freqs);
    }
    reader.finish()?;
    eprintln!("\n\n");

    let fname = infile.file_name().unwrap().to_string_lossy();
    let outname = if fname.ends_with(".bam") {
        fname.replace(".bam", ".nucleotide_stats.txt")
    } else if fname.ends_with(".cram") {
        fname.replace(".cram", ".nucleotide_stats.txt")
    } else {
        bail!("File needs to be in BAM or CRAM format (ending in .bam or .cram). Terminating process...");
    };
    let outpath = format!("{output_dir}{outname}");
    eprintln!("Printing nucleotide stats to >> {outpath} <<");
    let mut out = std::fs::File::create(&outpath).with_context(|| format!("failed to create {outpath}"))?;

    write_nucleotide_stats(&mut out, &read_freqs, &genomic_freqs)?;
    Ok(())
}

fn get_genomic_frequencies(
    genome: &Genome,
    genome_folder: &Path,
    output_dir: &str,
) -> Result<HashMap<String, u64>> {
    let cached = genome_folder.join("genomic_nucleotide_frequencies.txt");
    if cached.exists() {
        eprintln!("Detected file 'genomic_nucleotide_frequencies.txt' in the genome folder {}. Using nucleotide frequencies contained therein ...", genome_folder.display());
        eprintln!("{}", "=".repeat(188));
        let mut map = HashMap::new();
        for line in std::fs::read_to_string(&cached)?.lines() {
            let mut parts = line.splitn(2, '\t');
            if let (Some(k), Some(v)) = (parts.next(), parts.next()) {
                if let Ok(n) = v.parse::<u64>() { map.insert(k.to_string(), n); }
            }
        }
        return Ok(map);
    }

    eprintln!("Could not find genomic nucleotide frequency table in the genome folder, calculating genomic frequencies (this may take several minutes depending on genome size) ...");
    eprintln!("{}", "=".repeat(164));
    let mut freqs: HashMap<String, u64> = HashMap::new();
    for (chr, seq) in &genome.sequences {
        eprintln!("Processing chromosome >> {chr} <<");
        count_nucleotides(seq, &mut freqs);
    }

    // Try to write cache
    let mut keys: Vec<&String> = freqs.keys().collect();
    keys.sort();
    let write_cache = |path: &str| -> std::io::Result<()> {
        let mut f = std::fs::File::create(path)?;
        for k in &keys { writeln!(f, "{k}\t{}", freqs[*k])?; }
        Ok(())
    };
    if write_cache(&cached.to_string_lossy()).is_ok() {
        eprintln!("Writing genomic nucleotide frequencies to the file >{}<", cached.display());
    } else {
        let alt = format!("{output_dir}genomic_nucleotide_frequencies.txt");
        if write_cache(&alt).is_ok() {
            eprintln!("Writing genomic nucleotide frequencies to the file >{alt}<");
        }
    }

    Ok(freqs)
}

fn count_nucleotides(seq: &[u8], freqs: &mut HashMap<String, u64>) {
    for i in 0..seq.len() {
        if seq[i] != b'N' {
            *freqs.entry((seq[i] as char).to_string()).or_default() += 1;
        }
        if i + 1 < seq.len() {
            let di = [seq[i], seq[i + 1]];
            if !di.contains(&b'N') {
                *freqs.entry(String::from_utf8_lossy(&di).into_owned()).or_default() += 1;
            }
        }
    }
}

fn write_nucleotide_stats(
    out: &mut std::fs::File,
    read_freqs: &HashMap<String, u64>,
    genomic_freqs: &HashMap<String, u64>,
) -> Result<()> {
    eprintln!("Final Stage: Calculating averages\n{}", "=".repeat(33));
    let header = "(di-)nucleotide\tcount sample\tpercent sample\tcount genomic\tpercent genomic\tcoverage";
    eprintln!("{header}");
    writeln!(out, "{header}")?;

    let mono_r: u64 = MONO.iter().map(|k| read_freqs.get(*k).copied().unwrap_or(0)).sum();
    let mono_g: u64 = MONO.iter().map(|k| genomic_freqs.get(*k).copied().unwrap_or(0)).sum();
    for word in MONO {
        let rc = read_freqs.get(word).copied().unwrap_or(0);
        let gc = genomic_freqs.get(word).copied().unwrap_or(0);
        let line = stat_line(word, rc, gc, mono_r, mono_g);
        eprintln!("{line}"); writeln!(out, "{line}")?;
    }

    let di_r: u64 = DI.iter().map(|k| read_freqs.get(*k).copied().unwrap_or(0)).sum();
    let di_g: u64 = DI.iter().map(|k| genomic_freqs.get(*k).copied().unwrap_or(0)).sum();
    for word in DI {
        let rc = read_freqs.get(word).copied().unwrap_or(0);
        let gc = genomic_freqs.get(word).copied().unwrap_or(0);
        let line = stat_line(word, rc, gc, di_r, di_g);
        eprintln!("{line}"); writeln!(out, "{line}")?;
    }
    Ok(())
}

fn stat_line(word: &str, rc: u64, gc: u64, total_r: u64, total_g: u64) -> String {
    let rp = if total_r > 0 { format!("{:.2}", 100.0 * rc as f64 / total_r as f64) } else { "0.00".into() };
    let gp = if total_g > 0 { format!("{:.2}", 100.0 * gc as f64 / total_g as f64) } else { "0.00".into() };
    let cov = if gc > 0 { format!("{:.3}", rc as f64 / gc as f64) } else { "0.000".into() };
    format!("{word}\t{rc}\t{rp}\t{gc}\t{gp}\t{cov}")
}

fn needs_revcomp(flag: u16, is_paired: bool) -> bool {
    if is_paired { flag == 83 || flag == 163 } else { flag == 16 }
}

fn rev_comp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' => b'T', b'T' => b'A', b'G' => b'C', b'C' => b'G', _ => b'N'
    }).collect()
}

fn detect_is_paired(samtools: &str, path: &Path) -> Result<bool> {
    let output = std::process::Command::new(samtools)
        .args(["view", "-H"]).arg(path).output().context("samtools view -H")?;
    for chunk in output.stdout.split(|&b| b == b'\n') {
        if !chunk.starts_with(b"@PG") { continue; }
        let s = String::from_utf8_lossy(chunk);
        if s.contains(" -1 ") && s.contains(" -2 ") { return Ok(true); }
        return Ok(false);
    }
    Ok(false)
}

fn normalise_dir(s: &str) -> String {
    if s.is_empty() { return String::new(); }
    if s.ends_with('/') { s.to_string() } else { format!("{s}/") }
}

fn parse_u16(b: &[u8]) -> u16 {
    std::str::from_utf8(b).ok().and_then(|s| s.parse().ok()).unwrap_or(0)
}
fn parse_usize(b: &[u8]) -> usize {
    std::str::from_utf8(b).ok().and_then(|s| s.parse().ok()).unwrap_or(0)
}
