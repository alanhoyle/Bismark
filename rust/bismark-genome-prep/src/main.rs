use std::collections::HashSet;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{bail, Context, Result};
use bismark_lib::BISMARK_VERSION;
use clap::Parser;
use flate2::read::MultiGzDecoder;

#[derive(Parser)]
#[command(
    name = "bismark_genome_preparation",
    about = "Prepare bisulfite-converted genome indexes for Bismark",
    version = BISMARK_VERSION,
    disable_version_flag = true,
)]
struct Cli {
    /// Folder containing the reference genome FASTA files
    genome_folder: PathBuf,

    /// Use Bowtie 2 for indexing (default)
    #[arg(long = "bowtie2")]
    bowtie2: bool,

    /// Use HISAT2 for indexing
    #[arg(long = "hisat2")]
    hisat2: bool,

    /// Use Minimap2 for indexing
    #[arg(long = "minimap2", alias = "mm2")]
    minimap2: bool,

    /// Write each chromosome to an individual FASTA file instead of MFA
    #[arg(long = "single_fasta")]
    single_fasta: bool,

    /// Number of threads per indexer process
    #[arg(long = "parallel")]
    parallel: Option<u32>,

    /// Path to the aligner binary directory
    #[arg(long = "path_to_aligner")]
    path_to_aligner: Option<String>,

    /// Force large index for Bowtie 2/HISAT2
    #[arg(long = "large-index")]
    large_index: bool,

    /// SLAM-seq mode: convert T→C (CT) and A→G (GA) instead of C→T / G→A
    #[arg(long = "slam")]
    slam: bool,

    /// Only calculate genomic nucleotide composition (no indexing)
    #[arg(long = "genomic_composition")]
    genomic_composition: bool,

    /// Verbose output
    #[arg(long = "verbose")]
    verbose: bool,

    /// Print version and exit
    #[arg(long = "version")]
    version: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.version {
        println!(
            "\n          Bismark - Bisulfite Mapper and Methylation Caller.\n\n          Bismark Genome Preparation Version: {}\n        Copyright 2010-25, Felix Krueger, Altos Bioinformatics\n              \n               https://github.com/FelixKrueger/Bismark\n",
            BISMARK_VERSION
        );
        return Ok(());
    }

    // Validate aligner selection
    let n_aligners = cli.bowtie2 as u8 + cli.hisat2 as u8 + cli.minimap2 as u8;
    if n_aligners > 1 {
        bail!("You may not select more than one aligner. Please make your pick! (default is Bowtie 2)");
    }
    if cli.single_fasta && cli.minimap2 {
        bail!("Minimap2 mode does not work in conjunction with --single_fasta. Please respecify!");
    }
    if cli.slam && cli.minimap2 {
        bail!("Minimap2 mode does not work in conjunction with --slam. Please respecify!");
    }
    if cli.large_index && cli.minimap2 {
        bail!("Minimap2 mode does not work in conjunction with --large_index. Please drop one option and respecify!");
    }

    // Default to bowtie2 if nothing specified
    let use_bowtie2 = cli.bowtie2 || (!cli.hisat2 && !cli.minimap2);
    let use_hisat2 = cli.hisat2;
    let use_minimap2 = cli.minimap2;

    if let Some(p) = cli.parallel {
        if p < 2 {
            bail!("--parallel should have a value of 2 or more. Please respecify");
        }
    }

    let genome_folder = cli.genome_folder.canonicalize()
        .with_context(|| format!("Cannot access genome folder: {}", cli.genome_folder.display()))?;

    let multi_fasta = !cli.single_fasta;
    if cli.single_fasta {
        eprintln!("Writing individual genomes out into single-entry fasta files (one per chromosome)\n");
    } else {
        eprintln!("Writing bisulfite genomes out into a single MFA (multi FastA) file\n");
    }

    if cli.slam {
        eprintln!("Genome will be generated and indexed in with in-silico T->C transitions, and NOT in BISULFITE MODE");
    }

    if cli.large_index {
        eprintln!("Large-index specified. Forcing generated index to be 'large', even if reference has fewer than 4 billion nucleotides.\n");
    }

    // Step I: create folders and collect FASTA filenames
    let (ct_dir, ga_dir, fasta_files) = create_bisulfite_genome_folders(&genome_folder, cli.verbose)?;

    // Step II: convert and write sequences
    convert_genome(&genome_folder, &fasta_files, &ct_dir, &ga_dir, multi_fasta, cli.slam, cli.verbose)?;

    // Step III: launch indexer
    let aligner = resolve_aligner_path(cli.path_to_aligner.as_deref(), use_bowtie2, use_hisat2, use_minimap2);
    let large_flag = if cli.large_index { "--large-index" } else { "" };
    let threads_flag = cli.parallel.map(|p| p.to_string());

    if use_bowtie2 {
        eprintln!("Bismark Genome Preparation - Step III: Launching the Bowtie 2 indexer");
    } else if use_minimap2 {
        eprintln!("Bismark Genome Preparation - Step III: Launching the Minimap2 indexing process");
    } else {
        eprintln!("Bismark Genome Preparation - Step III: Launching the HISAT2 indexer");
    }
    eprintln!("Please be aware that this process can - depending on genome size - take several hours!");

    launch_indexers_parallel(
        &ct_dir, &ga_dir, &aligner,
        large_flag, threads_flag.as_deref(),
        use_minimap2, cli.verbose,
    )?;

    Ok(())
}

fn create_bisulfite_genome_folders(
    genome_folder: &Path,
    verbose: bool,
) -> Result<(PathBuf, PathBuf, Vec<PathBuf>)> {
    if verbose {
        eprintln!("Bismark Genome Preparation - Step I: Preparing folders\n");
    }

    // Collect FASTA files
    let fasta_files = find_fasta_files(genome_folder)?;
    if fasta_files.is_empty() {
        bail!(
            "The specified genome folder {} does not contain any sequence files in FastA format (with .fa, .fa.gz, .fasta or .fasta.gz file extensions)",
            genome_folder.display()
        );
    }

    eprintln!("Bisulfite Genome Indexer version {} (last modified: 19 May 2022)", BISMARK_VERSION);

    let bisulfite_dir = genome_folder.join("Bisulfite_Genome");
    if bisulfite_dir.exists() {
        println!("\nA directory called {} already exists. Already existing converted sequences and/or already existing Bowtie 2, HISAT2 or Minimap2) indices will be overwritten!\n", bisulfite_dir.display());
    } else {
        std::fs::create_dir(&bisulfite_dir)
            .with_context(|| format!("Unable to create directory {}", bisulfite_dir.display()))?;
        if verbose {
            eprintln!("Created Bisulfite Genome folder {}", bisulfite_dir.display());
        }
    }

    let ct_dir = bisulfite_dir.join("CT_conversion");
    let ga_dir = bisulfite_dir.join("GA_conversion");

    for dir in [&ct_dir, &ga_dir] {
        if !dir.exists() {
            std::fs::create_dir(dir)
                .with_context(|| format!("Unable to create directory {}", dir.display()))?;
            if verbose {
                eprintln!("Created Bisulfite Genome folder {}", dir.display());
            }
        }
    }

    eprintln!("\nStep I - Prepare genome folders - completed\n\n");
    Ok((ct_dir, ga_dir, fasta_files))
}

fn find_fasta_files(genome_folder: &Path) -> Result<Vec<PathBuf>> {
    for ext in &["*.fa", "*.fa.gz", "*.fasta", "*.fasta.gz"] {
        let mut files: Vec<PathBuf> = std::fs::read_dir(genome_folder)
            .with_context(|| format!("Cannot read {}", genome_folder.display()))?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| {
                let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("");
                match *ext {
                    "*.fa" => name.ends_with(".fa") && !name.ends_with(".fa.gz"),
                    "*.fa.gz" => name.ends_with(".fa.gz"),
                    "*.fasta" => name.ends_with(".fasta") && !name.ends_with(".fasta.gz"),
                    "*.fasta.gz" => name.ends_with(".fasta.gz"),
                    _ => false,
                }
            })
            .collect();
        files.sort();
        if !files.is_empty() {
            return Ok(files);
        }
    }
    Ok(Vec::new())
}

fn open_fasta(path: &Path) -> Result<Box<dyn BufRead>> {
    let f = std::fs::File::open(path)
        .with_context(|| format!("opening {}", path.display()))?;
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    if name.ends_with(".gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(f))))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}

fn convert_genome(
    _genome_folder: &Path,
    fasta_files: &[PathBuf],
    ct_dir: &Path,
    ga_dir: &Path,
    multi_fasta: bool,
    slam: bool,
    verbose: bool,
) -> Result<()> {
    if verbose {
        eprintln!("Bismark Genome Preparation - Step II: Bisulfite converting reference genome\n");
    }

    let mut total_ct: u64 = 0;
    let mut total_ga: u64 = 0;
    let mut seen_chrs: HashSet<String> = HashSet::new();

    let suffix_ct = if slam { "_TC_converted" } else { "_CT_converted" };
    let suffix_ga = if slam { "_AG_converted" } else { "_GA_converted" };

    // Use Box<dyn Write> so we can replace per-chromosome writers
    let mut ct_out: Box<dyn Write>;
    let mut ga_out: Box<dyn Write>;

    if multi_fasta {
        let ct_path = ct_dir.join("genome_mfa.CT_conversion.fa");
        let ga_path = ga_dir.join("genome_mfa.GA_conversion.fa");
        ct_out = Box::new(BufWriter::new(std::fs::File::create(&ct_path)
            .with_context(|| format!("creating {}", ct_path.display()))?));
        ga_out = Box::new(BufWriter::new(std::fs::File::create(&ga_path)
            .with_context(|| format!("creating {}", ga_path.display()))?));
    } else {
        // Placeholder; will be replaced before use
        ct_out = Box::new(std::io::sink());
        ga_out = Box::new(std::io::sink());
    }

    for fasta_path in fasta_files {
        let mut reader = open_fasta(fasta_path)?;
        let mut line = String::new();

        // Read first line — must be a FASTA header
        reader.read_line(&mut line)?;
        let first = line.trim_end_matches(['\n', '\r']).to_string();
        if !first.starts_with('>') {
            bail!("The specified file ({}) doesn't seem to be in FASTA format as required!", fasta_path.display());
        }
        let mut chr_name = extract_chr_name(&first);

        if !seen_chrs.insert(chr_name.clone()) {
            bail!("Exiting because chromosome name '{}' already exists. Please make sure all chromosomes have a unique name!", chr_name);
        }

        if !multi_fasta {
            let ct_path = ct_dir.join(format!("{}.CT_conversion.fa", chr_name));
            let ga_path = ga_dir.join(format!("{}.GA_conversion.fa", chr_name));
            ct_out = Box::new(BufWriter::new(std::fs::File::create(&ct_path)
                .with_context(|| format!("creating {}", ct_path.display()))?));
            ga_out = Box::new(BufWriter::new(std::fs::File::create(&ga_path)
                .with_context(|| format!("creating {}", ga_path.display()))?));
        }

        writeln!(ct_out, ">{chr_name}{suffix_ct}")?;
        writeln!(ga_out, ">{chr_name}{suffix_ga}")?;

        loop {
            line.clear();
            let n = reader.read_line(&mut line)?;
            if n == 0 { break; }
            let trimmed = line.trim_end_matches(['\n', '\r']).to_string();

            if trimmed.starts_with('>') {
                chr_name = extract_chr_name(&trimmed);
                if !seen_chrs.insert(chr_name.clone()) {
                    bail!("Exiting because chromosome name '{}' already exists. Please make sure all chromosomes have a unique name!", chr_name);
                }
                if !multi_fasta {
                    let ct_path = ct_dir.join(format!("{}.CT_conversion.fa", chr_name));
                    let ga_path = ga_dir.join(format!("{}.GA_conversion.fa", chr_name));
                    ct_out = Box::new(BufWriter::new(std::fs::File::create(&ct_path)
                        .with_context(|| format!("creating {}", ct_path.display()))?));
                    ga_out = Box::new(BufWriter::new(std::fs::File::create(&ga_path)
                        .with_context(|| format!("creating {}", ga_path.display()))?));
                }
                writeln!(ct_out, ">{chr_name}{suffix_ct}")?;
                writeln!(ga_out, ">{chr_name}{suffix_ga}")?;
            } else {
                let seq: Vec<u8> = trimmed.bytes().map(|b| {
                    let u = b.to_ascii_uppercase();
                    match u {
                        b'A' | b'T' | b'C' | b'G' | b'N' => u,
                        _ => b'N',
                    }
                }).collect();

                let (mut ct_seq, mut ga_seq) = (seq.clone(), seq);

                let ct_n: u64 = if slam {
                    ct_seq.iter_mut().filter(|b| **b == b'T').map(|b| { *b = b'C'; }).count() as u64
                } else {
                    ct_seq.iter_mut().filter(|b| **b == b'C').map(|b| { *b = b'T'; }).count() as u64
                };
                let ga_n: u64 = if slam {
                    ga_seq.iter_mut().filter(|b| **b == b'A').map(|b| { *b = b'G'; }).count() as u64
                } else {
                    ga_seq.iter_mut().filter(|b| **b == b'G').map(|b| { *b = b'A'; }).count() as u64
                };

                total_ct += ct_n;
                total_ga += ga_n;

                ct_out.write_all(&ct_seq)?;
                ct_out.write_all(b"\n")?;
                ga_out.write_all(&ga_seq)?;
                ga_out.write_all(b"\n")?;
            }
        }
    }

    println!("\nTotal number of conversions performed:");
    if slam {
        println!("T->C:\t{total_ct}");
        println!("A->G:\t{total_ga}");
        eprintln!("\nStep II - Genome SLAM conversions - completed\n\n");
    } else {
        println!("C->T:\t{total_ct}");
        println!("G->A:\t{total_ga}");
        eprintln!("\nStep II - Genome bisulfite conversions - completed\n\n");
    }

    Ok(())
}

fn extract_chr_name(header: &str) -> String {
    let without_gt = header.trim_start_matches('>');
    without_gt.split_whitespace().next().unwrap_or(without_gt).to_string()
}

fn resolve_aligner_path(
    user_path: Option<&str>,
    _bowtie2: bool,
    hisat2: bool,
    minimap2: bool,
) -> String {
    let binary = if minimap2 { "minimap2" } else if hisat2 { "hisat2-build" } else { "bowtie2-build" };
    if let Some(p) = user_path {
        let p = if p.ends_with('/') { p.to_string() } else { format!("{p}/") };
        format!("{p}{binary}")
    } else {
        binary.to_string()
    }
}

fn launch_indexers_parallel(
    ct_dir: &Path,
    ga_dir: &Path,
    aligner: &str,
    large_flag: &str,
    threads: Option<&str>,
    minimap2: bool,
    verbose: bool,
) -> Result<()> {
    let ct_dir = ct_dir.to_path_buf();
    let ga_dir = ga_dir.to_path_buf();
    let aligner = aligner.to_string();
    let large_flag = large_flag.to_string();
    let _threads_flag = threads.map(|t| t.to_string());

    // Run CT and GA indexing in parallel using threads
    let ct_dir_clone = ct_dir.clone();
    let aligner_ct = aligner.clone();
    let large_ct = large_flag.clone();
    let threads_ct = threads.map(|t| t.to_string());
    let minimap2_c = minimap2;
    let verbose_c = verbose;

    let ct_handle = std::thread::spawn(move || {
        let run = |dir: PathBuf| -> Result<()> {
            let fa_files: Vec<String> = std::fs::read_dir(&dir)?
                .filter_map(|e| e.ok())
                .filter(|e| {
                    let name = e.file_name();
                    let s = name.to_string_lossy();
                    s.ends_with(".fa") || s.ends_with(".fasta")
                })
                .map(|e| e.file_name().to_string_lossy().into_owned())
                .collect();
            if fa_files.is_empty() {
                bail!("No .fa files found in {}", dir.display());
            }
            let file_list = fa_files.join(",");
            let mut cmd = Command::new(&aligner_ct);
            cmd.current_dir(&dir);
            if minimap2_c {
                cmd.arg("-k").arg("20");
                if let Some(ref t) = threads_ct { cmd.arg("-t").arg(t); }
                cmd.arg("-d").arg("BS_CT.mmi").arg(&file_list);
            } else {
                if let Some(ref t) = threads_ct { cmd.arg("--threads").arg(t); }
                if !large_ct.is_empty() { cmd.arg(&large_ct); }
                cmd.arg("-f").arg(&file_list).arg("BS_CT");
            }
            if verbose_c {
                eprintln!("Parent process: Starting to index C->T converted genome: {:?}", cmd);
            }
            let status = cmd.status().context("launching CT indexer")?;
            if !status.success() { bail!("CT index build failed"); }
            Ok(())
        };
        run(ct_dir_clone)
    });

    // GA indexing in current thread
    {
        let fa_files: Vec<String> = std::fs::read_dir(&ga_dir)?
            .filter_map(|e| e.ok())
            .filter(|e| {
                let name = e.file_name();
                let s = name.to_string_lossy();
                s.ends_with(".fa") || s.ends_with(".fasta")
            })
            .map(|e| e.file_name().to_string_lossy().into_owned())
            .collect();

        if fa_files.is_empty() {
            bail!("No .fa files found in {}", ga_dir.display());
        }
        let file_list = fa_files.join(",");
        let mut cmd = Command::new(&aligner);
        cmd.current_dir(&ga_dir);
        if minimap2 {
            cmd.arg("-k").arg("20");
            if let Some(t) = threads { cmd.arg("-t").arg(t); }
            cmd.arg("-d").arg("BS_GA.mmi").arg(&file_list);
        } else {
            if let Some(t) = threads { cmd.arg("--threads").arg(t); }
            if !large_flag.is_empty() { cmd.arg(&large_flag); }
            cmd.arg("-f").arg(&file_list).arg("BS_GA");
        }
        if verbose {
            eprintln!("Child process: Starting to index G->A converted genome: {:?}", cmd);
        }
        let status = cmd.status().context("launching GA indexer")?;
        if !status.success() {
            bail!("GA index build failed");
        }
    }

    let ct_result = ct_handle.join().unwrap_or_else(|_| Err(anyhow::anyhow!("CT indexer thread panicked")));
    ct_result?;

    eprintln!("\n=========================================\n\nParallel genome indexing complete. Enjoy!\n");
    Ok(())
}
