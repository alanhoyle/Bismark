mod align;
mod call;
mod convert;
mod fastq;
mod report;

use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use bismark_lib::bam_io::{find_samtools, BamWriter};
use bismark_lib::fasta::Genome;
use bismark_lib::BISMARK_VERSION;
use clap::Parser;

use align::{
    align_strand_pe, align_strand_se, strand_configs_directional, strand_configs_nondirectional,
    strand_configs_pbat, StrandConfig,
};
use call::{select_best_pe, select_best_se, AlignedRead, Outcome};
use report::{write_pe_report, write_se_report, AlignStats};

#[derive(Parser)]
#[command(
    name = "bistromark",
    about = "Bisulfite-seq aligner (Rust reimplementation of bismark)",
    version = BISMARK_VERSION,
    disable_version_flag = true,
)]
struct Cli {
    // ── Required ──────────────────────────────────────────────────────────────
    /// Genome folder containing Bisulfite_Genome/ subdirectory.
    #[arg(long = "genome", required_unless_present = "version")]
    genome_folder: Option<PathBuf>,

    // ── Input mode ────────────────────────────────────────────────────────────
    /// R1 FASTQ for paired-end mode.
    #[arg(short = '1')]
    read1: Option<PathBuf>,

    /// R2 FASTQ for paired-end mode.
    #[arg(short = '2')]
    read2: Option<PathBuf>,

    /// FASTQ for single-end mode.
    #[arg(short = 'U', alias = "se")]
    unpaired: Option<PathBuf>,

    /// Input is FASTA format (quality scores are set to Phred 40).
    #[arg(short = 'f', long = "fasta")]
    fasta: bool,

    // ── Bisulfite strand mode ─────────────────────────────────────────────────
    /// Non-directional library: use all 4 bisulfite strands.
    #[arg(long = "non_directional")]
    non_directional: bool,

    /// PBAT mode: GA-converted reads aligned to both genomes.
    #[arg(long = "pbat")]
    pbat: bool,

    // ── Output ────────────────────────────────────────────────────────────────
    /// Output directory.
    #[arg(short = 'o', long = "output_dir", default_value = ".")]
    output_dir: PathBuf,

    /// Write SAM instead of BAM.
    #[arg(long = "sam")]
    sam_output: bool,

    /// Use STEM as the output file name stem (overrides auto-derived stem).
    #[arg(short = 'B', long = "basename")]
    basename: Option<String>,

    /// Prepend PREFIX_ to the auto-derived output stem.
    #[arg(long = "prefix")]
    prefix: Option<String>,

    /// Write unmapped reads to a FASTQ file.
    #[arg(long = "unmapped", alias = "un")]
    unmapped: bool,

    /// Write ambiguously-mapped reads to a FASTQ file.
    #[arg(long = "ambiguous")]
    ambiguous: bool,

    // ── Threads / paths ───────────────────────────────────────────────────────
    /// Number of threads passed to Bowtie2 with -p.
    #[arg(short = 'p', long = "threads", default_value_t = 1)]
    threads: usize,

    /// Path to bowtie2 binary (or directory containing it).
    #[arg(long = "path_to_bowtie2")]
    path_to_bowtie2: Option<PathBuf>,

    /// Path to samtools binary.
    #[arg(long = "samtools_path")]
    samtools_path: Option<String>,

    // ── Bowtie2 alignment options (passed through) ────────────────────────────
    /// Minimum alignment score function (e.g. L,0,-0.2).
    #[arg(long = "score_min")]
    score_min: Option<String>,

    /// Max seed mismatches (bowtie2 -N; 0 or 1, default 0).
    #[arg(short = 'n', long = "seedmms")]
    seed_mms: Option<u8>,

    /// Seed length (bowtie2 -L).
    #[arg(short = 'l', long = "seedlen")]
    seed_len: Option<u32>,

    /// Max consecutive failed extends before abandoning a seed (bowtie2 -D).
    #[arg(long = "D")]
    bowtie2_d: Option<u32>,

    /// Max re-seedings for repetitive seeds (bowtie2 -R).
    #[arg(long = "R")]
    bowtie2_r: Option<u32>,

    /// Read gap open,extend penalties (bowtie2 --rdg).
    #[arg(long = "rdg")]
    rdg: Option<String>,

    /// Reference gap open,extend penalties (bowtie2 --rfg).
    #[arg(long = "rfg")]
    rfg: Option<String>,

    /// Local alignment mode (default: end-to-end).
    #[arg(long = "local")]
    local: bool,

    /// Input qualities are Phred+64 (Illumina 1.3+).
    #[arg(long = "phred64")]
    phred64: bool,

    // ── PE pairing options ────────────────────────────────────────────────────
    /// Minimum insert size for valid PE alignment (bowtie2 -I).
    #[arg(short = 'I', long = "minins", default_value_t = 0)]
    minins: u32,

    /// Maximum insert size for valid PE alignment (bowtie2 -X).
    #[arg(short = 'X', long = "maxins", default_value_t = 500)]
    maxins: u32,

    /// Disable dovetail PE alignment (dovetail is on by default).
    #[arg(long = "no_dovetail")]
    no_dovetail: bool,

    // ── Read filtering ────────────────────────────────────────────────────────
    /// Skip the first N reads.
    #[arg(short = 's', long = "skip")]
    skip: Option<usize>,

    /// Process only the first N reads.
    #[arg(short = 'u', long = "upto")]
    upto: Option<usize>,

    // ── Read group tags ───────────────────────────────────────────────────────
    /// Attach RG:Z tags using --rg_id and --rg_sample.
    #[arg(long = "rg_tag")]
    rg_tag: bool,

    /// Read group ID for RG header and RG:Z tag (default: "1").
    #[arg(long = "rg_id")]
    rg_id: Option<String>,

    /// Sample name for RG header SM field (default: input file stem).
    #[arg(long = "rg_sample")]
    rg_sample: Option<String>,

    // ── Misc ──────────────────────────────────────────────────────────────────
    /// Print version and exit.
    #[arg(long = "version")]
    version: bool,

    /// Suppress progress messages.
    #[arg(long = "quiet")]
    quiet: bool,

    // ── Compat no-ops (accepted silently) ────────────────────────────────────
    #[arg(short = 'q', long = "fastq", hide = true)]
    fastq: bool,
    #[arg(long = "end_to_end", hide = true)]
    end_to_end: bool,
    #[arg(long = "phred33", hide = true)]
    phred33: bool,
    #[arg(long = "icpc", hide = true)]
    icpc: bool,
    #[arg(long = "no_unal", hide = true)]
    no_unal: bool,
    #[arg(long = "bowtie2", hide = true)]
    use_bowtie2: bool,
    #[arg(long = "bam", hide = true)]
    bam: bool,
    #[arg(long = "gzip", hide = true)]
    gzip: bool,
    #[arg(long = "bt2_large_index", hide = true)]
    bt2_large_index: bool,

    // ── Warn-and-continue (deferred features) ────────────────────────────────
    #[arg(long = "hisat2", hide = true)]
    hisat2: bool,
    #[arg(long = "path_to_hisat2", hide = true)]
    path_to_hisat2: Option<PathBuf>,
    #[arg(long = "hisat2_options_string", hide = true)]
    hisat2_options_string: Option<String>,
    #[arg(long = "minimap2", alias = "mm2", hide = true)]
    minimap2: bool,
    #[arg(long = "path_to_minimap2", hide = true)]
    path_to_minimap2: Option<PathBuf>,
    #[arg(long = "cram", hide = true)]
    cram: bool,
    #[arg(long = "cram_ref", hide = true)]
    cram_ref: Option<PathBuf>,
    #[arg(long = "slam", hide = true)]
    slam: bool,
    #[arg(long = "parallel", alias = "multicore", hide = true)]
    parallel: Option<u32>,
    #[arg(long = "non_bs_mm", hide = true)]
    non_bs_mm: bool,
    #[arg(long = "nucleotide_coverage", hide = true)]
    nucleotide_coverage: bool,
    #[arg(long = "old_flag", hide = true)]
    old_flag: bool,
    #[arg(long = "ambig_bam", hide = true)]
    ambig_bam: bool,
    #[arg(long = "temp_dir", hide = true)]
    temp_dir: Option<PathBuf>,
    #[arg(long = "strandID", hide = true)]
    strand_id: bool,
    #[arg(long = "sam_no_hd", hide = true)]
    sam_no_hd: bool,
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {:?}", e);
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let cli = Cli::parse();

    if cli.version {
        println!("Bismark (Rust bistromark) version {}", BISMARK_VERSION);
        return Ok(());
    }

    emit_deferred_warnings(&cli);

    let genome_folder = cli.genome_folder.as_deref().unwrap();

    let is_pe = cli.read1.is_some() && cli.read2.is_some();
    let is_se = cli.unpaired.is_some();
    if !is_pe && !is_se {
        bail!("Specify -1/-2 for paired-end or -U for single-end mode");
    }
    if is_pe && is_se {
        bail!("Cannot use -U together with -1/-2");
    }

    let bowtie2 = find_bowtie2(cli.path_to_bowtie2.as_deref())?;

    let samtools = if cli.sam_output {
        String::new()
    } else {
        find_samtools(cli.samtools_path.as_deref())?
    };

    let extra_args = build_base_extra_args(&cli);

    if !cli.quiet {
        eprintln!("Loading genome from: {}", genome_folder.display());
    }
    let genome = Genome::load(genome_folder)?;

    let strand_configs: Vec<StrandConfig> = if cli.non_directional {
        strand_configs_nondirectional(genome_folder)
    } else if cli.pbat {
        strand_configs_pbat(genome_folder).into()
    } else {
        strand_configs_directional(genome_folder).into()
    };

    std::fs::create_dir_all(&cli.output_dir)
        .with_context(|| format!("cannot create output dir {}", cli.output_dir.display()))?;

    if is_se {
        run_se(
            &cli,
            genome_folder,
            &genome,
            &strand_configs,
            &bowtie2,
            &samtools,
            &extra_args,
        )
    } else {
        run_pe(
            &cli,
            genome_folder,
            &genome,
            &strand_configs,
            &bowtie2,
            &samtools,
            &extra_args,
        )
    }
}

/// Emit stderr warnings for flags that are accepted but not yet implemented.
fn emit_deferred_warnings(cli: &Cli) {
    let warn = |msg: &str| eprintln!("WARNING: {msg}");

    if cli.hisat2 || cli.path_to_hisat2.is_some() {
        warn("--hisat2 is not implemented in bistromark; falling back to Bowtie2");
    }
    if cli.minimap2 || cli.path_to_minimap2.is_some() {
        warn("--minimap2 is not implemented in bistromark; falling back to Bowtie2");
    }
    if cli.cram || cli.cram_ref.is_some() {
        warn("--cram is not implemented; writing BAM output instead");
    }
    if cli.slam {
        warn("--slam is not implemented; running standard bisulfite mode");
    }
    if let Some(n) = cli.parallel {
        warn(&format!(
            "--parallel {n} is not implemented; running single-process (use -p for Bowtie2 threads)"
        ));
    }
    if cli.non_bs_mm {
        warn("--non_bs_mm is not implemented; all mismatches treated as bisulfite");
    }
    if cli.nucleotide_coverage {
        warn("--nucleotide_coverage is not implemented; skipping");
    }
    if cli.old_flag {
        warn("--old_flag is not implemented; using current SAM FLAG encoding");
    }
    if cli.ambig_bam {
        warn("--ambig_bam is not implemented; ambiguous reads not written to BAM");
    }
    if cli.strand_id {
        warn("--strandID is not implemented; YS:Z tag will not be added");
    }
    if cli.gzip {
        warn("--gzip not implemented for bistromark unmapped/ambiguous output files; writing plain FASTQ");
    }
}

/// Build the bowtie2 arguments common to both SE and PE modes.
fn build_base_extra_args(cli: &Cli) -> Vec<String> {
    let mut args: Vec<String> = Vec::new();

    // Input format
    if cli.fasta {
        args.push("-f".to_owned());
    } else {
        args.push("-q".to_owned());
    }
    args.push("--ignore-quals".to_owned());

    // Alignment mode
    if cli.local {
        args.push("--local".to_owned());
    } else {
        args.push("--end-to-end".to_owned());
    }

    // Score minimum function
    let score_min = cli
        .score_min
        .clone()
        .unwrap_or_else(|| {
            if cli.local {
                "L,0,0.5".to_owned()
            } else {
                "L,0,-0.2".to_owned()
            }
        });
    args.push("--score-min".to_owned());
    args.push(score_min);

    // Optional bowtie2 tuning flags
    if let Some(n) = cli.seed_mms {
        args.push("-N".to_owned());
        args.push(n.to_string());
    }
    if let Some(l) = cli.seed_len {
        args.push("-L".to_owned());
        args.push(l.to_string());
    }
    if let Some(d) = cli.bowtie2_d {
        args.push("-D".to_owned());
        args.push(d.to_string());
    }
    if let Some(r) = cli.bowtie2_r {
        args.push("-R".to_owned());
        args.push(r.to_string());
    }
    if let Some(rdg) = &cli.rdg {
        args.push("--rdg".to_owned());
        args.push(rdg.clone());
    }
    if let Some(rfg) = &cli.rfg {
        args.push("--rfg".to_owned());
        args.push(rfg.clone());
    }
    if cli.phred64 {
        args.push("--phred64".to_owned());
    }

    args
}

fn run_se(
    cli: &Cli,
    genome_folder: &Path,
    genome: &Genome,
    configs: &[StrandConfig],
    bowtie2: &str,
    samtools: &str,
    extra_args: &[String],
) -> Result<()> {
    let input = cli.unpaired.as_deref().unwrap();

    if !cli.quiet {
        eprintln!("Reading SE input: {}", input.display());
    }
    let mut reads = if cli.fasta {
        fastq::read_fasta(input)?
    } else {
        fastq::read_all(input)?
    };

    // Apply --skip / --upto
    if let Some(s) = cli.skip {
        if s >= reads.len() {
            reads.clear();
        } else {
            reads.drain(..s);
        }
    }
    if let Some(n) = cli.upto {
        reads.truncate(n);
    }

    if !cli.quiet {
        eprintln!("{} reads loaded", reads.len());
    }

    let reads_seq: Vec<Vec<u8>> = reads.iter().map(|r| r.seq.clone()).collect();
    let reads_qual: Vec<Vec<u8>> = reads.iter().map(|r| r.qual.clone()).collect();

    let mut strand_hits = Vec::new();
    for cfg in configs {
        if !cli.quiet {
            eprintln!(
                "Aligning to {} strand ({:?}/{:?})...",
                cfg.name, cfg.read_conv, cfg.genome_conv
            );
        }
        let hits = align_strand_se(bowtie2, cfg, &reads, extra_args, cli.threads)?;
        strand_hits.push(hits);
    }

    let aligned = select_best_se(&strand_hits, configs, &reads_seq, &reads_qual, genome);

    let stem = output_stem(cli, input);
    let bam_name = format!("{}_bismark_bt2.bam", stem);
    let bam_path = cli.output_dir.join(&bam_name);

    let rg = rg_info(cli, input);
    let rg_ref = rg.as_ref().map(|(id, s)| (id.as_str(), s.as_str()));
    let header = build_sam_header(genome, extra_args, rg_ref);

    let mut stats = AlignStats::default();
    write_alignments_se(
        &aligned,
        &header,
        &bam_path,
        samtools,
        cli.sam_output,
        rg_ref.map(|(id, _)| id),
        &mut stats,
    )?;

    if cli.unmapped {
        write_unmapped_se(&reads, &aligned, &cli.output_dir, &stem)?;
    }
    if cli.ambiguous {
        write_ambiguous_se(&reads, &aligned, &cli.output_dir, &stem)?;
    }

    let bowtie2_cmd_str = format!("bowtie2 -p {}", cli.threads);
    let report_path = cli
        .output_dir
        .join(format!("{}_bismark_bt2_SE_report.txt", stem));
    write_se_report(&stats, genome_folder, input, &report_path, &bowtie2_cmd_str)?;

    if !cli.quiet {
        eprintln!("Done. Output: {}", bam_path.display());
        print_stats_se(&stats);
    }
    Ok(())
}

fn run_pe(
    cli: &Cli,
    genome_folder: &Path,
    genome: &Genome,
    configs: &[StrandConfig],
    bowtie2: &str,
    samtools: &str,
    extra_args: &[String],
) -> Result<()> {
    let r1 = cli.read1.as_deref().unwrap();
    let r2 = cli.read2.as_deref().unwrap();

    if !cli.quiet {
        eprintln!("Reading PE inputs: {} + {}", r1.display(), r2.display());
    }
    let mut pairs = if cli.fasta {
        fastq::read_fasta_pe(r1, r2)?
    } else {
        fastq::read_all_pe(r1, r2)?
    };

    // Apply --skip / --upto
    if let Some(s) = cli.skip {
        if s >= pairs.len() {
            pairs.clear();
        } else {
            pairs.drain(..s);
        }
    }
    if let Some(n) = cli.upto {
        pairs.truncate(n);
    }

    if !cli.quiet {
        eprintln!("{} read pairs loaded", pairs.len());
    }

    let r1_seqs: Vec<Vec<u8>> = pairs.iter().map(|(r, _)| r.seq.clone()).collect();
    let r2_seqs: Vec<Vec<u8>> = pairs.iter().map(|(_, r)| r.seq.clone()).collect();
    let r1_quals: Vec<Vec<u8>> = pairs.iter().map(|(r, _)| r.qual.clone()).collect();
    let r2_quals: Vec<Vec<u8>> = pairs.iter().map(|(_, r)| r.qual.clone()).collect();

    // Build PE-specific extra args (extend base args with pairing flags).
    let mut pe_args = extra_args.to_vec();
    pe_args.push("--no-mixed".to_owned());
    pe_args.push("--no-discordant".to_owned());
    if !cli.no_dovetail {
        pe_args.push("--dovetail".to_owned());
    }
    if cli.minins > 0 {
        pe_args.push("--minins".to_owned());
        pe_args.push(cli.minins.to_string());
    }
    pe_args.push("--maxins".to_owned());
    pe_args.push(cli.maxins.to_string());

    let mut strand_hits = Vec::new();
    for cfg in configs {
        if !cli.quiet {
            eprintln!(
                "Aligning to {} strand ({:?}/{:?})...",
                cfg.name, cfg.read_conv, cfg.genome_conv
            );
        }
        let hits = align_strand_pe(bowtie2, cfg, &pairs, &pe_args, cli.threads)?;
        strand_hits.push(hits);
    }

    let aligned = select_best_pe(
        &strand_hits,
        configs,
        &r1_seqs,
        &r2_seqs,
        &r1_quals,
        &r2_quals,
        genome,
    );

    let stem = output_stem(cli, r1);
    let bam_name = format!("{}_bismark_bt2_pe.bam", stem);
    let bam_path = cli.output_dir.join(&bam_name);

    let rg = rg_info(cli, r1);
    let rg_ref = rg.as_ref().map(|(id, s)| (id.as_str(), s.as_str()));
    let header = build_sam_header(genome, &pe_args, rg_ref);

    let mut stats = AlignStats::default();
    write_alignments_pe(
        &aligned,
        &header,
        &bam_path,
        samtools,
        cli.sam_output,
        rg_ref.map(|(id, _)| id),
        &mut stats,
    )?;

    if cli.unmapped {
        write_unmapped_pe(&pairs, &aligned, &cli.output_dir, &stem)?;
    }
    if cli.ambiguous {
        write_ambiguous_pe(&pairs, &aligned, &cli.output_dir, &stem)?;
    }

    let bowtie2_cmd_str = format!("bowtie2 -p {}", cli.threads);
    let report_path = cli
        .output_dir
        .join(format!("{}_bismark_bt2_PE_report.txt", stem));
    write_pe_report(
        &stats,
        genome_folder,
        r1,
        r2,
        &report_path,
        &bowtie2_cmd_str,
    )?;

    if !cli.quiet {
        eprintln!("Done. Output: {}", bam_path.display());
        print_stats_pe(&stats);
    }
    Ok(())
}

// ── SAM header ───────────────────────────────────────────────────────────────

fn build_sam_header(
    genome: &Genome,
    extra_args: &[String],
    rg: Option<(&str, &str)>,
) -> Vec<u8> {
    let mut h = Vec::new();
    h.extend_from_slice(b"@HD\tVN:1.0\tSO:unsorted\n");
    for (name, seq) in &genome.sequences {
        let line = format!("@SQ\tSN:{}\tLN:{}\n", name, seq.len());
        h.extend_from_slice(line.as_bytes());
    }
    if let Some((id, sample)) = rg {
        let rg_line = format!("@RG\tID:{}\tSM:{}\tPL:ILLUMINA\n", id, sample);
        h.extend_from_slice(rg_line.as_bytes());
    }
    let pg = format!(
        "@PG\tID:bistromark\tPN:bistromark\tVN:{}\tCL:bistromark {}\n",
        BISMARK_VERSION,
        extra_args.join(" ")
    );
    h.extend_from_slice(pg.as_bytes());
    h
}

// ── Writing alignments ───────────────────────────────────────────────────────

fn write_alignments_se(
    aligned: &[AlignedRead],
    header: &[u8],
    out_path: &Path,
    samtools: &str,
    sam_output: bool,
    rg_id: Option<&str>,
    stats: &mut AlignStats,
) -> Result<()> {
    if sam_output {
        let mut w = BufWriter::new(std::fs::File::create(out_path.with_extension("sam"))?);
        w.write_all(header)?;
        for a in aligned {
            stats.tally_outcome(a.outcome);
            if a.outcome == Outcome::Unique {
                stats.tally_xm(&a.xm);
                w.write_all(&format_sam_record(a, rg_id))?;
                w.write_all(b"\n")?;
            }
        }
    } else {
        let mut bw = BamWriter::open(samtools, out_path)?;
        for line in header.split(|&b| b == b'\n') {
            if !line.is_empty() {
                bw.write_line(line)?;
            }
        }
        for a in aligned {
            stats.tally_outcome(a.outcome);
            if a.outcome == Outcome::Unique {
                stats.tally_xm(&a.xm);
                bw.write_line(&format_sam_record(a, rg_id))?;
            }
        }
        bw.finish()?;
    }
    Ok(())
}

fn write_alignments_pe(
    aligned: &[(AlignedRead, AlignedRead)],
    header: &[u8],
    out_path: &Path,
    samtools: &str,
    sam_output: bool,
    rg_id: Option<&str>,
    stats: &mut AlignStats,
) -> Result<()> {
    if sam_output {
        let mut w = BufWriter::new(std::fs::File::create(out_path.with_extension("sam"))?);
        w.write_all(header)?;
        for (r1, r2) in aligned {
            stats.tally_outcome(r1.outcome);
            if r1.outcome == Outcome::Unique {
                stats.tally_xm(&r1.xm);
                stats.tally_xm(&r2.xm);
                w.write_all(&format_sam_record(r1, rg_id))?;
                w.write_all(b"\n")?;
                w.write_all(&format_sam_record(r2, rg_id))?;
                w.write_all(b"\n")?;
            }
        }
    } else {
        let mut bw = BamWriter::open(samtools, out_path)?;
        for line in header.split(|&b| b == b'\n') {
            if !line.is_empty() {
                bw.write_line(line)?;
            }
        }
        for (r1, r2) in aligned {
            stats.tally_outcome(r1.outcome);
            if r1.outcome == Outcome::Unique {
                stats.tally_xm(&r1.xm);
                stats.tally_xm(&r2.xm);
                bw.write_line(&format_sam_record(r1, rg_id))?;
                bw.write_line(&format_sam_record(r2, rg_id))?;
            }
        }
        bw.finish()?;
    }
    Ok(())
}

/// Serialize one `AlignedRead` to a SAM line (without trailing newline).
fn format_sam_record(a: &AlignedRead, rg_id: Option<&str>) -> Vec<u8> {
    let xr = match a.xr {
        bismark_lib::sam::ReadConversion::CT => "CT",
        bismark_lib::sam::ReadConversion::GA => "GA",
    };
    let xg = match a.xg {
        bismark_lib::sam::GenomeConversion::CT => "CT",
        bismark_lib::sam::GenomeConversion::GA => "GA",
    };

    let mut out = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        String::from_utf8_lossy(&a.qname),
        a.flag,
        String::from_utf8_lossy(&a.rname),
        a.pos,
        a.mapq,
        String::from_utf8_lossy(&a.cigar),
        String::from_utf8_lossy(&a.rnext),
        a.pnext,
        a.tlen,
        String::from_utf8_lossy(&a.seq),
        String::from_utf8_lossy(&a.qual),
    );

    if !a.extra_tags.is_empty() {
        out.push('\t');
        out.push_str(&String::from_utf8_lossy(&a.extra_tags));
    }
    if let Some(id) = rg_id {
        out.push_str("\tRG:Z:");
        out.push_str(id);
    }
    out.push_str(&format!(
        "\tXM:Z:{}\tXR:Z:{}\tXG:Z:{}",
        String::from_utf8_lossy(&a.xm),
        xr,
        xg
    ));
    out.into_bytes()
}

// ── Unmapped / ambiguous FASTQ output ────────────────────────────────────────

fn write_unmapped_se(
    reads: &[fastq::FastqRecord],
    aligned: &[AlignedRead],
    out_dir: &Path,
    stem: &str,
) -> Result<()> {
    let path = out_dir.join(format!("{stem}_unmapped_reads.fq"));
    let mut w = BufWriter::new(std::fs::File::create(&path)?);
    for (rec, a) in reads.iter().zip(aligned.iter()) {
        if a.outcome == Outcome::Unmapped {
            fastq::write_record(&mut w, rec)?;
        }
    }
    Ok(())
}

fn write_ambiguous_se(
    reads: &[fastq::FastqRecord],
    aligned: &[AlignedRead],
    out_dir: &Path,
    stem: &str,
) -> Result<()> {
    let path = out_dir.join(format!("{stem}_ambiguous_reads.fq"));
    let mut w = BufWriter::new(std::fs::File::create(&path)?);
    for (rec, a) in reads.iter().zip(aligned.iter()) {
        if a.outcome == Outcome::Ambiguous {
            fastq::write_record(&mut w, rec)?;
        }
    }
    Ok(())
}

fn write_unmapped_pe(
    pairs: &[(fastq::FastqRecord, fastq::FastqRecord)],
    aligned: &[(AlignedRead, AlignedRead)],
    out_dir: &Path,
    stem: &str,
) -> Result<()> {
    let path1 = out_dir.join(format!("{stem}_unmapped_reads_1.fq"));
    let path2 = out_dir.join(format!("{stem}_unmapped_reads_2.fq"));
    let mut w1 = BufWriter::new(std::fs::File::create(&path1)?);
    let mut w2 = BufWriter::new(std::fs::File::create(&path2)?);
    for ((r1, r2), (a1, _)) in pairs.iter().zip(aligned.iter()) {
        if a1.outcome == Outcome::Unmapped {
            fastq::write_record(&mut w1, r1)?;
            fastq::write_record(&mut w2, r2)?;
        }
    }
    Ok(())
}

fn write_ambiguous_pe(
    pairs: &[(fastq::FastqRecord, fastq::FastqRecord)],
    aligned: &[(AlignedRead, AlignedRead)],
    out_dir: &Path,
    stem: &str,
) -> Result<()> {
    let path1 = out_dir.join(format!("{stem}_ambiguous_reads_1.fq"));
    let path2 = out_dir.join(format!("{stem}_ambiguous_reads_2.fq"));
    let mut w1 = BufWriter::new(std::fs::File::create(&path1)?);
    let mut w2 = BufWriter::new(std::fs::File::create(&path2)?);
    for ((r1, r2), (a1, _)) in pairs.iter().zip(aligned.iter()) {
        if a1.outcome == Outcome::Ambiguous {
            fastq::write_record(&mut w1, r1)?;
            fastq::write_record(&mut w2, r2)?;
        }
    }
    Ok(())
}

// ── Utility helpers ───────────────────────────────────────────────────────────

/// Compute the output stem from CLI overrides or the input file name.
fn output_stem(cli: &Cli, r1_path: &Path) -> String {
    if let Some(b) = &cli.basename {
        return b.clone();
    }
    let auto = fastq_stem(r1_path);
    if let Some(p) = &cli.prefix {
        format!("{}_{}", p, auto)
    } else {
        auto
    }
}

/// Build RG tag info if --rg_tag is set.
fn rg_info(cli: &Cli, r1_path: &Path) -> Option<(String, String)> {
    if !cli.rg_tag {
        return None;
    }
    let id = cli.rg_id.clone().unwrap_or_else(|| "1".to_owned());
    let sample = cli
        .rg_sample
        .clone()
        .unwrap_or_else(|| fastq_stem(r1_path));
    Some((id, sample))
}

/// Strip FASTQ/FASTA extensions from a filename to get the output stem.
fn fastq_stem(path: &Path) -> String {
    let name = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string();
    for ext in &[
        ".fastq.gz",
        ".fastq",
        ".fq.gz",
        ".fq",
        ".fasta.gz",
        ".fasta",
        ".fa.gz",
        ".fa",
    ] {
        if let Some(s) = name.strip_suffix(ext) {
            return s.to_owned();
        }
    }
    name
}

/// Locate the bowtie2 binary.
fn find_bowtie2(path_hint: Option<&Path>) -> Result<String> {
    if let Some(hint) = path_hint {
        if hint.is_dir() {
            let candidate = hint.join("bowtie2");
            if candidate.exists() {
                return Ok(candidate.to_string_lossy().into_owned());
            }
        } else if hint.exists() {
            return Ok(hint.to_string_lossy().into_owned());
        }
        bail!("bowtie2 not found at {}", hint.display());
    }
    if which("bowtie2").is_some() {
        return Ok("bowtie2".to_owned());
    }
    bail!("bowtie2 not found in PATH; use --path_to_bowtie2")
}

fn which(name: &str) -> Option<PathBuf> {
    std::env::var_os("PATH").and_then(|paths| {
        std::env::split_paths(&paths).find_map(|dir| {
            let full = dir.join(name);
            if full.is_file() {
                Some(full)
            } else {
                None
            }
        })
    })
}

fn print_stats_se(s: &AlignStats) {
    eprintln!("Sequences analysed: {}", s.total_reads);
    eprintln!(
        "Unique best-hit alignments: {} ({:.1}%)",
        s.unique,
        pct(s.unique, s.total_reads)
    );
    eprintln!("Unmapped: {}", s.unmapped);
    eprintln!("Ambiguous: {}", s.ambiguous);
}

fn print_stats_pe(s: &AlignStats) {
    eprintln!("Sequence pairs analysed: {}", s.total_reads);
    eprintln!(
        "Unique best-hit PE alignments: {} ({:.1}%)",
        s.unique,
        pct(s.unique, s.total_reads)
    );
    eprintln!("Unmapped pairs: {}", s.unmapped);
    eprintln!("Ambiguous pairs: {}", s.ambiguous);
}

fn pct(n: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        n as f64 / total as f64 * 100.0
    }
}
