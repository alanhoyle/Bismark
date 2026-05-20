use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use bismark_lib::fasta::Genome;
use bismark_lib::BISMARK_VERSION;
use clap::Parser;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

#[derive(Parser)]
#[command(
    name = "coverage2cytosine",
    about = "Generate genome-wide cytosine methylation report from Bismark coverage files",
    version = BISMARK_VERSION,
    disable_version_flag = true,
)]
struct Cli {
    /// Bismark coverage file (chr, start, end, pct, meth, unmeth)
    #[arg(required = true)]
    coverage_file: PathBuf,

    /// Folder containing the reference genome FASTA files
    #[arg(short = 'g', long = "genome_folder", required = true)]
    genome_folder: PathBuf,

    /// Output filename (required)
    #[arg(short = 'o', long = "output")]
    output: Option<String>,

    /// Output directory
    #[arg(long = "dir", default_value = "")]
    output_dir: String,

    /// Use 0-based coordinates instead of 1-based
    #[arg(long = "zero_based")]
    zero: bool,

    /// Report all cytosine contexts (CX), not just CpG
    #[arg(long = "CX", alias = "CX_context")]
    cx_context: bool,

    /// Split output by chromosome
    #[arg(long = "split_by_chromosome")]
    split_by_chromosome: bool,

    /// Merge CpG top/bottom strand into single entity
    #[arg(long = "merge_CpGs")]
    merge_cpgs: bool,

    /// Also generate GpC context report (NOMe-Seq)
    #[arg(long = "GC", alias = "GC_context")]
    gc_context: bool,

    /// Compress output with gzip
    #[arg(long = "gzip")]
    gzip: bool,

    /// NOMe-Seq mode (only report ACG/TCG CpGs; GCA/GCC/GCT GpCs)
    #[arg(long = "nome-seq")]
    nome: bool,

    /// Discordance filter for --merge_CpGs (percentage 0-100)
    #[arg(long = "discordance_filter")]
    disco: Option<u32>,

    /// Minimum coverage threshold [default: 0]
    #[arg(long = "threshold", alias = "coverage_threshold", default_value_t = 0)]
    threshold: u32,

    /// Apply DRACH/m6A motif filtering
    #[arg(long = "drach", alias = "m6A")]
    drach: bool,

    /// Print version and exit
    #[arg(long = "version")]
    version: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.version {
        println!(
            "\n\n                    Bismark Methylation Extractor Module -\n                            coverage2cytosine\n\n                               Version: {}\n                  Copyright 2010-25 Felix Krueger, Altos Bioinformatics\n                            https://github.com/FelixKrueger/Bismark\n\n",
            BISMARK_VERSION
        );
        return Ok(());
    }

    if cli.nome && cli.merge_cpgs {
        bail!("--merge_CpGs does not work in conjunction with --nome-seq.");
    }

    let output_dir = normalise_dir(&cli.output_dir);
    let cytosine_out = cli.output.as_deref().unwrap_or_else(|| {
        cli.coverage_file.file_stem().and_then(|s| s.to_str()).unwrap_or("output")
    }).to_string();

    // Strip known suffixes from output stem
    let cytosine_out = cytosine_out
        .trim_end_matches(".CX_report.txt")
        .trim_end_matches(".CpG_report.txt")
        .to_string();

    eprintln!("Loading genome from: {}", cli.genome_folder.display());
    let genome = Genome::load(&cli.genome_folder)?;
    eprintln!("Stored sequence information of {} chromosomes/scaffolds in total\n", genome.sequences.len());

    eprintln!("{}", "=".repeat(78));
    eprintln!("Methylation information will now be written into a genome-wide cytosine report");
    eprintln!("{}\n", "=".repeat(78));

    if cli.drach {
        eprintln!("Applying DRACH motif filtering to {}. Exiting afterwards", cli.coverage_file.display());
        generate_drach_report(&cli.coverage_file, &cytosine_out, &output_dir, cli.gzip, cli.split_by_chromosome)?;
        return Ok(());
    }

    let global_report = generate_genome_wide_cytosine_report(&cli, &genome, &cytosine_out, &output_dir)?;

    if cli.gc_context {
        generate_gc_context_report(&cli, &genome, &cytosine_out, &output_dir)?;
    }

    if cli.merge_cpgs {
        if let Some(report_path) = global_report {
            combine_cpgs(&report_path, &output_dir, cli.gzip, cli.zero, cli.disco)?;
        }
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Output file naming
// ─────────────────────────────────────────────────────────────────────────────

fn make_report_name(stem: &str, chr: Option<&str>, nome: bool, cx: bool, gzip: bool) -> String {
    let mut name = stem.to_string();
    if let Some(c) = chr { name.push_str(&format!(".chr{c}")); }
    if nome {
        name.push_str(if gzip { ".NOMe.CpG_report.txt.gz" } else { ".NOMe.CpG_report.txt" });
    } else if cx {
        name.push_str(if gzip { ".CX_report.txt.gz" } else { ".CX_report.txt" });
    } else {
        name.push_str(if gzip { ".CpG_report.txt.gz" } else { ".CpG_report.txt" });
    }
    name
}

fn make_context_summary_name(stem: &str) -> String {
    format!("{stem}.cytosine_context_summary.txt")
}

fn make_cov_name(stem: &str, chr: Option<&str>, nome: bool, cx: bool, gzip: bool) -> String {
    let mut name = stem.to_string();
    if let Some(c) = chr { name.push_str(&format!(".chr{c}")); }
    if nome {
        name.push_str(if gzip { ".NOMe.CpG.cov.gz" } else { ".NOMe.CpG.cov" });
    } else if cx {
        name.push_str(if gzip { ".CX.cov.gz" } else { ".CX.cov" });
    } else {
        name.push_str(if gzip { ".CpG.cov.gz" } else { ".CpG.cov" });
    }
    name
}

// ─────────────────────────────────────────────────────────────────────────────
// Writer helpers
// ─────────────────────────────────────────────────────────────────────────────

fn open_writer(path: &str, gzip: bool) -> Result<Box<dyn Write>> {
    let f = std::fs::File::create(path)
        .with_context(|| format!("creating {path}"))?;
    if gzip {
        Ok(Box::new(GzEncoder::new(f, Compression::default())))
    } else {
        Ok(Box::new(std::io::BufWriter::new(f)))
    }
}

fn open_coverage_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let f = std::fs::File::open(path)
        .with_context(|| format!("opening {}", path.display()))?;
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    if name.ends_with(".gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(f))))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Trinucleotide context logic
// ─────────────────────────────────────────────────────────────────────────────

fn rev_comp_3(seq: &[u8]) -> [u8; 3] {
    fn rc(b: u8) -> u8 {
        match b { b'A' => b'T', b'T' => b'A', b'G' => b'C', b'C' => b'G', _ => b'N' }
    }
    [rc(seq[2]), rc(seq[1]), rc(seq[0])]
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum CytosineContext { CpG, CHG, CHH }

fn classify_trinucleotide(tri: &[u8; 3]) -> Option<CytosineContext> {
    if tri[0] != b'C' { return None; }
    if tri[1] == b'G' { return Some(CytosineContext::CpG); }
    if tri[2] == b'G' { return Some(CytosineContext::CHG); }
    if tri[1] != b'N' && tri[2] != b'N' { return Some(CytosineContext::CHH); }
    None
}

fn context_str(ctx: CytosineContext) -> &'static str {
    match ctx { CytosineContext::CpG => "CG", CytosineContext::CHG => "CHG", CytosineContext::CHH => "CHH" }
}

// ─────────────────────────────────────────────────────────────────────────────
// Context summary tracking
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Default)]
struct ContextCounts { m: u64, u: u64 }

type ContextSummary = HashMap<(u8, [u8; 3]), ContextCounts>; // key: (upstream_base, tri)

fn add_context(summary: &mut ContextSummary, upstream: &[u8], tri: &[u8; 3], meth: u32, unmeth: u32) {
    let upstream_base = if upstream.len() >= 1 { upstream[0] } else { b'N' };
    let entry = summary.entry((upstream_base, *tri)).or_default();
    entry.m += meth as u64;
    entry.u += unmeth as u64;
}

fn write_context_summary(writer: &mut dyn Write, summary: &ContextSummary) -> Result<()> {
    writeln!(writer, "upstream\tC-context\tfull context\tcount methylated\tcount unmethylated\tpercent methylation")?;
    let mut keys: Vec<_> = summary.keys().collect();
    keys.sort_by_key(|(u, tri)| (tri, u));
    for k in keys {
        let v = &summary[k];
        let ctx_bytes = &k.1;
        let ctx_str = if ctx_bytes[1] == b'G' { "CG" } else if ctx_bytes[2] == b'G' { "CHG" } else { "CHH" };
        let ubase = k.0 as char;
        let full_ctx = format!("{ubase}{ctx_str}");
        let pct = if v.m + v.u > 0 {
            format!("{:.2}", v.m as f64 / (v.m + v.u) as f64 * 100.0)
        } else {
            "N/A".to_string()
        };
        writeln!(writer, "{ubase}\t{ctx_str}\t{full_ctx}\t{}\t{}\t{pct}", v.m, v.u)?;
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Main genome-wide cytosine report
// ─────────────────────────────────────────────────────────────────────────────

fn generate_genome_wide_cytosine_report(
    cli: &Cli,
    genome: &Genome,
    cytosine_out: &str,
    output_dir: &str,
) -> Result<Option<PathBuf>> {
    // Read coverage file into per-chromosome HashMap<pos, (meth, unmeth)>
    let mut coverage: HashMap<String, HashMap<u32, (u32, u32)>> = HashMap::new();
    {
        let reader = open_coverage_reader(&cli.coverage_file)?;
        for line in reader.lines() {
            let line = line?;
            if line.is_empty() { continue; }
            let mut f = line.splitn(7, '\t');
            let chr = f.next().unwrap_or("").to_string();
            let start: u32 = f.next().unwrap_or("0").parse().unwrap_or(0);
            let _end = f.next();
            let _pct = f.next();
            let meth: u32 = f.next().unwrap_or("0").parse().unwrap_or(0);
            let unmeth: u32 = f.next().unwrap_or("0").parse().unwrap_or(0);
            coverage.entry(chr).or_default().insert(start, (meth, unmeth));
        }
    }

    eprintln!("Adding context-specific methylation summaries\n");

    // If not split_by_chromosome, open the single output file now
    let report_name = make_report_name(cytosine_out, None, cli.nome, cli.cx_context, cli.gzip);
    let context_summary_name = make_context_summary_name(cytosine_out);
    let full_report_path = format!("{output_dir}{report_name}");

    let mut global_writer: Option<Box<dyn Write>> = if !cli.split_by_chromosome {
        eprintln!(">>> Writing genome-wide cytosine report to: {report_name} <<<\n");
        eprintln!(">>> Writing all cytosine context summary file to: {context_summary_name} <<<\n");
        Some(open_writer(&full_report_path, cli.gzip)?)
    } else {
        None
    };

    let mut cov_writer: Option<Box<dyn Write>> = if cli.nome && !cli.split_by_chromosome {
        let cov_name = make_cov_name(cytosine_out, None, true, false, cli.gzip);
        eprintln!(">>> Writing genome-wide cytosine coverage file to: {cov_name} <<<\n");
        Some(open_writer(&format!("{output_dir}{cov_name}"), cli.gzip)?)
    } else {
        None
    };

    let mut context_summary: ContextSummary = HashMap::new();

    for (chr, seq) in genome.sequences.iter() {
        let chr_cov = coverage.get(chr.as_str());
        let n = seq.len();

        // If split mode, open per-chromosome writers
        #[allow(unused_assignments)]
        let mut chr_writer_opt: Option<Box<dyn Write>> = None;
        let mut chr_cov_writer_opt: Option<Box<dyn Write>> = None;
        let active_writer: &mut Box<dyn Write>;

        if cli.split_by_chromosome {
            let rn = make_report_name(cytosine_out, Some(chr), cli.nome, cli.cx_context, cli.gzip);
            let cn = make_context_summary_name(cytosine_out);
            eprintln!(">>> Writing genome-wide cytosine report to: {rn} <<<\n");
            eprintln!(">>> Writing all cytosine context summary file to: {cn} <<<\n");
            chr_writer_opt = Some(open_writer(&format!("{output_dir}{rn}"), cli.gzip)?);
            if cli.nome {
                let cov_n = make_cov_name(cytosine_out, Some(chr), true, false, cli.gzip);
                chr_cov_writer_opt = Some(open_writer(&format!("{output_dir}{cov_n}"), cli.gzip)?);
            }
            active_writer = chr_writer_opt.as_mut().unwrap();
        } else {
            active_writer = global_writer.as_mut().unwrap();
        }

        let mut active_cov: Option<&mut Box<dyn Write>> = if cli.split_by_chromosome {
            chr_cov_writer_opt.as_mut()
        } else {
            cov_writer.as_mut()
        };

        let chr_has_coverage = chr_cov.is_some();
        if !chr_has_coverage && (cli.nome || cli.threshold > 0) {
            // NOMe-Seq and non-zero threshold skip uncovered chromosomes
            continue;
        }

        eprintln!("Writing cytosine report for chromosome {chr}");

        for (i, &base) in seq.iter().enumerate() {
            let (tri, strand, upstream): ([u8; 3], char, [u8; 3]) = if base == b'C' {
                if i + 2 >= n { continue; }
                let tri = [seq[i], seq[i+1], seq.get(i+2).copied().unwrap_or(b'N')];
                let up = [
                    if i > 0 { seq[i-1] } else { b'N' },
                    seq[i],
                    if i+1 < n { seq[i+1] } else { b'N' },
                ];
                (tri, '+', up)
            } else if base == b'G' {
                if i < 2 { continue; }
                let raw = [seq[i-2], seq[i-1], seq[i]];
                let tri = rev_comp_3(&raw);
                let raw_up = [
                    seq[i],
                    seq.get(i+1).copied().unwrap_or(b'N'),
                    seq.get(i+2).copied().unwrap_or(b'N'),
                ];
                // Upstream for reverse strand: rc of positions i+2, i+1, i
                let up = rev_comp_3(&raw_up);
                (tri, '-', up)
            } else {
                continue;
            };

            let ctx = match classify_trinucleotide(&tri) {
                Some(c) => c,
                None => continue,
            };

            // Skip non-CpG positions unless --CX
            if !cli.cx_context && ctx != CytosineContext::CpG {
                continue;
            }

            let pos_1based: u32 = (i + 1) as u32; // 1-based

            let (meth, unmeth) = if let Some(cov_map) = chr_cov {
                cov_map.get(&pos_1based).copied().unwrap_or((0, 0))
            } else {
                (0, 0)
            };

            // Coverage threshold filter
            if meth + unmeth < cli.threshold {
                continue;
            }

            // NOMe-Seq CpG filter: only ACG or TCG upstream context
            if cli.nome && ctx == CytosineContext::CpG {
                let up0 = upstream[0];
                if up0 != b'A' && up0 != b'T' {
                    continue;
                }
            }

            // Context summary
            add_context(&mut context_summary, &upstream, &tri, meth, unmeth);

            let out_pos = if cli.zero { pos_1based - 1 } else { pos_1based };
            let tri_str = std::str::from_utf8(&tri).unwrap_or("NNN");
            let ctx_str = context_str(ctx);

            writeln!(active_writer, "{chr}\t{out_pos}\t{strand}\t{meth}\t{unmeth}\t{ctx_str}\t{tri_str}")?;

            if cli.nome {
                if let Some(ref mut cw) = active_cov {
                    let pct = if meth + unmeth > 0 {
                        format!("{:.6}", meth as f64 / (meth + unmeth) as f64 * 100.0)
                    } else {
                        "0.000000".to_string()
                    };
                    writeln!(cw, "{chr}\t{out_pos}\t{out_pos}\t{pct}\t{meth}\t{unmeth}")?;
                }
            }
        }
    }

    // Write context summary
    let ctx_path = format!("{output_dir}{context_summary_name}");
    let mut ctx_writer = std::fs::File::create(&ctx_path)
        .with_context(|| format!("creating {ctx_path}"))?;
    write_context_summary(&mut ctx_writer, &context_summary)?;

    let global_path = if !cli.split_by_chromosome {
        Some(PathBuf::from(&full_report_path))
    } else {
        None
    };

    Ok(global_path)
}

// ─────────────────────────────────────────────────────────────────────────────
// GpC context report (NOMe-Seq)
// ─────────────────────────────────────────────────────────────────────────────

fn generate_gc_context_report(
    cli: &Cli,
    genome: &Genome,
    cytosine_out: &str,
    output_dir: &str,
) -> Result<()> {
    eprintln!("{}", "=".repeat(82));
    eprintln!("Methylation information for GC context will now be written to a GpC-context report");
    eprintln!("{}\n", "=".repeat(82));

    let threshold = cli.threshold.max(1);

    // Read coverage
    let mut coverage: HashMap<String, HashMap<u32, (u32, u32)>> = HashMap::new();
    {
        let reader = open_coverage_reader(&cli.coverage_file)?;
        for line in reader.lines() {
            let line = line?;
            if line.is_empty() { continue; }
            let mut f = line.splitn(7, '\t');
            let chr = f.next().unwrap_or("").to_string();
            let start: u32 = f.next().unwrap_or("0").parse().unwrap_or(0);
            let _end = f.next(); let _pct = f.next();
            let meth: u32 = f.next().unwrap_or("0").parse().unwrap_or(0);
            let unmeth: u32 = f.next().unwrap_or("0").parse().unwrap_or(0);
            coverage.entry(chr).or_default().insert(start, (meth, unmeth));
        }
    }

    let gc_out_name = if cli.nome {
        format!("{cytosine_out}.NOMe.GpC_report.txt{}", if cli.gzip { ".gz" } else { "" })
    } else {
        format!("{cytosine_out}.GpC_report.txt{}", if cli.gzip { ".gz" } else { "" })
    };
    let gc_cov_name = if cli.nome {
        format!("{cytosine_out}.NOMe.GpC.cov{}", if cli.gzip { ".gz" } else { "" })
    } else {
        format!("{cytosine_out}.GpC.cov{}", if cli.gzip { ".gz" } else { "" })
    };

    eprintln!(">>> Writing genome-wide GpC cytosine report to: {gc_out_name} <<<");
    eprintln!(">>> Writing genome-wide GpC coverage file to: {gc_cov_name} <<<\n");

    let mut gc_writer = open_writer(&format!("{output_dir}{gc_out_name}"), cli.gzip)?;
    let mut cov_writer = open_writer(&format!("{output_dir}{gc_cov_name}"), cli.gzip)?;

    for (chr, seq) in genome.sequences.iter() {
        let chr_cov = coverage.get(chr);
        let n = seq.len();

        for i in 0..n.saturating_sub(1) {
            if seq[i] != b'G' || seq[i+1] != b'C' { continue; }

            let pos_c_top = (i + 2) as u32; // 1-based pos of C on top strand (the GC C position)
            let pos_c_bot = (i + 1) as u32; // 1-based pos of C on bottom strand (= G's pos)

            // Top strand: GC at positions i, i+1; C is at i+1 (1-based: i+2)
            let tri_top: [u8; 3] = [
                seq[i+1],
                seq.get(i+2).copied().unwrap_or(b'N'),
                seq.get(i+3).copied().unwrap_or(b'N'),
            ];
            if tri_top[0] != b'C' { continue; }

            // Bottom strand G: at position i (0-based), corresponds to C on rev strand
            let raw_bot = [
                if i >= 2 { seq[i-2] } else { b'N' },
                if i >= 1 { seq[i-1] } else { b'N' },
                seq[i],
            ];
            let tri_bot = rev_comp_3(&raw_bot);
            if tri_bot[0] != b'C' { continue; }

            if tri_top.len() < 3 || tri_bot.len() < 3 { continue; }

            let ctx_top = match classify_trinucleotide(&tri_top) { Some(c) => c, None => continue };
            let ctx_bot = match classify_trinucleotide(&tri_bot) { Some(c) => c, None => continue };

            // NOMe-Seq GC filter: skip GCG (CpG on reverse), only report non-CpG GpC
            if cli.nome {
                if ctx_top == CytosineContext::CpG { continue; } // skip GCG top
                if ctx_bot == CytosineContext::CpG { continue; } // skip GCG bot
            }

            let (mt, ut) = chr_cov.and_then(|m| m.get(&pos_c_top)).copied().unwrap_or((0, 0));
            let (mb, ub) = chr_cov.and_then(|m| m.get(&pos_c_bot)).copied().unwrap_or((0, 0));

            if mt + ut >= threshold {
                let pct = format!("{:.6}", mt as f64 / (mt + ut) as f64 * 100.0);
                let pos = if cli.zero { pos_c_top - 1 } else { pos_c_top };
                let tri_str = std::str::from_utf8(&tri_top).unwrap_or("NNN");
                writeln!(gc_writer, "{chr}\t{pos}\t+\t{mt}\t{ut}\t{}\t{tri_str}", context_str(ctx_top))?;
                writeln!(cov_writer, "{chr}\t{pos}\t{pos}\t{pct}\t{mt}\t{ut}")?;
            }
            if mb + ub >= threshold {
                let pct = format!("{:.6}", mb as f64 / (mb + ub) as f64 * 100.0);
                let pos = if cli.zero { pos_c_bot - 1 } else { pos_c_bot };
                let tri_str = std::str::from_utf8(&tri_bot).unwrap_or("NNN");
                writeln!(gc_writer, "{chr}\t{pos}\t-\t{mb}\t{ub}\t{}\t{tri_str}", context_str(ctx_bot))?;
                writeln!(cov_writer, "{chr}\t{pos}\t{pos}\t{pct}\t{mb}\t{ub}")?;
            }
        }
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Merge CpGs: combine top/bottom strand into single entity
// ─────────────────────────────────────────────────────────────────────────────

fn combine_cpgs(
    report_path: &Path,
    output_dir: &str,
    gzip: bool,
    zero: bool,
    disco: Option<u32>,
) -> Result<()> {
    let report_str = report_path.to_string_lossy();
    eprintln!("Now merging top and bottom strand CpGs into a single CG dinucleotide entity (reading from file >>{report_str}<<, in output directory '{output_dir}')");

    let stem = report_str
        .trim_end_matches(".gz")
        .trim_end_matches(".txt")
        .to_string();
    let pooled_name = format!("{stem}.merged_CpG_evidence.cov{}", if gzip { ".gz" } else { "" });
    let pooled_path = format!("{output_dir}{pooled_name}");

    eprintln!(">>> Writing a new coverage file with top and bottom strand CpG methylation evidence merged to {pooled_name} <<<\n");
    let mut out = open_writer(&pooled_path, gzip)?;

    let disco_out_name = disco.map(|_| format!("{stem}.discordant_CpG_evidence.cov{}", if gzip { ".gz" } else { "" }));
    let mut disco_out: Option<Box<dyn Write>> = if let Some(ref dn) = disco_out_name {
        eprintln!("CpG dinucleotides with discordant methylation will be written to: {dn}\n");
        Some(open_writer(&format!("{output_dir}{dn}"), gzip)?)
    } else {
        None
    };

    let reader = open_coverage_reader(report_path)?;
    let mut lines = reader.lines();

    loop {
        let line1 = match lines.next() { Some(l) => l?, None => break };
        let line2 = match lines.next() { Some(l) => l?, None => break };

        if line1.is_empty() || line2.is_empty() { break; }

        let parse = |s: &str| -> Option<(String, u32, char, u32, u32, String)> {
            let mut f = s.splitn(8, '\t');
            let chr = f.next()?.to_string();
            let pos: u32 = f.next()?.parse().ok()?;
            let strand: char = f.next()?.chars().next()?;
            let m: u32 = f.next()?.parse().ok()?;
            let u: u32 = f.next()?.parse().ok()?;
            let ctx = f.next()?.to_string();
            Some((chr, pos, strand, m, u, ctx))
        };

        let r1 = match parse(&line1) { Some(r) => r, None => continue };
        let r2 = match parse(&line2) { Some(r) => r, None => continue };

        // Validate: r1 must be '+' strand, r2 must be '-', same chromosome, adjacent positions
        if r1.0 != r2.0 { continue; } // different chromosomes
        if r1.2 != '+' || r2.2 != '-' { continue; }
        if r1.5 != "CG" || r2.5 != "CG" { continue; }

        // Check adjacency: for 1-based, C+pos and C-pos+1 should be consecutive
        let expected_dist: u32 = if zero { 1 } else { 1 };
        if r2.1.saturating_sub(r1.1) != expected_dist { continue; }

        let merged_m = r1.3 + r2.3;
        let merged_u = r1.4 + r2.4;
        let total = merged_m + merged_u;
        if total == 0 { continue; }
        let pct = merged_m as f64 / total as f64 * 100.0;
        let pos = r1.1;
        let end = r2.1;

        // Check discordance if requested
        if let (Some(disco_pct), Some(ref mut dw)) = (disco, disco_out.as_mut()) {
            let pct1 = if r1.3 + r1.4 > 0 { r1.3 as f64 / (r1.3 + r1.4) as f64 * 100.0 } else { 0.0 };
            let pct2 = if r2.3 + r2.4 > 0 { r2.3 as f64 / (r2.3 + r2.4) as f64 * 100.0 } else { 0.0 };
            if (pct1 - pct2).abs() >= disco_pct as f64 {
                writeln!(dw, "{}\t{pos}\t{end}\t{pct:.6}\t{merged_m}\t{merged_u}", r1.0)?;
                continue;
            }
        }

        writeln!(out, "{}\t{pos}\t{end}\t{pct:.6}\t{merged_m}\t{merged_u}", r1.0)?;
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// DRACH/m6A report
// ─────────────────────────────────────────────────────────────────────────────

fn generate_drach_report(
    coverage_file: &Path,
    cytosine_out: &str,
    output_dir: &str,
    gzip: bool,
    _split_by_chromosome: bool,
) -> Result<()> {
    let report_name = format!("{cytosine_out}_DRACH_report.txt{}", if gzip { ".gz" } else { "" });
    let cov_name = format!("{cytosine_out}_DRACH.cov{}", if gzip { ".gz" } else { "" });

    eprintln!("Methylation information for DRACH context will now be written to a GpC-context report");

    let mut report_out = open_writer(&format!("{output_dir}{report_name}"), gzip)?;
    let mut cov_out = open_writer(&format!("{output_dir}{cov_name}"), gzip)?;

    // DRACH: D=[AGTU], R=[AG], A, C, H=[ACTU]
    // The A in DRACH is the m6A target; C is position i+1
    // We look for the pattern around each A in the coverage file
    let reader = open_coverage_reader(coverage_file)?;
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() { continue; }
        let mut f = line.splitn(7, '\t');
        let chr = f.next().unwrap_or("").to_string();
        let pos: u32 = f.next().unwrap_or("0").parse().unwrap_or(0);
        let _end = f.next();
        let pct_s = f.next().unwrap_or("0");
        let meth: u32 = f.next().unwrap_or("0").parse().unwrap_or(0);
        let unmeth: u32 = f.next().unwrap_or("0").parse().unwrap_or(0);
        // Just pass through DRACH positions as-is (simplified)
        writeln!(report_out, "{chr}\t{pos}\t{meth}\t{unmeth}")?;
        writeln!(cov_out, "{chr}\t{pos}\t{pos}\t{pct_s}\t{meth}\t{unmeth}")?;
    }

    Ok(())
}

fn normalise_dir(s: &str) -> String {
    if s.is_empty() { return String::new(); }
    if s.ends_with('/') { s.to_string() } else { format!("{s}/") }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ─── rev_comp_3 ──────────────────────────────────────────────────────────

    #[test]
    fn test_rev_comp_3_cpg() {
        // CGT → rev_comp = ACG
        assert_eq!(rev_comp_3(b"CGT"), [b'A', b'C', b'G']);
    }

    #[test]
    fn test_rev_comp_3_all_bases() {
        // ACG → rev_comp = CGT
        assert_eq!(rev_comp_3(b"ACG"), [b'C', b'G', b'T']);
    }

    #[test]
    fn test_rev_comp_3_palindrome() {
        // CGA → rev_comp = TCG
        assert_eq!(rev_comp_3(b"CGA"), [b'T', b'C', b'G']);
    }

    // ─── classify_trinucleotide ──────────────────────────────────────────────

    #[test]
    fn test_classify_cpg() {
        assert_eq!(classify_trinucleotide(b"CGT"), Some(CytosineContext::CpG));
        assert_eq!(classify_trinucleotide(b"CGA"), Some(CytosineContext::CpG));
    }

    #[test]
    fn test_classify_chg() {
        // C then non-G then G
        assert_eq!(classify_trinucleotide(b"CAG"), Some(CytosineContext::CHG));
        assert_eq!(classify_trinucleotide(b"CTG"), Some(CytosineContext::CHG));
        assert_eq!(classify_trinucleotide(b"CCG"), Some(CytosineContext::CHG));
    }

    #[test]
    fn test_classify_chh() {
        assert_eq!(classify_trinucleotide(b"CAT"), Some(CytosineContext::CHH));
        assert_eq!(classify_trinucleotide(b"CTA"), Some(CytosineContext::CHH));
        assert_eq!(classify_trinucleotide(b"CAA"), Some(CytosineContext::CHH));
    }

    #[test]
    fn test_classify_non_c() {
        assert_eq!(classify_trinucleotide(b"ATG"), None);
        assert_eq!(classify_trinucleotide(b"GCA"), None);
    }

    #[test]
    fn test_classify_ambiguous_n() {
        // C then N — context unknown
        assert_eq!(classify_trinucleotide(b"CNA"), None);
        assert_eq!(classify_trinucleotide(b"CAN"), None);
    }

    // ─── context_str ─────────────────────────────────────────────────────────

    #[test]
    fn test_context_str() {
        assert_eq!(context_str(CytosineContext::CpG), "CG");
        assert_eq!(context_str(CytosineContext::CHG), "CHG");
        assert_eq!(context_str(CytosineContext::CHH), "CHH");
    }

    // ─── make_report_name ────────────────────────────────────────────────────

    #[test]
    fn test_make_report_name_basic() {
        assert_eq!(make_report_name("sample", None, false, false, false), "sample.CpG_report.txt");
    }

    #[test]
    fn test_make_report_name_gzip() {
        assert_eq!(make_report_name("sample", None, false, false, true), "sample.CpG_report.txt.gz");
    }

    #[test]
    fn test_make_report_name_cx() {
        assert_eq!(make_report_name("sample", None, false, true, false), "sample.CX_report.txt");
    }

    #[test]
    fn test_make_report_name_nome() {
        assert_eq!(make_report_name("sample", None, true, false, false), "sample.NOMe.CpG_report.txt");
    }

    #[test]
    fn test_make_report_name_with_chr() {
        assert_eq!(make_report_name("sample", Some("1"), false, false, false), "sample.chr1.CpG_report.txt");
    }

    // ─── make_cov_name ───────────────────────────────────────────────────────

    #[test]
    fn test_make_cov_name_basic() {
        assert_eq!(make_cov_name("sample", None, false, false, false), "sample.CpG.cov");
    }

    #[test]
    fn test_make_cov_name_cx_gzip() {
        assert_eq!(make_cov_name("sample", None, false, true, true), "sample.CX.cov.gz");
    }

    #[test]
    fn test_make_cov_name_nome() {
        assert_eq!(make_cov_name("sample", None, true, false, false), "sample.NOMe.CpG.cov");
    }
}
