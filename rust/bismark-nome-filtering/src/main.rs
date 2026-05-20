/// NOMe_filtering: per-read CpG and GC methylation counting for NOMe-seq.
///
/// Input:  YACHT-format methylation file (tab-delimited):
///         ReadID  state(+/-)  chr  pos  context(Z/z/X/x/H/h)  start  end  strand
/// Output: Gzipped tab-delimited:
///         ReadID  Chr  Start  End  meth_CG  unmeth_CG  meth_GC  unmeth_GC
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use bismark_lib::fasta::Genome;
use bismark_lib::BISMARK_VERSION;
use clap::Parser;
use flate2::write::GzEncoder;
use flate2::Compression;

#[derive(Parser)]
#[command(
    name = "NOMe_filtering",
    about = "Per-read NOMe-seq CpG/GC methylation filtering (requires YACHT input)",
    version = BISMARK_VERSION,
    disable_version_flag = true,
)]
struct Cli {
    /// YACHT-format methylation file to process
    #[arg(required = true)]
    infile: PathBuf,

    /// Folder containing the reference genome FASTA files
    #[arg(short = 'g', long = "genome_folder", required = true)]
    genome_folder: PathBuf,

    /// Output directory
    #[arg(long = "dir", default_value = "")]
    output_dir: String,

    /// Print version and exit
    #[arg(long = "version")]
    version: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.version {
        println!("\n\tBismark NOMe Filtering version: {BISMARK_VERSION}\n\tCopyright 2010-25 Felix Krueger, Altos Bioinformatics\n\thttps://github.com/FelixKrueger/Bismark\n");
        return Ok(());
    }

    eprintln!("Methylation call infile:\t\t{}", cli.infile.display());
    eprintln!("Genome directory:\t\t\t>{}<", cli.genome_folder.display());
    eprintln!("Sample specified as NOMe-Seq\t\tyes (only reporting ACG and TCG context)");
    eprintln!("Optional GC context track:\t\tyes (NOMe-Seq; only reporting GCA, GCC and GCT context)");

    let genome = Genome::load(&cli.genome_folder)?;
    eprintln!("Stored sequence information of {} chromosomes/scaffolds in total\n", genome.sequences.len());

    process_file(&cli.infile, &genome, &cli.output_dir)?;
    Ok(())
}

fn process_file(infile: &Path, genome: &Genome, output_dir: &str) -> Result<()> {
    let reader: Box<dyn BufRead> = if infile.to_string_lossy().ends_with(".gz") {
        let f = std::fs::File::open(infile)?;
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::new(std::fs::File::open(infile)?))
    };

    // Derive output filename: replace existing extension with .manOwar.txt.gz
    let fname = infile.file_name().unwrap().to_string_lossy();
    let base = fname.trim_end_matches(".gz");
    let stem = if let Some(pos) = base.rfind('.') { &base[..pos] } else { base };
    let out_dir = if output_dir.is_empty() { String::new() }
                  else if output_dir.ends_with('/') { output_dir.to_string() }
                  else { format!("{output_dir}/") };
    let outpath = format!("{out_dir}{stem}.manOwar.txt.gz");

    let outfile = std::fs::File::create(&outpath)
        .with_context(|| format!("failed to create {outpath}"))?;
    let mut gz_out = GzEncoder::new(outfile, Compression::default());
    writeln!(gz_out, "ReadID\tChr\tStart\tEnd\tmeth_CG\tunmeth_CG\tmeth_GC\tunmeth_GC")?;
    eprintln!(">>> Writing genome-wide cytosine report to: {outpath} <<<\n");

    let mut current_id: Option<String> = None;
    let mut current_chr = String::new();
    let mut current_start: u64 = 0;
    let mut current_end: u64 = 0;
    // positions: BTreeMap<genomic_pos, (state b'+'/b'-', context_char)>
    let mut positions: BTreeMap<u64, (u8, u8)> = BTreeMap::new();

    let mut line_buf = String::new();
    let mut number_processed: u64 = 0;
    let mut line_reader = reader;

    loop {
        line_buf.clear();
        if line_reader.read_line(&mut line_buf)? == 0 { break; }
        let line = line_buf.trim_end_matches(['\n', '\r']);
        if line.starts_with("Bismark") || line.starts_with("ReadID") || line.is_empty() { continue; }

        let fields: Vec<&str> = line.splitn(9, '\t').collect();
        if fields.len() < 8 { continue; }

        let id      = fields[0];
        let state   = fields[1].as_bytes().first().copied().unwrap_or(b'-');
        let chr     = fields[2];
        let pos: u64 = fields[3].parse().unwrap_or(0);
        let ctx_char = fields[4].as_bytes().first().copied().unwrap_or(b'.');
        let start: u64 = fields[5].parse().unwrap_or(0);
        let end: u64   = fields[6].parse().unwrap_or(0);

        match &current_id {
            None => {
                current_id = Some(id.to_string());
                current_chr = chr.to_string();
                current_start = start;
                current_end = end;
                positions.insert(pos, (state, ctx_char));
            }
            Some(last) if last == id => {
                positions.insert(pos, (state, ctx_char));
            }
            Some(last) => {
                let last_id = last.clone();
                flush_read(&last_id, &current_chr, current_start, current_end,
                           &positions, genome, &mut gz_out)?;
                number_processed += 1;
                current_id = Some(id.to_string());
                current_chr = chr.to_string();
                current_start = start;
                current_end = end;
                positions.clear();
                positions.insert(pos, (state, ctx_char));
            }
        }
    }

    if let Some(last_id) = current_id {
        flush_read(&last_id, &current_chr, current_start, current_end,
                   &positions, genome, &mut gz_out)?;
        number_processed += 1;
    }

    gz_out.finish()?;

    eprintln!("Finished writing out NOMe-Seq specific filtering report (only reporting CGs in ACG and TCG context; reporting GCs only when not in CG context).");
    eprintln!("Processed {number_processed} reads in total.\n");
    Ok(())
}

/// Emit one line to the output for a completed read.
fn flush_read(
    id: &str,
    chr: &str,
    start: u64,
    end: u64,
    positions: &BTreeMap<u64, (u8, u8)>,
    genome: &Genome,
    out: &mut impl Write,
) -> Result<()> {
    let (len, offset): (u64, u64) = if end >= start {
        (end - start + 1, start)
    } else {
        (start - end + 1, end)
    };

    // Extended sequence: 2 extra bases on each side for trinucleotide context
    let ext_start = (offset as usize).saturating_sub(3); // 0-based
    let ext_len = len as usize + 6;
    let ext_seq = match genome.slice(chr, ext_start, ext_len) {
        Some(s) => s,
        None => return Ok(()),
    };

    let mut meth_cg: u32 = 0;
    let mut unmeth_cg: u32 = 0;
    let mut meth_gc: u32 = 0;
    let mut unmeth_gc: u32 = 0;

    for (i, &base) in ext_seq.iter().enumerate() {
        if base != b'C' && base != b'G' { continue; }
        if i + 3 > ext_seq.len() { continue; }

        // 1-based genomic position for this base
        let genomic_pos = ext_start as u64 + i as u64 + 1;

        let (tri_nt, upstream_ctx) = if base == b'C' {
            // Forward strand C: tri_nt[0..3] = this base + 2 downstream
            let t0 = ext_seq[i];
            let t1 = ext_seq.get(i+1).copied().unwrap_or(b'N');
            let t2 = ext_seq.get(i+2).copied().unwrap_or(b'N');
            // upstream_ctx: position i-1, i, i+1
            let u0 = if i > 0 { ext_seq[i-1] } else { b'N' };
            let u1 = t0;
            let u2 = t1;
            ([t0, t1, t2], [u0, u1, u2])
        } else {
            // Reverse strand G: reverse-complement of positions i+2, i+1, i
            let t0 = rc(ext_seq.get(i+2).copied().unwrap_or(b'N'));
            let t1 = rc(ext_seq.get(i+1).copied().unwrap_or(b'N'));
            let t2 = rc(ext_seq[i]);
            // upstream context = rc of positions i+3, i+2, i+1
            let u0 = rc(ext_seq.get(i+3).copied().unwrap_or(b'N'));
            let u1 = rc(ext_seq.get(i+2).copied().unwrap_or(b'N'));
            let u2 = rc(ext_seq.get(i+1).copied().unwrap_or(b'N'));
            ([t0, t1, t2], [u0, u1, u2])
        };

        // Determine cytosine context
        let is_cpg = tri_nt[1] == b'G';
        let is_chg = !is_cpg && tri_nt[2] == b'G';
        let is_chh = !is_cpg && !is_chg && !tri_nt.contains(&b'N');
        if !is_cpg && !is_chg && !is_chh { continue; }

        // Check if this genomic position has a methylation call in the read
        if let Some(&(state, _ctx_char)) = positions.get(&genomic_pos) {
            if is_cpg {
                // Only ACG or TCG upstream context
                if upstream_ctx[0] == b'A' || upstream_ctx[0] == b'T' {
                    if state == b'+' { meth_cg += 1; } else { unmeth_cg += 1; }
                }
            } else {
                // CHG or CHH: report GC context (upstream starts with GC)
                if upstream_ctx[0] == b'G' && upstream_ctx[1] == b'C' {
                    if state == b'+' { meth_gc += 1; } else { unmeth_gc += 1; }
                }
            }
        }
    }

    writeln!(out, "{id}\t{chr}\t{start}\t{end}\t{meth_cg}\t{unmeth_cg}\t{meth_gc}\t{unmeth_gc}")?;
    Ok(())
}

fn rc(b: u8) -> u8 {
    match b { b'A' => b'T', b'T' => b'A', b'G' => b'C', b'C' => b'G', _ => b'N' }
}
