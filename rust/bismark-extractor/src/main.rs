use std::collections::HashMap;
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use bismark_lib::bam_io::{bam_is_truncated, find_samtools, BamReader};
use bismark_lib::BISMARK_VERSION;
use clap::Parser;
use flate2::write::GzEncoder;
use flate2::Compression;
use rayon::prelude::*;

#[derive(Clone, Parser)]
#[command(
    name = "bismark_methylation_extractor",
    about = "Extract per-cytosine methylation information from Bismark SAM/BAM files",
    version = BISMARK_VERSION,
    disable_version_flag = true,
)]
struct Cli {
    /// Input Bismark SAM/BAM/CRAM file(s)
    files: Vec<PathBuf>,

    /// Single-end data (auto-detected if not specified)
    #[arg(short = 's', long = "single")]
    single: bool,

    /// Paired-end data (auto-detected if not specified)
    #[arg(short = 'p', long = "paired")]
    paired: bool,

    /// Comprehensive output: 3 files (CpG, CHG, CHH) instead of 12 strand-specific files
    #[arg(long = "comprehensive")]
    comprehensive: bool,

    /// Merge CHG and CHH into a single non-CpG context file
    #[arg(long = "merge_non_CpG")]
    merge_non_cpg: bool,

    /// Output directory
    #[arg(
        short = 'o',
        long = "output",
        alias = "output_dir",
        alias = "dir",
        default_value = ""
    )]
    output_dir: String,

    /// Omit header line from output files
    #[arg(long = "no_header")]
    no_header: bool,

    /// Ignore first N bp of each read (5' end)
    #[arg(long = "ignore", default_value_t = 0)]
    ignore: usize,

    /// Ignore last N bp of each read (3' end)
    #[arg(long = "ignore_3prime", default_value_t = 0)]
    ignore_3prime: usize,

    /// Ignore first N bp of Read 2 (5' end, PE only)
    #[arg(long = "ignore_r2", default_value_t = 0)]
    ignore_r2: usize,

    /// Ignore last N bp of Read 2 (3' end, PE only)
    #[arg(long = "ignore_3prime_r2", default_value_t = 0)]
    ignore_3prime_r2: usize,

    /// Do not extract overlapping PE reads twice (default: on)
    #[arg(long = "no_overlap")]
    no_overlap: bool,

    /// Extract overlapping PE reads from both mates
    #[arg(long = "include_overlap")]
    include_overlap: bool,

    /// Include all cytosine contexts in output (alias for comprehensive)
    #[arg(long = "CX", alias = "CX_context")]
    cx_context: bool,

    /// Compress output with gzip
    #[arg(long = "gzip")]
    gzip: bool,

    /// Yet Another Context Hunting Tool: output all Cs regardless of context
    #[arg(long = "yacht")]
    yacht: bool,

    /// Disable M-bias plot generation
    #[arg(long = "mbias_off")]
    mbias_off: bool,

    /// Only perform M-bias analysis, don't write methylation calls
    #[arg(long = "mbias_only")]
    mbias_only: bool,

    /// Path to samtools
    #[arg(long = "samtools_path")]
    samtools_path: Option<String>,

    /// Include read counts in bedGraph output (accepted for compatibility; counts are always included)
    #[arg(long = "counts")]
    counts: bool,

    /// Replace whitespace in read IDs with underscores (forwarded to bismark2bedGraph)
    #[arg(long = "remove_spaces")]
    remove_spaces: bool,

    /// Handle genomes with many scaffolds (forwarded to bismark2bedGraph)
    #[arg(long = "gazillion", alias = "scaffolds")]
    gazillion: bool,

    /// Sort in memory rather than using UNIX sort (forwarded to bismark2bedGraph)
    #[arg(long = "ample_memory")]
    ample_memory: bool,

    /// Print splitting report
    #[arg(long = "report")]
    report: bool,

    /// Number of parallel cores (accepted, currently processes single-threaded)
    #[arg(long = "multicore", alias = "parallel", default_value_t = 1)]
    multicore: u32,

    /// Also generate bedGraph and bismark.cov files via bismark2bedGraph
    #[arg(long = "bedGraph")]
    bedgraph: bool,

    /// Also generate genome-wide cytosine report via coverage2cytosine (requires --bedGraph and --genome_folder)
    #[arg(long = "cytosine_report")]
    cytosine_report: bool,

    /// Genome folder for --cytosine_report
    #[arg(short = 'g', long = "genome_folder")]
    genome_folder: Option<PathBuf>,

    /// Split cytosine report by chromosome (forwarded to coverage2cytosine)
    #[arg(long = "split_by_chromosome")]
    split_by_chromosome: bool,

    /// Use 0-based coordinates in cytosine report (forwarded to coverage2cytosine)
    #[arg(long = "zero_based")]
    zero_based: bool,

    /// Buffer size for sort (forwarded to bismark2bedGraph) [default: 2G]
    #[arg(long = "buffer_size", default_value = "2G")]
    buffer_size: String,

    /// Write UCSC-compatible bedGraph (forwarded to bismark2bedGraph)
    #[arg(long = "ucsc")]
    ucsc: bool,

    /// Minimum coverage to call methylation in bedGraph (forwarded to bismark2bedGraph) [default: 1]
    #[arg(long = "coverage_threshold", default_value_t = 1)]
    coverage_threshold: u32,

    /// Print version and exit
    #[arg(long = "version")]
    version: bool,
}

// ─── Cytosine context ────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum CytosineContext {
    CpG,
    CHG,
    CHH,
}

// ─── M-bias tracking ────────────────────────────────────────────────────────

#[derive(Default)]
struct MbiasPos {
    meth: u64,
    unmeth: u64,
}

type MbiasTable = HashMap<CytosineContext, Vec<MbiasPos>>; // indexed by read position (0-based)

fn mbias_add(table: &mut MbiasTable, ctx: CytosineContext, pos: usize, methylated: bool) {
    let vec = table.entry(ctx).or_default();
    if pos >= vec.len() {
        vec.resize_with(pos + 1, MbiasPos::default);
    }
    if methylated {
        vec[pos].meth += 1;
    } else {
        vec[pos].unmeth += 1;
    }
}

// ─── Output file handles ────────────────────────────────────────────────────

struct OutputFiles {
    // strand_specific[strand_idx][context_idx]: strand-specific mode
    strand_specific: [[Option<OutputTarget>; 3]; 4],
    // comprehensive[context_idx]: comprehensive/CX mode
    comprehensive: [Option<OutputTarget>; 3],
    // yacht/any_c
    any_c: Option<OutputTarget>,
    mode: OutputMode,
}

enum OutputTarget {
    File(Box<dyn Write + Send>),
    Buffer(Vec<u8>),
}

impl OutputTarget {
    fn buffer() -> Self {
        OutputTarget::Buffer(Vec::new())
    }

    fn buffer_bytes(&self) -> Option<&[u8]> {
        match self {
            OutputTarget::Buffer(v) => Some(v.as_slice()),
            OutputTarget::File(_) => None,
        }
    }
}

impl Write for OutputTarget {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            OutputTarget::File(w) => w.write(buf),
            OutputTarget::Buffer(v) => v.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            OutputTarget::File(w) => w.flush(),
            OutputTarget::Buffer(v) => v.flush(),
        }
    }
}

#[derive(Clone, Copy)]
enum OutputMode {
    StrandSpecific,     // 12 files (4 strands × 3 contexts)
    Comprehensive,      // 3 files (3 contexts across all strands)
    MergeNonCpG,        // 8 files (4 strands × [CpG, non-CpG])
    ComprehensiveMerge, // 2 files (CpG, non-CpG)
    Yacht,              // 1 file (all C)
}

impl OutputFiles {
    fn write_call(
        &mut self,
        id: &[u8],
        strand: &str,
        chr: &[u8],
        pos: u64,
        call: u8,
        read_start: u64,
        read_end: u64,
        sam_strand: u8,
    ) -> Result<()> {
        let ctx = context_of(call);
        let methylated = call.is_ascii_uppercase();
        let meth_strand = if methylated { b'+' } else { b'-' };
        let strand_idx = match strand {
            "OT" => 0,
            "CTOT" => 1,
            "CTOB" => 2,
            _ => 3,
        };
        let ctx_idx = match ctx {
            CytosineContext::CpG => 0,
            CytosineContext::CHG => 1,
            CytosineContext::CHH => 2,
        };

        match self.mode {
            OutputMode::StrandSpecific => {
                if let Some(ref mut w) = self.strand_specific[strand_idx][ctx_idx] {
                    write_line(w, id, meth_strand, chr, pos, call)?;
                }
            }
            OutputMode::Comprehensive => {
                if let Some(ref mut w) = self.comprehensive[ctx_idx] {
                    write_line(w, id, meth_strand, chr, pos, call)?;
                }
            }
            OutputMode::MergeNonCpG => {
                let ci = if ctx == CytosineContext::CpG { 0 } else { 1 };
                if let Some(ref mut w) = self.strand_specific[strand_idx][ci] {
                    write_line(w, id, meth_strand, chr, pos, call)?;
                }
            }
            OutputMode::ComprehensiveMerge => {
                let ci = if ctx == CytosineContext::CpG { 0 } else { 1 };
                if let Some(ref mut w) = self.comprehensive[ci] {
                    write_line(w, id, meth_strand, chr, pos, call)?;
                }
            }
            OutputMode::Yacht => {
                if let Some(ref mut w) = self.any_c {
                    w.write_all(id)?;
                    w.write_all(b"\t")?;
                    w.write_all(&[meth_strand, b'\t'])?;
                    w.write_all(chr)?;
                    w.write_all(b"\t")?;
                    write!(w, "{pos}\t")?;
                    w.write_all(&[call])?;
                    write!(w, "\t{read_start}\t{read_end}\t{}\n", sam_strand as char)?;
                }
            }
        }
        Ok(())
    }

    fn append_from(&mut self, other: &OutputFiles) -> Result<()> {
        for si in 0..4 {
            for ci in 0..3 {
                append_target(
                    &mut self.strand_specific[si][ci],
                    &other.strand_specific[si][ci],
                )?;
            }
        }
        for ci in 0..3 {
            append_target(&mut self.comprehensive[ci], &other.comprehensive[ci])?;
        }
        append_target(&mut self.any_c, &other.any_c)?;
        Ok(())
    }
}

fn append_target(dst: &mut Option<OutputTarget>, src: &Option<OutputTarget>) -> Result<()> {
    if let (Some(d), Some(s)) = (dst, src) {
        if let Some(bytes) = s.buffer_bytes() {
            d.write_all(bytes)?;
        }
    }
    Ok(())
}

fn write_line(
    w: &mut OutputTarget,
    id: &[u8],
    meth: u8,
    chr: &[u8],
    pos: u64,
    call: u8,
) -> Result<()> {
    w.write_all(id)?;
    w.write_all(b"\t")?;
    w.write_all(&[meth, b'\t'])?;
    w.write_all(chr)?;
    w.write_all(b"\t")?;
    write!(w, "{pos}\t")?;
    w.write_all(&[call, b'\n'])?;
    Ok(())
}

fn context_of(call: u8) -> CytosineContext {
    match call.to_ascii_uppercase() {
        b'Z' => CytosineContext::CpG,
        b'X' => CytosineContext::CHG,
        _ => CytosineContext::CHH,
    }
}

// ─── CIGAR helpers ──────────────────────────────────────────────────────────

fn expand_cigar(cigar: &[u8]) -> Vec<u8> {
    let mut ops = Vec::new();
    let mut n: usize = 0;
    for &b in cigar {
        if b.is_ascii_digit() {
            n = n * 10 + (b - b'0') as usize;
        } else {
            for _ in 0..n {
                ops.push(b);
            }
            n = 0;
        }
    }
    ops
}

fn apply_ignore_5prime(
    xm: &mut Vec<u8>,
    start: &mut u64,
    cigar_ops: &mut Vec<u8>,
    ignore: usize,
    forward: bool,
) {
    if ignore == 0 {
        return;
    }
    let trim = ignore.min(xm.len());
    if forward {
        let trimmed_ops: Vec<u8> = cigar_ops.drain(..trim.min(cigar_ops.len())).collect();
        let mut d = 0i64;
        for &op in &trimmed_ops {
            match op {
                b'D' | b'N' => d += 1,
                b'I' => d -= 1,
                _ => {}
            }
        }
        *start = (*start as i64 + ignore as i64 + d) as u64;
        xm.drain(..trim);
    } else {
        cigar_ops.truncate(cigar_ops.len().saturating_sub(trim));
        xm.drain(..trim);
    }
}

fn apply_ignore_3prime(xm: &mut Vec<u8>, cigar_ops: &mut Vec<u8>, ignore: usize) {
    if ignore == 0 {
        return;
    }
    let trim = ignore.min(xm.len());
    xm.truncate(xm.len().saturating_sub(trim));
    cigar_ops.truncate(cigar_ops.len().saturating_sub(trim));
}

fn mdn_count(cigar_ops: &[u8]) -> u64 {
    cigar_ops
        .iter()
        .filter(|&&b| matches!(b, b'M' | b'D' | b'N' | b'=' | b'X'))
        .count() as u64
}

fn apply_trimming(
    xm: &mut Vec<u8>,
    start: &mut u64,
    cigar_ops: &mut Vec<u8>,
    ignore_5prime: usize,
    ignore_3prime: usize,
    forward: bool,
) {
    let mut reverse_start_adjusted = false;

    if ignore_5prime > 0 {
        apply_ignore_5prime(xm, start, cigar_ops, ignore_5prime, forward);
        if !forward {
            *start += mdn_count(cigar_ops).saturating_sub(1);
            reverse_start_adjusted = true;
        }
    }

    if ignore_3prime > 0 {
        if !forward && !reverse_start_adjusted {
            *start += ignore_3prime.min(xm.len()) as u64;
        }
        apply_ignore_3prime(xm, cigar_ops, ignore_3prime);
    }

    if !forward && !reverse_start_adjusted {
        *start += mdn_count(cigar_ops).saturating_sub(1);
    }
}

// ─── Core extraction ────────────────────────────────────────────────────────

struct Counts {
    meth_cpg: u64,
    unmeth_cpg: u64,
    meth_chg: u64,
    unmeth_chg: u64,
    meth_chh: u64,
    unmeth_chh: u64,
    total: u64,
}

impl Counts {
    fn new() -> Self {
        Counts {
            meth_cpg: 0,
            unmeth_cpg: 0,
            meth_chg: 0,
            unmeth_chg: 0,
            meth_chh: 0,
            unmeth_chh: 0,
            total: 0,
        }
    }
    fn add_call(&mut self, call: u8) {
        match call {
            b'Z' => self.meth_cpg += 1,
            b'z' => self.unmeth_cpg += 1,
            b'X' => self.meth_chg += 1,
            b'x' => self.unmeth_chg += 1,
            b'H' => self.meth_chh += 1,
            b'h' => self.unmeth_chh += 1,
            _ => return,
        }
        self.total += 1;
    }

    fn merge(&mut self, other: &Counts) {
        self.meth_cpg += other.meth_cpg;
        self.unmeth_cpg += other.unmeth_cpg;
        self.meth_chg += other.meth_chg;
        self.unmeth_chg += other.unmeth_chg;
        self.meth_chh += other.meth_chh;
        self.unmeth_chh += other.unmeth_chh;
        self.total += other.total;
    }
}

fn merge_mbias(dst: &mut MbiasTable, src: &MbiasTable) {
    for (&ctx, src_vec) in src {
        let dst_vec = dst.entry(ctx).or_default();
        if dst_vec.len() < src_vec.len() {
            dst_vec.resize_with(src_vec.len(), MbiasPos::default);
        }
        for (i, src_pos) in src_vec.iter().enumerate() {
            dst_vec[i].meth += src_pos.meth;
            dst_vec[i].unmeth += src_pos.unmeth;
        }
    }
}

struct RecordGroup {
    first: Vec<u8>,
    second: Option<Vec<u8>>,
}

struct ChunkResult {
    out: OutputFiles,
    mbias1: MbiasTable,
    mbias2: MbiasTable,
    counts: Counts,
}

fn extract_calls(
    xm: &[u8],
    cigar: &[u8],
    start: u64,
    chr: &[u8],
    id: &[u8],
    bismark_strand: &str, // "OT" "CTOT" "CTOB" "OB"
    forward: bool,        // alignment is on + strand
    read_identity: u8,    // 1 or 2 for PE
    no_overlap: bool,
    overlap_limit: u64, // if no_overlap, stop when pos >= overlap_limit
    out: &mut OutputFiles,
    mbias: &mut MbiasTable,
    mbias_only: bool,
    counts: &mut Counts,
) -> Result<()> {
    if xm.is_empty() {
        return Ok(());
    }

    // For expanded CIGAR (all-M) the offset loop is a no-op; only needed for indels.
    // `cigar` here is the already-expanded array (e.g. [M, M, M, ...]).
    // We re-expand it via expand_cigar only if it actually contains ops that need tracking.
    let has_indels = cigar
        .iter()
        .any(|&b| matches!(b, b'I' | b'D' | b'N' | b'S'));
    let ops = if has_indels {
        cigar.to_vec()
    } else {
        Vec::new()
    };

    let mut pos_offset: i64 = 0;
    let mut ci: usize = 0;
    let read_span = mdn_count(cigar);
    let read_start = start;
    let read_end = if forward {
        start + read_span.saturating_sub(1)
    } else {
        start.saturating_sub(read_span.saturating_sub(1))
    };
    let sam_strand = if forward { b'+' } else { b'-' };

    for (i, &call) in xm.iter().enumerate() {
        // Advance cigar offset for indels.
        // For reverse reads (XM and CIGAR already reversed), I/D/N adjustments are mirrored.
        if !ops.is_empty() && ci < ops.len() {
            let op = ops[ci];
            ci += 1;
            if forward {
                match op {
                    b'I' => pos_offset -= 1,
                    b'D' | b'N' => pos_offset += 1,
                    _ => {}
                }
            } else {
                match op {
                    b'I' => pos_offset += 1,
                    b'D' | b'N' => pos_offset -= 1,
                    _ => {}
                }
            }
        }

        if !matches!(call, b'Z' | b'z' | b'X' | b'x' | b'H' | b'h') {
            continue;
        }

        // Forward: start = POS, pos goes up. Reverse: start = end_pos, pos goes down.
        let pos = if forward {
            (start as i64 + i as i64 + pos_offset) as u64
        } else {
            (start as i64 - i as i64 + pos_offset) as u64
        };

        // no_overlap check
        if no_overlap && read_identity == 2 {
            if forward {
                if pos >= overlap_limit {
                    break;
                }
            } else {
                if pos <= overlap_limit {
                    break;
                }
            }
        }

        counts.add_call(call);

        let ctx = context_of(call);
        let mbias_pos = i;
        if read_identity == 1 {
            mbias_add(mbias, ctx, mbias_pos, call.is_ascii_uppercase());
        } else {
            mbias_add(mbias, ctx, mbias_pos, call.is_ascii_uppercase());
        }

        if !mbias_only {
            out.write_call(
                id,
                bismark_strand,
                chr,
                pos,
                call,
                read_start,
                read_end,
                sam_strand,
            )?;
        }
    }
    Ok(())
}

// ─── SAM tag extraction ──────────────────────────────────────────────────────

fn find_tag<'a>(fields: &[&'a [u8]], tag: &[u8]) -> Option<&'a [u8]> {
    for &f in fields.iter().skip(11) {
        if f.starts_with(tag) {
            return Some(&f[tag.len()..]);
        }
    }
    None
}

fn determine_strand(xr: &[u8], xg: &[u8]) -> Option<(&'static str, bool)> {
    match (xr, xg) {
        (b"CT", b"CT") => Some(("OT", true)),
        (b"GA", b"CT") => Some(("CTOT", false)),
        (b"GA", b"GA") => Some(("CTOB", true)),
        (b"CT", b"GA") => Some(("OB", false)),
        _ => None,
    }
}

// ─── Output file opening ─────────────────────────────────────────────────────

fn new_writer(path: &str, gzip: bool, no_header: bool) -> Result<OutputTarget> {
    let f = std::fs::File::create(path).with_context(|| format!("creating {path}"))?;
    let mut w = if gzip {
        OutputTarget::File(Box::new(GzEncoder::new(f, Compression::default())))
    } else {
        OutputTarget::File(Box::new(BufWriter::new(f)))
    };
    if !no_header {
        writeln!(w, "Bismark methylation extractor version {BISMARK_VERSION}")?;
    }
    eprintln!("Writing result file: {path}");
    Ok(w)
}

fn new_buffer(no_header: bool) -> Result<OutputTarget> {
    let mut w = OutputTarget::buffer();
    if !no_header {
        writeln!(w, "Bismark methylation extractor version {BISMARK_VERSION}")?;
    }
    Ok(w)
}

/// Strip path components and the file extension, returning the bare base name.
/// E.g. "/data/sample.bam" → "sample", "run.sam.gz" → "run".
fn strip_ext(filename: &str) -> &str {
    let base = filename.split('/').last().unwrap_or(filename);
    let base = base.trim_end_matches(".gz");
    base.trim_end_matches(".bam")
        .trim_end_matches(".cram")
        .trim_end_matches(".sam")
        .trim_end_matches(".txt")
}

/// Prepend output_dir to the bare base name to form the output path prefix.
fn make_stem(filename: &str, output_dir: &str) -> String {
    format!("{}{}", output_dir, strip_ext(filename))
}

fn open_outputs(stem: &str, mode: OutputMode, gzip: bool, no_header: bool) -> Result<OutputFiles> {
    let ext = if gzip { ".txt.gz" } else { ".txt" };
    let contexts = ["CpG", "CHG", "CHH"];
    let strand_names_long = ["OT", "CTOT", "CTOB", "OB"];

    // stem = "<output_dir><basename>"; split so context prefix lands after the dir separator
    let (dir, base) = match stem.rfind('/') {
        Some(i) => (&stem[..=i], &stem[i + 1..]),
        None => ("", stem),
    };

    // Initialize with None
    const NONE_WRITER: Option<OutputTarget> = None;
    let mut ss: [[Option<OutputTarget>; 3]; 4] = [
        [NONE_WRITER; 3],
        [NONE_WRITER; 3],
        [NONE_WRITER; 3],
        [NONE_WRITER; 3],
    ];
    let mut comp: [Option<OutputTarget>; 3] = [NONE_WRITER; 3];
    let mut any_c: Option<OutputTarget> = None;

    match mode {
        OutputMode::StrandSpecific => {
            for si in 0..4 {
                for ci in 0..3 {
                    let path = format!(
                        "{dir}{}_{}_{base}{ext}",
                        contexts[ci], strand_names_long[si]
                    );
                    ss[si][ci] = Some(new_writer(&path, gzip, no_header)?);
                }
            }
        }
        OutputMode::Comprehensive => {
            for ci in 0..3 {
                let path = format!("{dir}{}_context_{base}{ext}", contexts[ci]);
                comp[ci] = Some(new_writer(&path, gzip, no_header)?);
            }
        }
        OutputMode::MergeNonCpG => {
            for si in 0..4 {
                let p0 = format!("{dir}CpG_{}_{base}{ext}", strand_names_long[si]);
                let p1 = format!("{dir}Non_CpG_{}_{base}{ext}", strand_names_long[si]);
                ss[si][0] = Some(new_writer(&p0, gzip, no_header)?);
                ss[si][1] = Some(new_writer(&p1, gzip, no_header)?);
            }
        }
        OutputMode::ComprehensiveMerge => {
            comp[0] = Some(new_writer(
                &format!("{dir}CpG_context_{base}{ext}"),
                gzip,
                no_header,
            )?);
            comp[1] = Some(new_writer(
                &format!("{dir}Non_CpG_context_{base}{ext}"),
                gzip,
                no_header,
            )?);
        }
        OutputMode::Yacht => {
            any_c = Some(new_writer(
                &format!("{dir}any_C_context_{base}{ext}"),
                gzip,
                no_header,
            )?);
        }
    }

    Ok(OutputFiles {
        strand_specific: ss,
        comprehensive: comp,
        any_c,
        mode,
    })
}

fn open_buffer_outputs(mode: OutputMode) -> Result<OutputFiles> {
    const NONE_WRITER: Option<OutputTarget> = None;
    let mut ss: [[Option<OutputTarget>; 3]; 4] = [
        [NONE_WRITER; 3],
        [NONE_WRITER; 3],
        [NONE_WRITER; 3],
        [NONE_WRITER; 3],
    ];
    let mut comp: [Option<OutputTarget>; 3] = [NONE_WRITER; 3];
    let mut any_c: Option<OutputTarget> = None;

    match mode {
        OutputMode::StrandSpecific => {
            for si in 0..4 {
                for ci in 0..3 {
                    ss[si][ci] = Some(new_buffer(true)?);
                }
            }
        }
        OutputMode::Comprehensive => {
            for ci in 0..3 {
                comp[ci] = Some(new_buffer(true)?);
            }
        }
        OutputMode::MergeNonCpG => {
            for si in 0..4 {
                ss[si][0] = Some(new_buffer(true)?);
                ss[si][1] = Some(new_buffer(true)?);
            }
        }
        OutputMode::ComprehensiveMerge => {
            comp[0] = Some(new_buffer(true)?);
            comp[1] = Some(new_buffer(true)?);
        }
        OutputMode::Yacht => {
            any_c = Some(new_buffer(true)?);
        }
    }

    Ok(OutputFiles {
        strand_specific: ss,
        comprehensive: comp,
        any_c,
        mode,
    })
}

// ─── bedGraph / cytosine-report dispatch ─────────────────────────────────────

/// Locate a Bismark tool binary: prefer one sitting next to the current exe,
/// fall back to searching PATH.
fn find_bismark_tool(name: &str) -> PathBuf {
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            let candidate = dir.join(name);
            if candidate.exists() {
                return candidate;
            }
        }
    }
    PathBuf::from(name)
}

/// Return the list of output files that `open_outputs` would have created for
/// a given stem/mode/gzip combination, so we can pass them to bismark2bedGraph.
fn collect_output_files(stem: &str, mode: OutputMode, gzip: bool) -> Vec<PathBuf> {
    let ext = if gzip { ".txt.gz" } else { ".txt" };
    let contexts = ["CpG", "CHG", "CHH"];
    let strands = ["OT", "CTOT", "CTOB", "OB"];

    let (dir, base) = match stem.rfind('/') {
        Some(i) => (&stem[..=i], &stem[i + 1..]),
        None => ("", stem),
    };

    let mut files = Vec::new();
    match mode {
        OutputMode::StrandSpecific => {
            for si in 0..4 {
                for ci in 0..3 {
                    files.push(PathBuf::from(format!("{dir}{}_{}_{base}{ext}", contexts[ci], strands[si])));
                }
            }
        }
        OutputMode::Comprehensive => {
            for ci in 0..3 {
                files.push(PathBuf::from(format!("{dir}{}_context_{base}{ext}", contexts[ci])));
            }
        }
        OutputMode::MergeNonCpG => {
            for si in 0..4 {
                files.push(PathBuf::from(format!("{dir}CpG_{}_{base}{ext}", strands[si])));
                files.push(PathBuf::from(format!("{dir}Non_CpG_{}_{base}{ext}", strands[si])));
            }
        }
        OutputMode::ComprehensiveMerge => {
            files.push(PathBuf::from(format!("{dir}CpG_context_{base}{ext}")));
            files.push(PathBuf::from(format!("{dir}Non_CpG_context_{base}{ext}")));
        }
        OutputMode::Yacht => {
            files.push(PathBuf::from(format!("{dir}any_C_context_{base}{ext}")));
        }
    }
    files
}

/// Call bismark2bedGraph on the collected extractor output files.
/// Returns the path to the resulting bismark.cov.gz file.
fn run_bedgraph(
    cli: &Cli,
    output_dir: &str,
    all_files: &[PathBuf],
    bare_stem: &str,
) -> Result<PathBuf> {
    let tool = find_bismark_tool("bismark2bedGraph");
    let bedgraph_name = format!("{bare_stem}.bedGraph");

    eprintln!("\n\nNow generating a bedGraph file from the methylation extractor output...\n");

    let mut cmd = std::process::Command::new(&tool);
    cmd.arg("--output").arg(&bedgraph_name);
    if !output_dir.is_empty() {
        cmd.arg("--dir").arg(output_dir.trim_end_matches('/'));
    }
    if cli.cx_context {
        cmd.arg("--CX_context");
    }
    cmd.arg("--cutoff").arg(cli.coverage_threshold.to_string());
    if cli.buffer_size != "2G" {
        cmd.arg("--buffer_size").arg(&cli.buffer_size);
    }
    if cli.ucsc {
        cmd.arg("--ucsc");
    }
    if cli.zero_based {
        cmd.arg("--zero_based");
    }
    if cli.remove_spaces {
        cmd.arg("--remove_spaces");
    }
    if cli.gazillion {
        cmd.arg("--gazillion");
    }
    if cli.ample_memory {
        cmd.arg("--ample_memory");
    }
    for f in all_files {
        cmd.arg(f);
    }

    let status = cmd
        .status()
        .with_context(|| format!("failed to run {}", tool.display()))?;
    if !status.success() {
        bail!("bismark2bedGraph failed with non-zero exit status");
    }

    // Derive the coverage file path that bismark2bedGraph will have written.
    let coverage_name = format!("{bare_stem}.bismark.cov.gz");
    let coverage_path = if output_dir.is_empty() {
        PathBuf::from(&coverage_name)
    } else {
        PathBuf::from(format!("{}/{coverage_name}", output_dir.trim_end_matches('/')))
    };
    Ok(coverage_path)
}

/// Call coverage2cytosine on the bismark.cov.gz produced by bismark2bedGraph.
fn run_cytosine_report(
    cli: &Cli,
    output_dir: &str,
    coverage_file: &Path,
    bare_stem: &str,
) -> Result<()> {
    let tool = find_bismark_tool("coverage2cytosine");

    let cytosine_out = if cli.cx_context {
        format!("{bare_stem}.CX_report.txt")
    } else {
        format!("{bare_stem}.CpG_report.txt")
    };

    let genome_folder = cli
        .genome_folder
        .as_ref()
        .context("--genome_folder is required when --cytosine_report is set")?;

    eprintln!("\n\nNow generating a genome-wide cytosine methylation report...\n");

    let mut cmd = std::process::Command::new(&tool);
    cmd.arg("--output").arg(&cytosine_out);
    if !output_dir.is_empty() {
        cmd.arg("--dir").arg(output_dir.trim_end_matches('/'));
    }
    cmd.arg("--genome_folder").arg(genome_folder);
    if cli.zero_based {
        cmd.arg("--zero_based");
    }
    if cli.cx_context {
        cmd.arg("--CX_context");
    }
    if cli.split_by_chromosome {
        cmd.arg("--split_by_chromosome");
    }
    if cli.gzip {
        cmd.arg("--gzip");
    }
    cmd.arg(coverage_file);

    let status = cmd
        .status()
        .with_context(|| format!("failed to run {}", tool.display()))?;
    if !status.success() {
        bail!("coverage2cytosine failed with non-zero exit status");
    }
    Ok(())
}

// ─── Main processing ─────────────────────────────────────────────────────────

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.version {
        println!(
            "\n\n                    Bismark - Bisulfite Mapper and Methylation Caller.\n\n                         Bismark Methylation Extractor Version: {}\n                 Copyright 2010-22 Felix Krueger, Altos Bioinformatics\n                          https://github.com/FelixKrueger/Bismark\n\n",
            BISMARK_VERSION
        );
        return Ok(());
    }

    if cli.files.is_empty() {
        bail!("No input files specified. Use --help for usage.");
    }

    if cli.single && cli.paired {
        bail!("Cannot specify both --single and --paired");
    }

    if cli.mbias_only && cli.mbias_off {
        bail!("Options '--mbias_only' and '--mbias_off' are not compatible");
    }

    if cli.mbias_only && cli.yacht {
        bail!("The option '--yacht' does not work together with '--mbias_only'");
    }

    if cli.multicore == 0 {
        bail!("Core usage needs to be set to 1 or more");
    }

    if cli.cytosine_report && !cli.bedgraph {
        bail!("--cytosine_report requires --bedGraph");
    }
    if cli.cytosine_report && cli.genome_folder.is_none() {
        bail!("--cytosine_report requires --genome_folder");
    }

    let output_dir = normalise_dir(&cli.output_dir);

    let samtools = find_samtools(cli.samtools_path.as_deref())?;

    // Determine output mode
    let mode = if cli.yacht {
        OutputMode::Yacht
    } else if cli.comprehensive && cli.merge_non_cpg {
        OutputMode::ComprehensiveMerge
    } else if cli.comprehensive || cli.cx_context {
        OutputMode::Comprehensive
    } else if cli.merge_non_cpg {
        OutputMode::MergeNonCpG
    } else {
        OutputMode::StrandSpecific
    };

    let mut all_output_files: Vec<PathBuf> = Vec::new();
    let mut first_bare_stem: Option<String> = None;

    for file in &cli.files.clone() {
        let bare = strip_ext(file.to_str().unwrap_or("")).to_string();
        if first_bare_stem.is_none() {
            first_bare_stem = Some(bare.clone());
        }
        process_file(file, &cli, &samtools, &output_dir, mode)?;
        if cli.bedgraph {
            let stem = format!("{output_dir}{bare}");
            all_output_files.extend(collect_output_files(&stem, mode, cli.gzip));
        }
    }

    if cli.bedgraph {
        if let Some(ref bare_stem) = first_bare_stem {
            let coverage_path =
                run_bedgraph(&cli, output_dir.trim_end_matches('/'), &all_output_files, bare_stem)?;
            if cli.cytosine_report {
                run_cytosine_report(
                    &cli,
                    output_dir.trim_end_matches('/'),
                    &coverage_path,
                    bare_stem,
                )?;
            }
        }
    }

    Ok(())
}

fn process_file(
    path: &Path,
    cli: &Cli,
    samtools: &str,
    output_dir: &str,
    mode: OutputMode,
) -> Result<()> {
    let filename = path.to_str().unwrap_or("");

    if bam_is_truncated(samtools, path) {
        bail!("File {} appears truncated", path.display());
    }

    // Auto-detect SE/PE
    let is_paired = if cli.single {
        false
    } else if cli.paired {
        true
    } else {
        detect_is_paired(samtools, path)?
    };
    if cli.include_overlap && !is_paired {
        bail!("The option '--include_overlap' can only be specified for paired-end input");
    }
    let no_overlap = is_paired && !cli.include_overlap;

    let bare = strip_ext(filename);               // base name without path or extension
    let stem = format!("{output_dir}{bare}");     // full output path prefix
    let mut out = open_outputs(&stem, mode, cli.gzip, cli.no_header)?;
    let mut mbias1: MbiasTable = HashMap::new();
    let mut mbias2: MbiasTable = HashMap::new();
    let mut counts = Counts::new();

    eprintln!("Now reading in Bismark result file {filename}");

    let groups = read_record_groups(samtools, path, is_paired)?;
    let line_count = groups.len() as u64;

    if cli.multicore > 1 && groups.len() > 1 {
        let threads = cli.multicore as usize;
        eprintln!("Processing {line_count} record groups with {threads} worker threads");
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .context("building Rayon thread pool")?;
        let chunk_results: Vec<Result<ChunkResult>> = pool.install(|| {
            (0..threads)
                .into_par_iter()
                .map(|worker_idx| {
                    let mut chunk_out = open_buffer_outputs(mode)?;
                    let mut chunk_mbias1: MbiasTable = HashMap::new();
                    let mut chunk_mbias2: MbiasTable = HashMap::new();
                    let mut chunk_counts = Counts::new();
                    for group in groups.iter().skip(worker_idx).step_by(threads) {
                        process_record_group(
                            group,
                            is_paired,
                            cli,
                            no_overlap,
                            &mut chunk_out,
                            &mut chunk_mbias1,
                            &mut chunk_mbias2,
                            &mut chunk_counts,
                        )?;
                    }
                    Ok(ChunkResult {
                        out: chunk_out,
                        mbias1: chunk_mbias1,
                        mbias2: chunk_mbias2,
                        counts: chunk_counts,
                    })
                })
                .collect()
        });

        for chunk_result in chunk_results {
            let chunk_result = chunk_result?;
            out.append_from(&chunk_result.out)?;
            merge_mbias(&mut mbias1, &chunk_result.mbias1);
            merge_mbias(&mut mbias2, &chunk_result.mbias2);
            counts.merge(&chunk_result.counts);
        }
    } else {
        for (idx, group) in groups.iter().enumerate() {
            let processed = idx as u64 + 1;
            if processed % 500_000 == 0 {
                eprintln!("Processed {processed} lines");
            }
            process_record_group(
                group,
                is_paired,
                cli,
                no_overlap,
                &mut out,
                &mut mbias1,
                &mut mbias2,
                &mut counts,
            )?;
        }
    }

    // Write M-bias report
    if !cli.mbias_off {
        write_mbias_report(&stem, &mbias1, &mbias2, is_paired)?;
    }

    // Write splitting report
    if cli.report {
        write_splitting_report(&stem, &counts, is_paired)?;
    }

    eprintln!("\nTotal {} sequences analysed.", line_count);
    eprintln!("Methylated CpGs: {}", counts.meth_cpg);
    eprintln!("Unmethylated CpGs: {}", counts.unmeth_cpg);

    Ok(())
}

fn read_record_groups(samtools: &str, path: &Path, is_paired: bool) -> Result<Vec<RecordGroup>> {
    let mut reader = BamReader::open(samtools, path, &[])?;
    let mut buf = Vec::with_capacity(8192);
    let mut groups = Vec::new();

    loop {
        buf.clear();
        let n = reader.lines().read_until(b'\n', &mut buf)?;
        if n == 0 {
            break;
        }
        while buf.last() == Some(&b'\n') || buf.last() == Some(&b'\r') {
            buf.pop();
        }
        if buf.starts_with(b"@") || buf.is_empty() {
            continue;
        }

        if is_paired {
            let mut r2_buf = Vec::with_capacity(8192);
            loop {
                r2_buf.clear();
                let n = reader.lines().read_until(b'\n', &mut r2_buf)?;
                if n == 0 {
                    groups.push(RecordGroup {
                        first: buf.clone(),
                        second: None,
                    });
                    reader.finish()?;
                    return Ok(groups);
                }
                while r2_buf.last() == Some(&b'\n') || r2_buf.last() == Some(&b'\r') {
                    r2_buf.pop();
                }
                if r2_buf.starts_with(b"@") || r2_buf.is_empty() {
                    continue;
                }
                break;
            }
            groups.push(RecordGroup {
                first: buf.clone(),
                second: Some(r2_buf),
            });
        } else {
            groups.push(RecordGroup {
                first: buf.clone(),
                second: None,
            });
        }
    }
    reader.finish()?;
    Ok(groups)
}

fn process_record_group(
    group: &RecordGroup,
    is_paired: bool,
    cli: &Cli,
    no_overlap: bool,
    out: &mut OutputFiles,
    mbias1: &mut MbiasTable,
    mbias2: &mut MbiasTable,
    counts: &mut Counts,
) -> Result<()> {
    let buf = &group.first;

    let fields: Vec<&[u8]> = buf.split(|&b| b == b'\t').collect();
    if fields.len() < 11 {
        return Ok(());
    }

    let id = fields[0];
    let chr = fields[2];
    let pos: u64 = std::str::from_utf8(fields[3])
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    let cigar = fields[5];

    let xm = match find_tag(&fields, b"XM:Z:") {
        Some(x) => x,
        None => return Ok(()),
    };
    let xr = match find_tag(&fields, b"XR:Z:") {
        Some(x) => x,
        None => return Ok(()),
    };
    let xg = match find_tag(&fields, b"XG:Z:") {
        Some(x) => x,
        None => return Ok(()),
    };

    let (bismark_strand, forward) = match determine_strand(xr, xg) {
        Some(s) => s,
        None => return Ok(()),
    };

    let mut xm_vec: Vec<u8> = xm.to_vec();
    if !forward {
        xm_vec.reverse();
    }

    let mut cigar_ops = if cigar != b"*" {
        expand_cigar(cigar)
    } else {
        Vec::new()
    };
    if !forward {
        cigar_ops.reverse();
    }

    let mut start = pos;

    if is_paired {
        process_pair(
            group.second.as_deref(),
            id,
            chr,
            start,
            cigar_ops,
            xm_vec,
            bismark_strand,
            forward,
            cli,
            no_overlap,
            out,
            mbias1,
            mbias2,
            counts,
        )?;
    } else {
        apply_trimming(
            &mut xm_vec,
            &mut start,
            &mut cigar_ops,
            cli.ignore,
            cli.ignore_3prime,
            forward,
        );
        extract_calls(
            &xm_vec,
            &cigar_ops,
            start,
            chr,
            id,
            bismark_strand,
            forward,
            1,
            false,
            0,
            out,
            mbias1,
            cli.mbias_only,
            counts,
        )?;
    }
    Ok(())
}

fn process_pair(
    r2_buf: Option<&[u8]>,
    id1: &[u8],
    chr: &[u8],
    start_r1: u64,
    mut cigar_ops_r1: Vec<u8>,
    mut xm_r1: Vec<u8>,
    bismark_strand: &str,
    forward: bool, // R1 forward
    cli: &Cli,
    no_overlap: bool,
    out: &mut OutputFiles,
    mbias1: &mut MbiasTable,
    mbias2: &mut MbiasTable,
    counts: &mut Counts,
) -> Result<()> {
    let Some(r2_buf) = r2_buf else {
        return Ok(());
    };

    let fields2: Vec<&[u8]> = r2_buf.split(|&b| b == b'\t').collect();
    if fields2.len() < 11 {
        return Ok(());
    }

    let id2 = fields2[0];
    let chr2 = fields2[2];
    let pos2: u64 = std::str::from_utf8(fields2[3])
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    let cigar2 = fields2[5];

    let xm2 = match find_tag(&fields2, b"XM:Z:") {
        Some(x) => x,
        None => return Ok(()),
    };
    let xr2 = match find_tag(&fields2, b"XR:Z:") {
        Some(x) => x,
        None => return Ok(()),
    };
    let xg2 = match find_tag(&fields2, b"XG:Z:") {
        Some(x) => x,
        None => return Ok(()),
    };

    let (strand2, forward2) = match determine_strand(xr2, xg2) {
        Some(s) => s,
        None => return Ok(()),
    };

    let mut xm2_vec: Vec<u8> = xm2.to_vec();
    if !forward2 {
        xm2_vec.reverse();
    }
    let mut cigar2_ops = if cigar2 != b"*" {
        expand_cigar(cigar2)
    } else {
        Vec::new()
    };
    if !forward2 {
        cigar2_ops.reverse();
    }
    let mut start_r2 = pos2;

    let mut start_r1 = start_r1;
    apply_trimming(
        &mut xm_r1,
        &mut start_r1,
        &mut cigar_ops_r1,
        cli.ignore,
        cli.ignore_3prime,
        forward,
    );
    apply_trimming(
        &mut xm2_vec,
        &mut start_r2,
        &mut cigar2_ops,
        cli.ignore_r2,
        cli.ignore_3prime_r2,
        forward2,
    );

    let mdn1 = mdn_count(&cigar_ops_r1);
    let (r1_start_eff, r2_start_eff, end_r1) = if forward {
        (start_r1, start_r2, start_r1 + mdn1.saturating_sub(1))
    } else {
        (
            start_r1,
            start_r2,
            start_r1.saturating_sub(mdn1.saturating_sub(1)),
        )
    };

    // Extract R1
    extract_calls(
        &xm_r1,
        &cigar_ops_r1,
        r1_start_eff,
        chr,
        id1,
        bismark_strand,
        forward,
        1,
        false,
        0,
        out,
        mbias1,
        cli.mbias_only,
        counts,
    )?;

    // Extract R2 (with optional no_overlap)
    extract_calls(
        &xm2_vec,
        &cigar2_ops,
        r2_start_eff,
        chr2,
        id2,
        strand2,
        forward2,
        2,
        no_overlap,
        end_r1,
        out,
        mbias2,
        cli.mbias_only,
        counts,
    )?;

    Ok(())
}

// ─── M-bias report ────────────────────────────────────────────────────────────

fn write_mbias_report(
    stem: &str,
    mbias1: &MbiasTable,
    mbias2: &MbiasTable,
    is_paired: bool,
) -> Result<()> {
    let path = format!("{stem}.M-bias.txt");
    let mut f = std::fs::File::create(&path).with_context(|| format!("creating {path}"))?;

    writeln!(f, "CpG context\tRead 1")?;
    writeln!(
        f,
        "position\tcount methylated\tcount unmethylated\t% methylation\tcoverage"
    )?;
    write_mbias_context(&mut f, mbias1, CytosineContext::CpG)?;

    if is_paired {
        writeln!(f, "\nCpG context\tRead 2")?;
        writeln!(
            f,
            "position\tcount methylated\tcount unmethylated\t% methylation\tcoverage"
        )?;
        write_mbias_context(&mut f, mbias2, CytosineContext::CpG)?;
    }

    Ok(())
}

fn write_mbias_context(
    f: &mut std::fs::File,
    table: &MbiasTable,
    ctx: CytosineContext,
) -> Result<()> {
    if let Some(vec) = table.get(&ctx) {
        for (i, p) in vec.iter().enumerate() {
            let cov = p.meth + p.unmeth;
            let pct = if cov > 0 {
                format!("{:.2}", p.meth as f64 / cov as f64 * 100.0)
            } else {
                "0.00".into()
            };
            writeln!(f, "{}\t{}\t{}\t{pct}\t{cov}", i + 1, p.meth, p.unmeth)?;
        }
    }
    Ok(())
}

// ─── Splitting report ─────────────────────────────────────────────────────────

fn write_splitting_report(stem: &str, counts: &Counts, is_paired: bool) -> Result<()> {
    let path = format!("{stem}_splitting_report.txt");
    let mut f = std::fs::File::create(&path).with_context(|| format!("creating {path}"))?;

    let mode = if is_paired {
        "paired-end"
    } else {
        "single-end"
    };
    writeln!(f, "Bismark Extractor Version: {BISMARK_VERSION}")?;
    writeln!(f, "Bismark result file: {mode} (SAM format)")?;
    writeln!(f)?;
    writeln!(f, "Final Cytosine Methylation Report")?;
    writeln!(f, "=================================")?;
    let total = counts.meth_cpg
        + counts.unmeth_cpg
        + counts.meth_chg
        + counts.unmeth_chg
        + counts.meth_chh
        + counts.unmeth_chh;
    writeln!(f, "Total number of C's analysed:\t{total}")?;
    writeln!(f)?;
    writeln!(
        f,
        "Total methylated C's in CpG context:\t{}",
        counts.meth_cpg
    )?;
    writeln!(
        f,
        "Total methylated C's in CHG context:\t{}",
        counts.meth_chg
    )?;
    writeln!(
        f,
        "Total methylated C's in CHH context:\t{}",
        counts.meth_chh
    )?;
    writeln!(f)?;
    writeln!(
        f,
        "Total unmethylated C's in CpG context:\t{}",
        counts.unmeth_cpg
    )?;
    writeln!(
        f,
        "Total unmethylated C's in CHG context:\t{}",
        counts.unmeth_chg
    )?;
    writeln!(
        f,
        "Total unmethylated C's in CHH context:\t{}",
        counts.unmeth_chh
    )?;
    writeln!(f)?;
    let cpg_total = counts.meth_cpg + counts.unmeth_cpg;
    let chg_total = counts.meth_chg + counts.unmeth_chg;
    let chh_total = counts.meth_chh + counts.unmeth_chh;
    let pct_cpg = if cpg_total > 0 {
        format!("{:.1}", counts.meth_cpg as f64 / cpg_total as f64 * 100.0)
    } else {
        "N/A".into()
    };
    let pct_chg = if chg_total > 0 {
        format!("{:.1}", counts.meth_chg as f64 / chg_total as f64 * 100.0)
    } else {
        "N/A".into()
    };
    let pct_chh = if chh_total > 0 {
        format!("{:.1}", counts.meth_chh as f64 / chh_total as f64 * 100.0)
    } else {
        "N/A".into()
    };
    writeln!(f, "C methylated in CpG context:\t{pct_cpg}%")?;
    writeln!(f, "C methylated in CHG context:\t{pct_chg}%")?;
    writeln!(f, "C methylated in CHH context:\t{pct_chh}%")?;
    Ok(())
}

// ─── Library-type detection ───────────────────────────────────────────────────

fn detect_is_paired(samtools: &str, path: &Path) -> Result<bool> {
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
        if s.contains(" -1 ") && s.contains(" -2 ") {
            return Ok(true);
        }
        return Ok(false);
    }
    Ok(false)
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

    // ─── expand_cigar ────────────────────────────────────────────────────────

    #[test]
    fn test_expand_cigar_simple_match() {
        assert_eq!(expand_cigar(b"10M"), vec![b'M'; 10]);
    }

    #[test]
    fn test_expand_cigar_mixed() {
        // 3M1I2D = MMM I DD
        let ops = expand_cigar(b"3M1I2D");
        assert_eq!(ops, b"MMMIDD");
    }

    #[test]
    fn test_expand_cigar_empty() {
        assert!(expand_cigar(b"").is_empty());
    }

    #[test]
    fn test_expand_cigar_two_digit() {
        assert_eq!(expand_cigar(b"12M"), vec![b'M'; 12]);
    }

    // ─── mdn_count ───────────────────────────────────────────────────────────

    #[test]
    fn test_mdn_count_all_match() {
        let ops = expand_cigar(b"10M");
        assert_eq!(mdn_count(&ops), 10);
    }

    #[test]
    fn test_mdn_count_with_insertion() {
        // 8M2I: insertions don't consume reference
        let ops = expand_cigar(b"8M2I");
        assert_eq!(mdn_count(&ops), 8);
    }

    #[test]
    fn test_mdn_count_deletion() {
        // 5M2D3M: M and D both consume reference
        let ops = expand_cigar(b"5M2D3M");
        assert_eq!(mdn_count(&ops), 10);
    }

    #[test]
    fn test_mdn_count_with_softclip() {
        // 2S8M: S does not consume reference
        let ops = expand_cigar(b"2S8M");
        assert_eq!(mdn_count(&ops), 8);
    }

    // ─── context_of ──────────────────────────────────────────────────────────

    #[test]
    fn test_context_of_meth_cpg() {
        assert_eq!(context_of(b'Z'), CytosineContext::CpG);
    }

    #[test]
    fn test_context_of_unmeth_cpg() {
        assert_eq!(context_of(b'z'), CytosineContext::CpG);
    }

    #[test]
    fn test_context_of_meth_chg() {
        assert_eq!(context_of(b'X'), CytosineContext::CHG);
    }

    #[test]
    fn test_context_of_unmeth_chg() {
        assert_eq!(context_of(b'x'), CytosineContext::CHG);
    }

    #[test]
    fn test_context_of_chh() {
        assert_eq!(context_of(b'H'), CytosineContext::CHH);
        assert_eq!(context_of(b'h'), CytosineContext::CHH);
        assert_eq!(context_of(b'U'), CytosineContext::CHH);
        assert_eq!(context_of(b'.'), CytosineContext::CHH);
    }

    // ─── determine_strand ────────────────────────────────────────────────────

    #[test]
    fn test_determine_strand_ot() {
        assert_eq!(determine_strand(b"CT", b"CT"), Some(("OT", true)));
    }

    #[test]
    fn test_determine_strand_ctot() {
        assert_eq!(determine_strand(b"GA", b"CT"), Some(("CTOT", false)));
    }

    #[test]
    fn test_determine_strand_ctob() {
        assert_eq!(determine_strand(b"GA", b"GA"), Some(("CTOB", true)));
    }

    #[test]
    fn test_determine_strand_ob() {
        assert_eq!(determine_strand(b"CT", b"GA"), Some(("OB", false)));
    }

    #[test]
    fn test_determine_strand_invalid() {
        assert_eq!(determine_strand(b"XX", b"YY"), None);
    }

    // ─── strip_ext ───────────────────────────────────────────────────────────

    #[test]
    fn test_strip_ext_plain_bam() {
        assert_eq!(strip_ext("sample.bam"), "sample");
    }

    #[test]
    fn test_strip_ext_gzipped_sam() {
        // .gz stripped first, then .sam — result is bare name
        assert_eq!(strip_ext("sample.sam.gz"), "sample");
    }

    #[test]
    fn test_strip_ext_path_component() {
        assert_eq!(strip_ext("/data/run/sample.bam"), "sample");
    }

    // ─── make_stem ───────────────────────────────────────────────────────────

    #[test]
    fn test_make_stem_prepends_output_dir() {
        assert_eq!(make_stem("sample.bam", "/out/"), "/out/sample");
    }

    #[test]
    fn test_make_stem_empty_dir() {
        assert_eq!(make_stem("sample.bam", ""), "sample");
    }

    // ─── normalise_dir ───────────────────────────────────────────────────────

    #[test]
    fn test_normalise_dir_empty() {
        assert_eq!(normalise_dir(""), "");
    }

    #[test]
    fn test_normalise_dir_no_slash() {
        assert_eq!(normalise_dir("/out"), "/out/");
    }

    #[test]
    fn test_normalise_dir_already_slash() {
        assert_eq!(normalise_dir("/out/"), "/out/");
    }

    // ─── apply_ignore_5prime / apply_ignore_3prime ───────────────────────────

    #[test]
    fn test_apply_ignore_5prime_forward() {
        let mut xm = b"ZzXxHh".to_vec();
        let mut start = 100u64;
        let mut ops = vec![b'M'; 6];
        apply_ignore_5prime(&mut xm, &mut start, &mut ops, 2, true);
        assert_eq!(xm, b"XxHh");
        assert_eq!(start, 102);
        assert_eq!(ops.len(), 4);
    }

    #[test]
    fn test_apply_ignore_5prime_reverse() {
        // For reverse reads, trimming 5' means truncating the end of the (already-reversed) arrays
        let mut xm = b"hHxXzZ".to_vec();
        let mut start = 200u64;
        let mut ops = vec![b'M'; 6];
        apply_ignore_5prime(&mut xm, &mut start, &mut ops, 2, false);
        // drain(..2) from front of reversed arrays; start unchanged for reverse
        assert_eq!(xm, b"xXzZ");
        assert_eq!(ops.len(), 4);
    }

    #[test]
    fn test_apply_ignore_3prime() {
        let mut xm = b"ZzXxHh".to_vec();
        let mut ops = vec![b'M'; 6];
        apply_ignore_3prime(&mut xm, &mut ops, 2);
        assert_eq!(xm, b"ZzXx");
        assert_eq!(ops.len(), 4);
    }
}
