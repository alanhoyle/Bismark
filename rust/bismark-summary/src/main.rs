use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use bismark_lib::BISMARK_VERSION;
use clap::Parser;

const SUMMARY_TEMPLATE: &str = include_str!("../../../plotly/bismark_summary_template.tpl");
const PLOT_LY:          &str = include_str!("../../../plotly/plot.ly");
const BISMARK_LOGO:     &str = include_str!("../../../plotly/bismark.logo");

#[derive(Parser)]
#[command(
    name = "bismark2summary",
    about = "Generate HTML summary report across multiple Bismark BAM files",
    version = BISMARK_VERSION,
    disable_version_flag = true,
)]
struct Cli {
    /// BAM files to process (auto-detected if not specified)
    #[arg()]
    bam_files: Vec<PathBuf>,

    #[arg(short = 'o', long = "basename", default_value = "bismark_summary_report")]
    basename: String,

    #[arg(long = "title", default_value = "Bismark Summary Report")]
    title: String,

    #[arg(long = "verbose")]
    verbose: bool,

    #[arg(long = "version")]
    version: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.version {
        println!(
            "\n\n                         Bismark HTML Summary Report\n\n                         bismark2summary version: {}\n                   Copyright 2010-25 Felix Krueger, Altos Bioinformatics\n                          https://github.com/FelixKrueger/Bismark\n\n",
            BISMARK_VERSION
        );
        return Ok(());
    }

    let bam_files: Vec<PathBuf> = if !cli.bam_files.is_empty() {
        cli.bam_files.clone()
    } else {
        auto_detect_bam_files()?
    };

    if bam_files.is_empty() {
        bail!("No Bismark BAM files found. Please specify BAM files or run in a Bismark output directory.");
    }

    eprintln!("Generating Bismark summary report from {} BAM file(s)...", bam_files.len());

    let mut samples: Vec<Sample> = Vec::new();
    let mut csv_rows: Vec<String> = Vec::new();

    csv_rows.push("File\tTotal Reads\tAligned Reads\tUnaligned Reads\tAmbiguously Aligned Reads\tNo Genomic Sequence\tDuplicate Reads (removed)\tUnique Reads (remaining)\tTotal Cs\tMethylated CpGs\tUnmethylated CpGs\tMethylated chgs\tUnmethylated chgs\tMethylated CHHs\tUnmethylated CHHs".to_string());

    for bam in &bam_files {
        let sample = read_sample(bam, cli.verbose)?;
        let row = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            bam.display(),
            sample.total_reads, sample.aligned_reads, sample.unaligned,
            sample.ambig_reads, sample.no_seq_reads, sample.dup_reads, sample.unique_reads,
            sample.total_c,
            sample.meth_cpg, sample.unmeth_cpg,
            sample.meth_chg, sample.unmeth_chg,
            sample.meth_chh, sample.unmeth_chh,
        );
        csv_rows.push(row);

        let skip = (sample.meth_cpg == 0 && sample.unmeth_cpg == 0)
                || (sample.meth_chg == 0 && sample.unmeth_chg == 0)
                || (sample.meth_chh == 0 && sample.unmeth_chh == 0);
        if skip {
            eprintln!("Excluding sample {} from plotting (no calls in some context)", sample.name);
        } else {
            samples.push(sample);
        }
    }

    let txt_fn = format!("{}.txt", cli.basename);
    let mut f = fs::File::create(&txt_fn).with_context(|| format!("creating {txt_fn}"))?;
    for row in &csv_rows { writeln!(f, "{row}")?; }
    drop(f);

    let html = build_html(&samples, &cli.title, cli.verbose)?;
    let html_fn = format!("{}.html", cli.basename);
    let mut hf = fs::File::create(&html_fn).with_context(|| format!("creating {html_fn}"))?;
    hf.write_all(html.as_bytes())?;

    eprintln!("\nWrote Bismark project summary to >> {html_fn} <<\n");
    Ok(())
}

// ─── Per-sample data ─────────────────────────────────────────────────────────

struct Sample {
    name:          String,
    total_reads:   u64,
    aligned_reads: u64,
    unaligned:     u64,
    ambig_reads:   u64,
    no_seq_reads:  u64,
    dup_reads:     u64,
    unique_reads:  u64,
    total_c:       u64,
    meth_cpg:      u64,
    unmeth_cpg:    u64,
    meth_chg:      u64,
    unmeth_chg:    u64,
    meth_chh:      u64,
    unmeth_chh:    u64,
}

fn read_sample(bam: &Path, verbose: bool) -> Result<Sample> {
    let bam_str = bam.to_str().unwrap_or("");
    let base = bam_str.trim_end_matches(".bam");
    let paired_end = base.ends_with("_pe");
    let base_stripped = if paired_end { &base[..base.len() - 3] } else { base };

    let report_path = if paired_end {
        PathBuf::from(format!("{base_stripped}_PE_report.txt"))
    } else {
        PathBuf::from(format!("{base_stripped}_SE_report.txt"))
    };

    if !report_path.exists() {
        bail!("Could not find Bismark report: {}", report_path.display());
    }
    eprintln!(">> Reading from Bismark report: {}", report_path.display());

    let mut total_reads = 0u64;
    let mut aligned_reads = 0u64;
    let mut unaligned = 0u64;
    let mut ambig_reads = 0u64;
    let mut no_seq_reads = 0u64;
    let mut total_c = 0u64;
    let mut meth_cpg = 0u64;
    let mut unmeth_cpg = 0u64;
    let mut meth_chg = 0u64;
    let mut unmeth_chg = 0u64;
    let mut meth_chh = 0u64;
    let mut unmeth_chh = 0u64;

    let text = fs::read_to_string(&report_path)
        .with_context(|| format!("reading {}", report_path.display()))?;

    for line in text.lines() {
        let line = line.trim_end_matches('\r');
        if paired_end {
            if let Some(v) = tab_u64(line, "Sequence pairs analysed in total:") { total_reads = v; }
            if let Some(v) = tab_u64(line, "Sequence pairs with no alignments under any condition:") { unaligned = v; }
            if let Some(v) = tab_u64(line, "Sequence pairs did not map uniquely:") { ambig_reads = v; }
            if let Some(v) = tab_u64(line, "Sequence pairs which were discarded because genomic sequence could not be extracted:") { no_seq_reads = v; }
            if let Some(v) = tab_u64(line, "Number of paired-end alignments with a unique best hit:") { aligned_reads = v; }
        } else {
            if let Some(v) = tab_u64(line, "Sequences analysed in total:") { total_reads = v; }
            if let Some(v) = tab_u64(line, "Sequences with no alignments under any condition:") { unaligned = v; }
            if let Some(v) = tab_u64(line, "Sequences did not map uniquely:") { ambig_reads = v; }
            if let Some(v) = tab_u64(line, "Sequences which were discarded because genomic sequence could not be extracted:") { no_seq_reads = v; }
            if let Some(v) = tab_u64(line, "Number of alignments with a unique best hit from the different alignments:") { aligned_reads = v; }
        }
        if let Some(v) = tab_u64(line, "Total number of C's analysed:") { total_c = v; }
        if let Some(v) = tab_u64(line, "Total methylated C's in CpG context:") { meth_cpg = v; }
        if let Some(v) = tab_u64(line, "Total methylated C's in CHG context:") { meth_chg = v; }
        if let Some(v) = tab_u64(line, "Total methylated C's in CHH context:") { meth_chh = v; }
        if let Some(v) = tab_u64(line, "Total unmethylated C's in CpG context:") { unmeth_cpg = v; }
        if let Some(v) = tab_u64(line, "Total unmethylated C's in CHG context:") { unmeth_chg = v; }
        if let Some(v) = tab_u64(line, "Total unmethylated C's in CHH context:") { unmeth_chh = v; }
    }

    // Dedup report (optional)
    let dedup_path = if paired_end {
        PathBuf::from(format!("{base_stripped}_pe.deduplication_report.txt"))
    } else {
        PathBuf::from(format!("{base_stripped}.deduplication_report.txt"))
    };

    let mut dup_reads = 0u64;
    let mut unique_reads = 0u64;
    let has_dedup = dedup_path.exists();

    if has_dedup {
        let text2 = fs::read_to_string(&dedup_path)?;
        for line in text2.lines() {
            let line = line.trim_end_matches('\r');
            if line.starts_with("Total number of alignments analysed in ") {
                if let Some(v) = line.split('\t').nth(1).and_then(|s| s.parse().ok()) {
                    aligned_reads = v;
                }
            }
            if line.starts_with("Total number duplicated alignments removed:") {
                if let Some(v) = line.split('\t').nth(1).and_then(|s| s.split_whitespace().next()).and_then(|s| s.parse().ok()) {
                    dup_reads = v;
                }
            }
            if line.starts_with("Total count of deduplicated leftover sequences:") {
                if let Some(v) = line.split('\t').nth(1).and_then(|s| s.split_whitespace().next()).and_then(|s| s.parse().ok()) {
                    unique_reads = v;
                }
            }
        }
    } else {
        eprintln!("No deduplication report present, skipping...");
    }

    // Splitting report (optional)
    let split_path = if paired_end {
        if has_dedup { PathBuf::from(format!("{base_stripped}_pe.deduplicated_splitting_report.txt")) }
        else         { PathBuf::from(format!("{base_stripped}_pe_splitting_report.txt")) }
    } else if has_dedup {
        PathBuf::from(format!("{base_stripped}.deduplicated_splitting_report.txt"))
    } else {
        PathBuf::from(format!("{base_stripped}_splitting_report.txt"))
    };

    if split_path.exists() {
        let text3 = fs::read_to_string(&split_path)?;
        for line in text3.lines() {
            let line = line.trim_end_matches('\r');
            if let Some(v) = tab_u64(line, "Total number of C's analysed:") { total_c = v; }
            if let Some(v) = tab_u64(line, "Total methylated C's in CpG context:") { meth_cpg = v; }
            if let Some(v) = tab_u64(line, "Total methylated C's in CHG context:") { meth_chg = v; }
            if let Some(v) = tab_u64(line, "Total methylated C's in CHH context:") { meth_chh = v; }
            if let Some(v) = tab_u64(line, "Total C to T conversions in CpG context:") { unmeth_cpg = v; }
            if let Some(v) = tab_u64(line, "Total C to T conversions in CHG context:") { unmeth_chg = v; }
            if let Some(v) = tab_u64(line, "Total C to T conversions in CHH context:") { unmeth_chh = v; }
        }
    } else {
        eprintln!("No methylation extractor report present, skipping...");
    }

    let mut name = bam.file_stem().and_then(|s| s.to_str()).unwrap_or("sample").to_string();
    for suffix in &["_bismark", ".fq.gz", "_trimmed", "_1", "_2"] {
        if name.ends_with(suffix) { name.truncate(name.len() - suffix.len()); }
    }

    if verbose {
        eprintln!("  {name}: total={total_reads} aligned={aligned_reads} meth_cpg={meth_cpg}");
    }

    Ok(Sample {
        name, total_reads, aligned_reads, unaligned, ambig_reads, no_seq_reads,
        dup_reads, unique_reads, total_c,
        meth_cpg, unmeth_cpg, meth_chg, unmeth_chg, meth_chh, unmeth_chh,
    })
}

fn tab_u64(line: &str, prefix: &str) -> Option<u64> {
    if !line.starts_with(prefix) { return None; }
    line.split('\t').nth(1)?.split_whitespace().next()?.parse().ok()
}

// ─── HTML generation ─────────────────────────────────────────────────────────

fn build_html(samples: &[Sample], title: &str, _verbose: bool) -> Result<String> {
    let mut doc = SUMMARY_TEMPLATE.to_string();

    doc = replace_section(&doc, "plotly_goes_here",       PLOT_LY);
    doc = replace_section(&doc, "bismark_logo_goes_here", BISMARK_LOGO);
    doc = remove_section(&doc,  "bioinf_logo_goes_here");

    sub(&mut doc, "report_timestamp", &current_timestamp());
    sub(&mut doc, "page_title",       title);
    sub(&mut doc, "bismark_version",  BISMARK_VERSION);

    let n = samples.len();
    sub(&mut doc, "num_samples", &n.to_string());

    let x_str: String = (1..=n).map(|i| i.to_string()).collect::<Vec<_>>().join(",");
    sub(&mut doc, "x_values_alignment",   &x_str);
    sub(&mut doc, "x_values_methylation", &x_str);

    let categories: Vec<String> = samples.iter().map(|s| format!("'{}'", s.name)).collect();
    sub(&mut doc, "filenames_replace", &categories.join(","));

    let has_dedup = samples.iter().any(|s| s.dup_reads > 0 || s.unique_reads > 0);

    let joined = |vals: &[u64]| -> String { vals.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",") };
    let only_empty = |s: &str| s.is_empty() || s.chars().all(|c| c == ',');

    let aligned_str = if has_dedup {
        String::new()
    } else {
        joined(&samples.iter().map(|s| s.aligned_reads).collect::<Vec<_>>())
    };

    sub(&mut doc, "aligned_seq",       &aligned_str);
    sub(&mut doc, "dup_alignments",    &joined(&samples.iter().map(|s| s.dup_reads).collect::<Vec<_>>()));
    sub(&mut doc, "unique_alignments", &joined(&samples.iter().map(|s| s.unique_reads).collect::<Vec<_>>()));
    sub(&mut doc, "not_aligned",       &joined(&samples.iter().map(|s| s.unaligned).collect::<Vec<_>>()));
    sub(&mut doc, "ambig_aligned",     &joined(&samples.iter().map(|s| s.ambig_reads).collect::<Vec<_>>()));
    sub(&mut doc, "no_seq",            &joined(&samples.iter().map(|s| s.no_seq_reads).collect::<Vec<_>>()));

    let dup_str = joined(&samples.iter().map(|s| s.dup_reads).collect::<Vec<_>>());
    if only_empty(&dup_str) {
        doc = remove_section(&doc, "deduplicated_unique_reads_section");
        doc = remove_section(&doc, "duplicated_reads_section");
        doc = doc.replace("{{raw_aligned_reads_section}}", "");
    } else {
        doc = remove_section(&doc, "raw_aligned_reads_section");
        doc = doc.replace("{{deduplicated_unique_reads_section}}", "");
        doc = doc.replace("{{duplicated_reads_section}}", "");
    }

    sub(&mut doc, "meth_cpg_string",   &joined(&samples.iter().map(|s| s.meth_cpg).collect::<Vec<_>>()));
    sub(&mut doc, "unmeth_cpg_string", &joined(&samples.iter().map(|s| s.unmeth_cpg).collect::<Vec<_>>()));
    sub(&mut doc, "meth_chg_string",   &joined(&samples.iter().map(|s| s.meth_chg).collect::<Vec<_>>()));
    sub(&mut doc, "unmeth_chg_string", &joined(&samples.iter().map(|s| s.unmeth_chg).collect::<Vec<_>>()));
    sub(&mut doc, "meth_chh_string",   &joined(&samples.iter().map(|s| s.meth_chh).collect::<Vec<_>>()));
    sub(&mut doc, "unmeth_chh_string", &joined(&samples.iter().map(|s| s.unmeth_chh).collect::<Vec<_>>()));

    // Alignment percentages
    let mut p_al: Vec<String>    = Vec::new();
    let mut p_dedup: Vec<String> = Vec::new();
    let mut p_dup: Vec<String>   = Vec::new();
    let mut p_unal: Vec<String>  = Vec::new();
    let mut p_noseq: Vec<String> = Vec::new();
    let mut p_ambig: Vec<String> = Vec::new();

    for s in samples {
        let total = if has_dedup {
            (s.unique_reads + s.dup_reads + s.no_seq_reads + s.unaligned + s.ambig_reads) as f64
        } else {
            (s.aligned_reads + s.no_seq_reads + s.unaligned + s.ambig_reads) as f64
        };
        let fmt = |n: u64| -> String {
            if total == 0.0 { "0.00".into() } else { format!("{:.2}", n as f64 / total * 100.0) }
        };
        if has_dedup { p_dedup.push(fmt(s.unique_reads)); p_dup.push(fmt(s.dup_reads)); }
        else         { p_al.push(fmt(s.aligned_reads)); }
        p_unal.push(fmt(s.unaligned));
        p_noseq.push(fmt(s.no_seq_reads));
        p_ambig.push(fmt(s.ambig_reads));
    }

    if has_dedup {
        doc = remove_section(&doc, "raw_unique_reads_percentage_section");
        doc = doc.replace("{{duplicated_reads_percentage_section}}", "");
        doc = doc.replace("{{deduplicated_unique_reads_percentage_section}}", "");
        sub(&mut doc, "p_deduplicated_unique_alignments", &p_dedup.join(","));
        sub(&mut doc, "p_duplicated_alignments",          &p_dup.join(","));
    } else {
        doc = remove_section(&doc, "deduplicated_unique_reads_percentage_section");
        doc = remove_section(&doc, "duplicated_reads_percentage_section");
        doc = doc.replace("{{raw_unique_reads_percentage_section}}", "");
        sub(&mut doc, "p_aligned_replace", &p_al.join(","));
    }

    sub(&mut doc, "p_no_seq_replace", &p_noseq.join(","));
    sub(&mut doc, "p_unal_replace",   &p_unal.join(","));
    sub(&mut doc, "p_ambig_replace",  &p_ambig.join(","));

    // Methylation percentages
    let meth_pcts = |m_arr: &[u64], u_arr: &[u64], zero_on_empty: bool| -> (String, String) {
        let mut pm: Vec<String> = Vec::new();
        let mut pu: Vec<String> = Vec::new();
        for (&m, &u) in m_arr.iter().zip(u_arr.iter()) {
            let t = (m + u) as f64;
            if t == 0.0 {
                let na = if zero_on_empty { "0" } else { "NA" };
                pm.push(na.into()); pu.push(na.into());
            } else {
                let pmc = format!("{:.2}", m as f64 / t * 100.0);
                let puc = format!("{:.2}", u as f64 / t * 100.0);
                pm.push(pmc); pu.push(puc);
            }
        }
        (pm.join(","), pu.join(","))
    };

    let mc: Vec<u64> = samples.iter().map(|s| s.meth_cpg).collect();
    let uc: Vec<u64> = samples.iter().map(|s| s.unmeth_cpg).collect();
    let (pm, pu) = meth_pcts(&mc, &uc, false);
    sub(&mut doc, "p_CpG_m_replace", &pm);
    sub(&mut doc, "p_CpG_u_replace", &pu);

    let mc: Vec<u64> = samples.iter().map(|s| s.meth_chg).collect();
    let uc: Vec<u64> = samples.iter().map(|s| s.unmeth_chg).collect();
    let (pm, pu) = meth_pcts(&mc, &uc, true);
    sub(&mut doc, "p_CHG_m_replace", &pm);
    sub(&mut doc, "p_CHG_u_replace", &pu);

    let mc: Vec<u64> = samples.iter().map(|s| s.meth_chh).collect();
    let uc: Vec<u64> = samples.iter().map(|s| s.unmeth_chh).collect();
    let (pm, pu) = meth_pcts(&mc, &uc, true);
    sub(&mut doc, "p_CHH_m_replace", &pm);
    sub(&mut doc, "p_CHH_u_replace", &pu);

    Ok(doc)
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

fn replace_section(doc: &str, tag: &str, content: &str) -> String {
    let marker = format!("{{{{{tag}}}}}");
    if let Some(start) = doc.find(&marker) {
        if let Some(end_rel) = doc[start + marker.len()..].find(&marker) {
            let end = start + marker.len() + end_rel;
            let mut result = String::with_capacity(doc.len() + content.len());
            result.push_str(&doc[..start]);
            result.push_str(content);
            result.push_str(&doc[end + marker.len()..]);
            return result;
        }
    }
    doc.to_string()
}

fn remove_section(doc: &str, tag: &str) -> String {
    let marker = format!("{{{{{tag}}}}}");
    if let Some(start) = doc.find(&marker) {
        if let Some(end_rel) = doc[start + marker.len()..].find(&marker) {
            let end = start + marker.len() + end_rel;
            let mut result = String::with_capacity(doc.len());
            result.push_str(&doc[..start]);
            result.push_str(&doc[end + marker.len()..]);
            return result;
        }
    }
    doc.replace(&marker, "")
}

fn sub(doc: &mut String, tag: &str, value: &str) {
    let key = format!("{{{{{tag}}}}}");
    *doc = doc.replace(&key, value);
}

fn current_timestamp() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now().duration_since(UNIX_EPOCH).unwrap_or_default().as_secs();
    let s = secs % 60;
    let m = (secs / 60) % 60;
    let h = (secs / 3600) % 24;
    let days = secs / 86400;
    let (yr, mo, day) = days_to_date(days);
    format!("{yr:04}-{mo:02}-{day:02} {h:02}:{m:02}:{s:02} UTC")
}

fn days_to_date(mut days: u64) -> (u64, u64, u64) {
    let mut year = 1970u64;
    loop {
        let ydays = if is_leap(year) { 366 } else { 365 };
        if days < ydays { break; }
        days -= ydays;
        year += 1;
    }
    let month_days = [31u64, if is_leap(year) { 29 } else { 28 }, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    let mut month = 1u64;
    for &md in &month_days {
        if days < md { break; }
        days -= md;
        month += 1;
    }
    (year, month, days + 1)
}

fn is_leap(y: u64) -> bool { (y % 4 == 0 && y % 100 != 0) || y % 400 == 0 }

fn auto_detect_bam_files() -> Result<Vec<PathBuf>> {
    let cwd = std::env::current_dir()?;
    let suffixes = ["bismark_bt2.bam", "bismark_bt2_pe.bam", "bismark_hisat2.bam", "bismark_hisat2_pe.bam"];
    let mut results = Vec::new();
    for suffix in &suffixes {
        let mut found: Vec<PathBuf> = Vec::new();
        for entry in fs::read_dir(&cwd)? {
            let e = entry?;
            if e.file_name().to_string_lossy().ends_with(suffix) {
                found.push(e.path());
            }
        }
        if !found.is_empty() {
            eprintln!("Found {} file(s) matching *{suffix}", found.len());
            found.sort();
            results.extend(found);
        }
    }
    Ok(results)
}
