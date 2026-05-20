use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use bismark_lib::BISMARK_VERSION;
use clap::Parser;

// Embedded templates (compile-time paths relative to src/main.rs)
const PLOTLY_TEMPLATE: &str = include_str!("../../../plotly/plotly_template.tpl");
const PLOT_LY:         &str = include_str!("../../../plotly/plot.ly");
const BISMARK_LOGO:    &str = include_str!("../../../plotly/bismark.logo");
const BIOINF_LOGO:     &str = include_str!("../../../plotly/bioinf.logo");

#[derive(Parser)]
#[command(
    name = "bismark2report",
    about = "Generate graphical HTML report from Bismark reports",
    version = BISMARK_VERSION,
    disable_version_flag = true,
)]
struct Cli {
    #[arg(long = "alignment_report")]
    alignment_report: Option<PathBuf>,

    #[arg(long = "dedup_report")]
    dedup_report: Option<String>,

    #[arg(long = "splitting_report")]
    splitting_report: Option<String>,

    #[arg(long = "mbias_report")]
    mbias_report: Option<String>,

    #[arg(long = "nucleotide_report")]
    nucleotide_report: Option<String>,

    #[arg(long = "dir", default_value = "")]
    dir: String,

    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    #[arg(long = "verbose")]
    verbose: bool,

    #[arg(long = "version")]
    version: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.version {
        println!(
            "\n\n                    Bismark HTML Report Module\n\n                    bismark2report version: {}\n                Copyright 2010-25 Felix Krueger, Altos Bioinformatics\n                       https://github.com/FelixKrueger/Bismark\n\n",
            BISMARK_VERSION
        );
        return Ok(());
    }

    let output_dir = normalise_dir(&cli.dir);

    // Collect all alignment reports to process
    let alignment_reports: Vec<PathBuf> = if let Some(ref path) = cli.alignment_report {
        vec![path.clone()]
    } else {
        glob_reports("*E_report.txt")?
    };

    if alignment_reports.is_empty() {
        bail!("Found no alignment reports in current directory. Specify one with --alignment_report");
    }

    if alignment_reports.len() > 1 && cli.output.is_some() {
        bail!("Cannot use -o/--output with more than one alignment report");
    }

    eprintln!("Found {} alignment report(s)", alignment_reports.len());

    for aln_path in &alignment_reports {
        let basename = extract_basename(aln_path);

        let dedup_path  = resolve_optional(&cli.dedup_report,      &basename, "*deduplication_report.txt")?;
        let split_path  = resolve_optional(&cli.splitting_report,   &basename, "*splitting_report.txt")?;
        let mbias_path  = resolve_optional(&cli.mbias_report,       &basename, "*M-bias.txt")?;
        let nuc_path    = resolve_optional(&cli.nucleotide_report,  &basename, "*nucleotide_stats.txt")?;

        // Output filename
        let html_out = if let Some(ref out) = cli.output {
            format!("{}{}", output_dir, out.to_string_lossy())
        } else {
            let stem = aln_path.file_stem().and_then(|s| s.to_str()).unwrap_or("bismark_report");
            format!("{}{}.html", output_dir, stem)
        };

        eprintln!("\nWriting Bismark HTML report to >> {} <<\n", html_out);

        let doc = build_report(aln_path, dedup_path.as_deref(), split_path.as_deref(),
                               mbias_path.as_deref(), nuc_path.as_deref(), cli.verbose)?;

        let mut f = fs::File::create(&html_out)
            .with_context(|| format!("creating {html_out}"))?;
        f.write_all(doc.as_bytes())?;
    }

    Ok(())
}

fn build_report(
    aln_path:   &Path,
    dedup_path: Option<&Path>,
    split_path: Option<&Path>,
    mbias_path: Option<&Path>,
    nuc_path:   Option<&Path>,
    verbose:    bool,
) -> Result<String> {
    let mut doc = PLOTLY_TEMPLATE.to_string();

    // Inject plot.ly
    doc = replace_section(&doc, "plotly_goes_here", PLOT_LY);
    // Inject logos
    doc = replace_section(&doc, "bismark_logo_goes_here", BISMARK_LOGO);
    doc = replace_section(&doc, "bioinf_logo_goes_here",  BIOINF_LOGO);

    // Timestamp
    let now = chrono_now();
    doc = doc.replace("{{date}}", &now.0);
    doc = doc.replace("{{time}}", &now.1);

    // Alignment report (mandatory)
    doc = read_alignment_report(aln_path, doc, verbose)?;

    // Deduplication (optional)
    if let Some(p) = dedup_path {
        eprintln!("Using deduplication report: {}", p.display());
        doc = doc.replace("{{deduplication_section}}", "").replace("{{deduplication_section}}", "");
        doc = read_dedup_report(p, doc, verbose)?;
    } else {
        eprintln!("No deduplication report, skipping");
        doc = remove_section(&doc, "deduplication_section");
    }

    // Splitting report (optional)
    if let Some(p) = split_path {
        eprintln!("Using splitting report: {}", p.display());
        doc = doc.replace("{{cytosine_methylation_post_deduplication_section}}", "");
        doc = read_splitting_report(p, doc, verbose)?;
    } else {
        eprintln!("No splitting report, skipping");
        doc = remove_section(&doc, "cytosine_methylation_post_deduplication_section");
    }

    // M-bias report (optional)
    if let Some(p) = mbias_path {
        eprintln!("Using M-bias report: {}", p.display());
        doc = doc.replace("{{mbias_r1_section}}", "");
        let (state, new_doc) = read_mbias_report(p, doc, verbose)?;
        doc = new_doc;
        if state == "single" {
            doc = remove_section(&doc, "mbias_r2_section");
        } else {
            doc = doc.replace("{{mbias_r2_section}}", "");
        }
    } else {
        eprintln!("No M-bias report, skipping");
        doc = remove_section(&doc, "mbias_r1_section");
        doc = remove_section(&doc, "mbias_r2_section");
    }

    // Nucleotide coverage (optional)
    if let Some(p) = nuc_path {
        eprintln!("Using nucleotide coverage report: {}", p.display());
        doc = doc.replace("{{nucleotide_coverage_section}}", "");
        doc = read_nuc_report(p, doc, verbose)?;
    } else {
        eprintln!("No nucleotide coverage report, skipping");
        doc = remove_section(&doc, "nucleotide_coverage_section");
    }

    Ok(doc)
}

// ─── Template helpers ────────────────────────────────────────────────────────

fn replace_section(doc: &str, tag: &str, content: &str) -> String {
    let open  = format!("{{{{{tag}}}}}");
    let close = format!("{{{{{tag}}}}}");
    if let Some(start) = doc.find(&open) {
        if let Some(end_rel) = doc[start + open.len()..].find(&close) {
            let end = start + open.len() + end_rel;
            let mut result = String::with_capacity(doc.len() + content.len());
            result.push_str(&doc[..start]);
            result.push_str(content);
            result.push_str(&doc[end + close.len()..]);
            return result;
        }
    }
    doc.to_string()
}

fn remove_section(doc: &str, tag: &str) -> String {
    let open  = format!("{{{{{tag}}}}}");
    let close = format!("{{{{{tag}}}}}");
    if let Some(start) = doc.find(&open) {
        if let Some(end_rel) = doc[start + open.len()..].find(&close) {
            let end = start + open.len() + end_rel;
            let mut result = String::with_capacity(doc.len());
            result.push_str(&doc[..start]);
            result.push_str(&doc[end + close.len()..]);
            return result;
        }
    }
    // Tag not found or no pair — strip any remaining single occurrences
    doc.replace(&open, "")
}

fn sub(doc: &mut String, tag: &str, value: &str) {
    let key = format!("{{{{{tag}}}}}");
    *doc = doc.replace(&key, value);
}

// ─── Report parsers ──────────────────────────────────────────────────────────

fn read_alignment_report(path: &Path, mut doc: String, verbose: bool) -> Result<String> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("reading {}", path.display()))?;

    let mut vals: HashMap<&str, String> = HashMap::new();
    let mut total_seq_text = String::new();
    let mut unique_text    = String::new();
    let mut no_aln_text    = String::new();
    let mut multiple_text  = String::new();
    let mut meth_unknown: Option<String> = None;
    let mut unmeth_unknown: Option<String> = None;
    let mut perc_unknown: Option<String>  = None;

    for line in text.lines() {
        let line = line.trim_end_matches('\r');
        let tab_val = || line.split('\t').nth(1).unwrap_or("").to_string();

        if line.starts_with("Sequence pairs analysed in total:") {
            vals.insert("total_sequences_alignments", tab_val());
            total_seq_text = "Sequence pairs analysed in total".into();
        } else if line.starts_with("Sequences analysed in total:") {
            vals.insert("total_sequences_alignments", tab_val());
            total_seq_text = "Sequences analysed in total".into();
        } else if line.starts_with("Bismark report for:") {
            // "Bismark report for: FILE (version: VER)"
            if let Some(rest) = line.strip_prefix("Bismark report for: ") {
                if let Some(idx) = rest.rfind(" (version: ") {
                    vals.insert("filename", rest[..idx].to_string());
                    let ver = rest[idx + " (version: ".len()..].trim_end_matches(')');
                    vals.insert("bismark_version", ver.to_string());
                }
            }
        } else if line.starts_with("Number of paired-end alignments with a unique best hit:") {
            vals.insert("unique_seqs", tab_val());
            unique_text = "Paired-end alignments with a unique best hit".into();
        } else if line.starts_with("Number of alignments with a unique best hit from") {
            vals.insert("unique_seqs", tab_val());
            unique_text = "Single-end alignments with a unique best hit".into();
        } else if line.starts_with("Sequence pairs with no alignments under any condition:") {
            vals.insert("no_alignments", tab_val());
            no_aln_text = "Pairs without alignments under any condition".into();
        } else if line.starts_with("Sequences with no alignments under any condition:") {
            vals.insert("no_alignments", tab_val());
            no_aln_text = "Sequences without alignments under any condition".into();
        } else if line.starts_with("Sequence pairs did not map uniquely:") {
            vals.insert("multiple_alignments", tab_val());
            multiple_text = "Pairs that did not map uniquely".into();
        } else if line.starts_with("Sequences did not map uniquely:") {
            vals.insert("multiple_alignments", tab_val());
            multiple_text = "Sequences that did not map uniquely".into();
        } else if line.starts_with("Sequence pairs which were discarded because genomic sequence could not be extracted:")
               || line.starts_with("Sequences which were discarded because genomic sequence could not be extracted:") {
            vals.insert("no_genomic", tab_val());
        } else if line.starts_with("Total number of C") {
            vals.insert("total_C_count", tab_val());
        } else if line.starts_with("Total methylated C's in CpG context:") {
            vals.insert("meth_CpG", tab_val());
        } else if line.starts_with("Total methylated C's in CHG context:") {
            vals.insert("meth_CHG", tab_val());
        } else if line.starts_with("Total methylated C's in CHH context:") {
            vals.insert("meth_CHH", tab_val());
        } else if line.starts_with("Total methylated C's in Unknown context:") {
            meth_unknown = Some(tab_val());
        } else if line.starts_with("Total unmethylated C's in CpG context:")
               || line.starts_with("Total C to T conversions in CpG context:") {
            vals.insert("unmeth_CpG", tab_val());
        } else if line.starts_with("Total unmethylated C's in CHG context:")
               || line.starts_with("Total C to T conversions in CHG context:") {
            vals.insert("unmeth_CHG", tab_val());
        } else if line.starts_with("Total unmethylated C's in CHH context:")
               || line.starts_with("Total C to T conversions in CHH context:") {
            vals.insert("unmeth_CHH", tab_val());
        } else if line.starts_with("Total unmethylated C's in Unknown context:")
               || line.starts_with("Total C to T conversions in Unknown context:") {
            unmeth_unknown = Some(tab_val());
        } else if line.starts_with("C methylated in CpG context:") {
            vals.insert("perc_CpG", tab_val().trim_end_matches('%').to_string());
        } else if line.starts_with("C methylated in CHG context:") {
            vals.insert("perc_CHG", tab_val().trim_end_matches('%').to_string());
        } else if line.starts_with("C methylated in CHH context:") {
            vals.insert("perc_CHH", tab_val().trim_end_matches('%').to_string());
        } else if line.starts_with("C methylated in Unknown context") {
            perc_unknown = Some(tab_val().trim_end_matches('%').to_string());
        } else if line.starts_with("CT/GA/CT:") {
            vals.insert("number_OT", tab_val());
        } else if line.starts_with("CT/CT:") {
            vals.insert("number_OT", tab_val());
        } else if line.starts_with("GA/CT/CT:") {
            vals.insert("number_CTOT", tab_val());
        } else if line.starts_with("GA/CT:") {
            vals.insert("number_CTOT", tab_val());
        } else if line.starts_with("GA/CT/GA:") {
            vals.insert("number_CTOB", tab_val());
        } else if line.starts_with("GA/GA:") {
            vals.insert("number_CTOB", tab_val());
        } else if line.starts_with("CT/GA/GA:") {
            vals.insert("number_OB", tab_val());
        } else if line.starts_with("CT/GA:") {
            vals.insert("number_OB", tab_val());
        }
    }

    if verbose {
        eprintln!("Alignment report: {:?}", vals);
    }

    sub(&mut doc, "sequences_analysed_in_total", &total_seq_text);
    sub(&mut doc, "unique_seqs_text",            &unique_text);
    sub(&mut doc, "no_alignments_text",          &no_aln_text);
    sub(&mut doc, "multiple_alignments_text",    &multiple_text);

    for (k, v) in &vals {
        sub(&mut doc, k, v);
    }

    // Plotly arrays
    let unique   = vals.get("unique_seqs").map(|s| s.as_str()).unwrap_or("0");
    let no_aln   = vals.get("no_alignments").map(|s| s.as_str()).unwrap_or("0");
    let multiple = vals.get("multiple_alignments").map(|s| s.as_str()).unwrap_or("0");
    let no_geno  = vals.get("no_genomic").map(|s| s.as_str()).unwrap_or("0");
    sub(&mut doc, "alignment_stats_plotly", &format!("{unique},{no_aln},{multiple},{no_geno}"));

    let ot   = vals.get("number_OT").map(|s| s.as_str()).unwrap_or("0");
    let ctot = vals.get("number_CTOT").map(|s| s.as_str()).unwrap_or("0");
    let ctob = vals.get("number_CTOB").map(|s| s.as_str()).unwrap_or("0");
    let ob   = vals.get("number_OB").map(|s| s.as_str()).unwrap_or("0");
    sub(&mut doc, "strand_alignment_plotly", &format!("{ot},{ctot},{ctob},{ob}"));

    let perc_cpg_g = perc_graph(vals.get("perc_CpG").map(|s| s.as_str()).unwrap_or("N/A"));
    let perc_chg_g = perc_graph(vals.get("perc_CHG").map(|s| s.as_str()).unwrap_or("N/A"));
    let perc_chh_g = perc_graph(vals.get("perc_CHH").map(|s| s.as_str()).unwrap_or("N/A"));
    sub(&mut doc, "cytosine_methylation_plotly", &format!("{perc_cpg_g},{perc_chg_g},{perc_chh_g}"));

    // Unknown context (bowtie2 only)
    let (mu, uu, pu) = build_unknown_html(meth_unknown.as_deref(), unmeth_unknown.as_deref(), perc_unknown.as_deref());
    sub(&mut doc, "meth_unknown",   &mu);
    sub(&mut doc, "unmeth_unknown", &uu);
    sub(&mut doc, "perc_unknown",   &pu);

    // Fallback N/A for percentage fields
    for k in &["perc_CpG", "perc_CHG", "perc_CHH"] {
        if doc.contains(&format!("{{{{{k}}}}}")) {
            sub(&mut doc, k, "N/A");
        }
    }

    Ok(doc)
}

fn read_dedup_report(path: &Path, mut doc: String, _verbose: bool) -> Result<String> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("reading {}", path.display()))?;

    let mut total: Option<String> = None;
    let mut dups:  Option<String> = None;
    let mut diff:  Option<String> = None;
    let mut leftover: Option<String> = None;

    for line in text.lines() {
        let line = line.trim_end_matches('\r');
        let tab_val = || line.split('\t').nth(1).unwrap_or("").split_whitespace().next().unwrap_or("").to_string();

        if line.starts_with("Total number of alignments") {
            total = Some(line.split('\t').nth(1).unwrap_or("0").to_string());
        } else if line.starts_with("Total number duplicated") {
            dups = Some(tab_val());
        } else if line.starts_with("Duplicated alignments were found at") {
            diff = Some(tab_val());
        } else if let Some(rest) = line.strip_prefix("Total count of deduplicated leftover sequences: ") {
            leftover = Some(rest.split_whitespace().next().unwrap_or("0").to_string());
        }
    }

    // Compute leftover if not present
    if leftover.is_none() {
        if let (Some(ref t), Some(ref d)) = (&total, &dups) {
            let t_n: u64 = t.parse().unwrap_or(0);
            let d_n: u64 = d.parse().unwrap_or(0);
            leftover = Some((t_n.saturating_sub(d_n)).to_string());
        }
    }

    let total   = total.unwrap_or_default();
    let dups    = dups.unwrap_or_default();
    let diff    = diff.unwrap_or_default();
    let leftover = leftover.unwrap_or_default();

    sub(&mut doc, "seqs_total_duplicates",         &total);
    sub(&mut doc, "unique_alignments_duplicates",  &leftover);
    sub(&mut doc, "duplicate_alignments_duplicates", &dups);
    sub(&mut doc, "different_positions_duplicates", &diff);
    sub(&mut doc, "duplication_stats_plotly",       &format!("{leftover},{dups}"));

    Ok(doc)
}

fn read_splitting_report(path: &Path, mut doc: String, _verbose: bool) -> Result<String> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("reading {}", path.display()))?;

    let mut meth_cpg: Option<String> = None;
    let mut meth_chg: Option<String> = None;
    let mut meth_chh: Option<String> = None;
    let mut meth_unk: Option<String> = None;
    let mut unmeth_cpg: Option<String> = None;
    let mut unmeth_chg: Option<String> = None;
    let mut unmeth_chh: Option<String> = None;
    let mut unmeth_unk: Option<String> = None;
    let mut perc_cpg:  Option<String> = None;
    let mut perc_chg:  Option<String> = None;
    let mut perc_chh:  Option<String> = None;
    let mut perc_unk:  Option<String> = None;
    let mut total_c:   Option<String> = None;

    for line in text.lines() {
        let line = line.trim_end_matches('\r');
        let tab_val = || line.split('\t').nth(1).unwrap_or("").to_string();

        if line.starts_with("Total number of C") {
            total_c = Some(tab_val());
        } else if line.starts_with("Total methylated C's in CpG context:") {
            meth_cpg = Some(tab_val());
        } else if line.starts_with("Total methylated C's in CHG context:") {
            meth_chg = Some(tab_val());
        } else if line.starts_with("Total methylated C's in CHH context:") {
            meth_chh = Some(tab_val());
        } else if line.starts_with("Total methylated C's in Unknown context:") {
            meth_unk = Some(tab_val());
        } else if line.starts_with("Total C to T conversions in CpG context:")
               || line.starts_with("Total unmethylated C's in CpG context:") {
            unmeth_cpg = Some(tab_val());
        } else if line.starts_with("Total C to T conversions in CHG context:")
               || line.starts_with("Total unmethylated C's in CHG context:") {
            unmeth_chg = Some(tab_val());
        } else if line.starts_with("Total C to T conversions in CHH context:")
               || line.starts_with("Total unmethylated C's in CHH context:") {
            unmeth_chh = Some(tab_val());
        } else if line.starts_with("Total C to T conversions in Unknown context:")
               || line.starts_with("Total unmethylated C's in Unknown context:") {
            unmeth_unk = Some(tab_val());
        } else if line.starts_with("C methylated in CpG context:") {
            perc_cpg = Some(tab_val().trim_end_matches('%').to_string());
        } else if line.starts_with("C methylated in CHG context:") {
            perc_chg = Some(tab_val().trim_end_matches('%').to_string());
        } else if line.starts_with("C methylated in CHH context:") {
            perc_chh = Some(tab_val().trim_end_matches('%').to_string());
        } else if line.starts_with("C methylated in Unknown context:") {
            perc_unk = Some(tab_val().trim_end_matches('%').to_string());
        }
    }

    sub(&mut doc, "total_C_count_splitting", &total_c.unwrap_or_default());
    sub(&mut doc, "meth_CpG_splitting",   &meth_cpg.as_deref().unwrap_or("0").to_string());
    sub(&mut doc, "meth_CHG_splitting",   &meth_chg.as_deref().unwrap_or("0").to_string());
    sub(&mut doc, "meth_CHH_splitting",   &meth_chh.as_deref().unwrap_or("0").to_string());
    sub(&mut doc, "unmeth_CpG_splitting", &unmeth_cpg.as_deref().unwrap_or("0").to_string());
    sub(&mut doc, "unmeth_CHG_splitting", &unmeth_chg.as_deref().unwrap_or("0").to_string());
    sub(&mut doc, "unmeth_CHH_splitting", &unmeth_chh.as_deref().unwrap_or("0").to_string());

    let pcg = perc_cpg.as_deref().unwrap_or("N/A");
    let pchg = perc_chg.as_deref().unwrap_or("N/A");
    let pchh = perc_chh.as_deref().unwrap_or("N/A");
    let _punk = perc_unk.as_deref().unwrap_or("N/A");

    sub(&mut doc, "perc_CpG_splitting", pcg);
    sub(&mut doc, "perc_CHG_splitting", pchg);
    sub(&mut doc, "perc_CHH_splitting", pchh);

    sub(&mut doc, "cytosine_methylation_post_duplication_plotly",
        &format!("{},{},{}", perc_graph(pcg), perc_graph(pchg), perc_graph(pchh)));

    let (mu, uu, pu) = build_unknown_html_split(
        meth_unk.as_deref(), unmeth_unk.as_deref(), perc_unk.as_deref());
    sub(&mut doc, "meth_unknown_splitting",   &mu);
    sub(&mut doc, "unmeth_unknown_splitting", &uu);
    sub(&mut doc, "perc_unknown_splitting",   &pu);

    Ok(doc)
}

fn read_mbias_report(path: &Path, mut doc: String, _verbose: bool) -> Result<(String, String)> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("reading {}", path.display()))?;

    // mbias_1/2: context -> arrays of (pos, meth, unmeth, perc, coverage)
    let mut r1: HashMap<String, MbiasCtx> = HashMap::new();
    let mut r2: HashMap<String, MbiasCtx> = HashMap::new();

    let mut context = String::new();
    let mut read_id = 1u32;
    let mut state = "single".to_string();

    for line in text.lines() {
        let line = line.trim_end_matches('\r');
        if let Some(ctx_part) = parse_mbias_header(line) {
            context = ctx_part;
            if line.contains("R2") || line.contains("Read 2") {
                read_id = 2;
                state = "paired".to_string();
            } else {
                read_id = 1;
            }
        } else if line.starts_with(|c: char| c.is_ascii_digit()) {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 5 {
                let pos = parts[0].to_string();
                let _meth = parts[1];
                let _unmeth = parts[2];
                let perc = parts[3].to_string();
                let cov  = parts[4].to_string();
                let tbl = if read_id == 1 { &mut r1 } else { &mut r2 };
                let entry = tbl.entry(context.clone()).or_default();
                entry.cov_x.push(pos.clone());
                entry.cov_y.push(cov);
                entry.perc_x.push(pos);
                entry.perc_y.push(perc);
            }
        }
    }

    for ctx in &["CpG", "CHG", "CHH"] {
        inject_mbias(&mut doc, "mbias1", ctx, r1.get(*ctx).map(|e| e as &MbiasCtx).unwrap_or(&MbiasCtx::default()));
    }
    if !r2.is_empty() {
        for ctx in &["CpG", "CHG", "CHH"] {
            inject_mbias(&mut doc, "mbias2", ctx, r2.get(*ctx).map(|e| e as &MbiasCtx).unwrap_or(&MbiasCtx::default()));
        }
    }

    Ok((state, doc))
}

#[derive(Default)]
struct MbiasCtx {
    cov_x:  Vec<String>,
    cov_y:  Vec<String>,
    perc_x: Vec<String>,
    perc_y: Vec<String>,
}

fn parse_mbias_header(line: &str) -> Option<String> {
    if line.starts_with("CpG context") { return Some("CpG".into()); }
    if line.starts_with("CHG context") { return Some("CHG".into()); }
    if line.starts_with("CHH context") { return Some("CHH".into()); }
    None
}

fn inject_mbias(doc: &mut String, read: &str, ctx: &str, e: &MbiasCtx) {
    sub(doc, &format!("{read}_{ctx}_meth_x"),     &e.perc_x.join(","));
    sub(doc, &format!("{read}_{ctx}_meth_y"),     &e.perc_y.join(","));
    sub(doc, &format!("{read}_{ctx}_coverage_x"), &e.cov_x.join(","));
    sub(doc, &format!("{read}_{ctx}_coverage_y"), &e.cov_y.join(","));
}

fn read_nuc_report(path: &Path, mut doc: String, verbose: bool) -> Result<String> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("reading {}", path.display()))?;

    let mut nucs: HashMap<String, [String; 5]> = HashMap::new(); // [p_obs, p_exp, c_obs, c_exp, cov]
    let mut linecount = 0usize;

    for line in text.lines() {
        let line = line.trim_end_matches('\r');
        let parts: Vec<&str> = line.split('\t').collect();
        if linecount == 0 {
            linecount += 1;
            continue; // header
        }
        if parts.len() >= 6 {
            let elem = parts[0].to_string();
            nucs.insert(elem, [
                parts[2].to_string(), // p_obs
                parts[4].to_string(), // p_exp
                parts[1].to_string(), // c_obs
                parts[3].to_string(), // c_exp
                parts[5].to_string(), // coverage
            ]);
        }
        linecount += 1;
    }

    let nuc_order = ["A","T","C","G","AC","CA","TC","CT","CC","CG","GC","GG","AG","GA","TG","GT","TT","TA","AT","AA"];

    let mut y_arr:   Vec<String> = Vec::new();
    let mut x_samp:  Vec<String> = Vec::new();
    let mut x_geno:  Vec<String> = Vec::new();

    for nuc in &nuc_order {
        let entry = nucs.get(*nuc).cloned().unwrap_or_else(|| ["0","0","0","0","0"].map(|s| s.to_string()));
        let [p_obs, p_exp, c_obs, c_exp, cov] = &entry;

        if verbose {
            eprintln!("{nuc}: obs={p_obs} exp={p_exp}");
        }

        sub(&mut doc, &format!("nuc_{nuc}_p_obs"),     p_obs);
        sub(&mut doc, &format!("nuc_{nuc}_p_exp"),     p_exp);
        sub(&mut doc, &format!("nuc_{nuc}_counts_obs"), c_obs);
        sub(&mut doc, &format!("nuc_{nuc}_counts_exp"), c_exp);
        sub(&mut doc, &format!("nuc_{nuc}_coverage"),   cov);

        y_arr.push(format!("'{nuc}'"));
        x_samp.push(p_obs.clone());
        x_geno.push(p_exp.clone());
    }

    sub(&mut doc, "nucleo_sample_y",  &y_arr.join(","));
    sub(&mut doc, "nucleo_genomic_y", &y_arr.join(","));
    sub(&mut doc, "nucleo_sample_x",  &x_samp.join(" , "));
    sub(&mut doc, "nucleo_genomic_x", &x_geno.join(" , "));

    Ok(doc)
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

fn perc_graph(s: &str) -> &str {
    if s == "N/A" { "0" } else { s }
}

fn build_unknown_html(meth: Option<&str>, unmeth: Option<&str>, perc: Option<&str>) -> (String, String, String) {
    if let (Some(m), Some(u), Some(p)) = (meth, unmeth, perc) {
        let m_html = format!(
            "     <tr>\n                                <th>Methylated C's in Unknown context</th>\n    \t\t\t<td>{m}</td>\n    \t\t</tr>"
        );
        let u_html = format!(
            "     <tr>\n                                <th>Unmethylated C's in Unknown context</th>\n    \t\t\t<td>{u}</td>\n    \t\t</tr>"
        );
        let p_html = format!(
            "     <tr>\n                                <th>Methylated C's in Unknown context</th>\n    \t\t\t<td>{p}%</td>\n    \t\t</tr>"
        );
        (m_html, u_html, p_html)
    } else {
        (String::new(), String::new(), String::new())
    }
}

fn build_unknown_html_split(meth: Option<&str>, unmeth: Option<&str>, perc: Option<&str>) -> (String, String, String) {
    build_unknown_html(meth, unmeth, perc)
}

fn chrono_now() -> (String, String) {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now().duration_since(UNIX_EPOCH).unwrap_or_default().as_secs();
    // Simple UTC decomposition (no external crate needed)
    let s = secs % 60;
    let m = (secs / 60) % 60;
    let h = (secs / 3600) % 24;
    let days = secs / 86400;
    // Days since epoch to calendar date (Gregorian)
    let (year, month, day) = days_to_date(days);
    (
        format!("{year:04}-{month:02}-{day:02}"),
        format!("{h:02}:{m:02}:{s:02}"),
    )
}

fn days_to_date(mut days: u64) -> (u64, u64, u64) {
    let mut year = 1970u64;
    loop {
        let leap = is_leap(year);
        let ydays = if leap { 366 } else { 365 };
        if days < ydays { break; }
        days -= ydays;
        year += 1;
    }
    let leap = is_leap(year);
    let month_days = [31u64, if leap {29} else {28}, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    let mut month = 1u64;
    for &md in &month_days {
        if days < md { break; }
        days -= md;
        month += 1;
    }
    (year, month, days + 1)
}

fn is_leap(y: u64) -> bool { (y % 4 == 0 && y % 100 != 0) || y % 400 == 0 }

fn normalise_dir(s: &str) -> String {
    if s.is_empty() { return String::new(); }
    if s.ends_with('/') { s.to_string() } else { format!("{s}/") }
}

fn extract_basename(path: &Path) -> String {
    let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
    // Strip trailing _PE_report or _SE_report suffix
    if let Some(b) = stem.strip_suffix("_PE_report").or_else(|| stem.strip_suffix("_SE_report")) {
        b.to_string()
    } else {
        stem.to_string()
    }
}

fn glob_reports(pattern: &str) -> Result<Vec<PathBuf>> {
    let cwd = std::env::current_dir()?;
    let suffix = pattern.trim_start_matches('*');
    let mut results = Vec::new();
    for entry in fs::read_dir(&cwd).with_context(|| "reading current dir")? {
        let e = entry?;
        let name = e.file_name();
        let name_str = name.to_string_lossy();
        if name_str.ends_with(suffix) {
            results.push(e.path());
        }
    }
    results.sort();
    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ─── replace_section ─────────────────────────────────────────────────────

    #[test]
    fn test_replace_section_basic() {
        let doc = "before {{FOO}}old content{{FOO}} after";
        let result = replace_section(doc, "FOO", "NEW");
        assert_eq!(result, "before NEW after");
    }

    #[test]
    fn test_replace_section_empty_content() {
        let doc = "{{SEC}}stuff{{SEC}}";
        assert_eq!(replace_section(doc, "SEC", ""), "");
    }

    #[test]
    fn test_replace_section_tag_not_found() {
        let doc = "no tags here";
        assert_eq!(replace_section(doc, "MISSING", "x"), "no tags here");
    }

    // ─── remove_section ──────────────────────────────────────────────────────

    #[test]
    fn test_remove_section_basic() {
        let doc = "before {{BAR}}content to drop{{BAR}} after";
        assert_eq!(remove_section(doc, "BAR"), "before  after");
    }

    #[test]
    fn test_remove_section_missing_tag() {
        let doc = "nothing to remove";
        assert_eq!(remove_section(doc, "X"), "nothing to remove");
    }

    #[test]
    fn test_remove_section_single_orphan_tag() {
        // Only one occurrence — stripped by the fallback replace
        let doc = "text {{LONE}} more";
        assert_eq!(remove_section(doc, "LONE"), "text  more");
    }

    // ─── perc_graph ──────────────────────────────────────────────────────────

    #[test]
    fn test_perc_graph_na() {
        assert_eq!(perc_graph("N/A"), "0");
    }

    #[test]
    fn test_perc_graph_numeric() {
        assert_eq!(perc_graph("42.5"), "42.5");
    }

    #[test]
    fn test_perc_graph_zero() {
        assert_eq!(perc_graph("0"), "0");
    }

    // ─── extract_basename ────────────────────────────────────────────────────

    #[test]
    fn test_extract_basename_se_report() {
        let p = PathBuf::from("sample_SE_report.txt");
        assert_eq!(extract_basename(&p), "sample");
    }

    #[test]
    fn test_extract_basename_pe_report() {
        let p = PathBuf::from("sample_PE_report.txt");
        assert_eq!(extract_basename(&p), "sample");
    }

    #[test]
    fn test_extract_basename_no_suffix() {
        let p = PathBuf::from("mysample.txt");
        assert_eq!(extract_basename(&p), "mysample");
    }

    #[test]
    fn test_extract_basename_with_path() {
        let p = PathBuf::from("/data/run/sample_SE_report.txt");
        assert_eq!(extract_basename(&p), "sample");
    }

    // ─── days_to_date ────────────────────────────────────────────────────────

    #[test]
    fn test_days_to_date_epoch() {
        // Day 0 = 1970-01-01
        assert_eq!(days_to_date(0), (1970, 1, 1));
    }

    #[test]
    fn test_days_to_date_end_of_1970() {
        // Day 364 = 1970-12-31
        assert_eq!(days_to_date(364), (1970, 12, 31));
    }

    #[test]
    fn test_days_to_date_leap_year_day() {
        // 1972 is a leap year. 1970=365, 1971=365, so day 730 = 1972-01-01
        assert_eq!(days_to_date(730), (1972, 1, 1));
        // 1972-02-29 = day 730 + 59 = 789
        assert_eq!(days_to_date(789), (1972, 2, 29));
    }

    #[test]
    fn test_days_to_date_known_date() {
        // 2000-01-01: days from epoch
        // 1970..1999 = 30 years with 7 leap years (72,76,80,84,88,92,96) = 23*365 + 7*366 = 8395 + 2562 = 10957
        assert_eq!(days_to_date(10957), (2000, 1, 1));
    }
}

fn resolve_optional(
    cli_arg: &Option<String>,
    basename: &str,
    glob_suffix: &str,
) -> Result<Option<PathBuf>> {
    if let Some(ref s) = cli_arg {
        if s.to_lowercase() == "none" || s.is_empty() {
            return Ok(None);
        }
        return Ok(Some(PathBuf::from(s)));
    }

    // Auto-discover: look for files matching <basename><glob_suffix>
    // glob_suffix is like "*deduplication_report.txt" — strip the leading *
    let suffix = glob_suffix.trim_start_matches('*');
    let cwd = std::env::current_dir()?;
    let mut matches = Vec::new();
    for entry in fs::read_dir(&cwd).with_context(|| "reading current dir")? {
        let e = entry?;
        let name = e.file_name();
        let name_str = name.to_string_lossy();
        if name_str.starts_with(basename) && name_str.ends_with(suffix) {
            matches.push(e.path());
        }
    }

    match matches.len() {
        0 => Ok(None),
        1 => Ok(Some(matches.remove(0))),
        n => {
            eprintln!("Warning: found {n} matches for {basename}{suffix}, using first");
            Ok(Some(matches.remove(0)))
        }
    }
}
