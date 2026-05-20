// Integration tests for coverage2cytosine.
// Creates a minimal genome FASTA and coverage file, then checks output format.

use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_coverage2cytosine"))
}

/// Minimal FASTA: chromosome "chr1" with known CpG sites.
/// Sequence: ACGTCGACGT (10 bp) — CpG at positions 4-5 (1-based: C at 5, G at 6).
fn write_genome(dir: &std::path::Path) {
    let fasta = dir.join("genome.fa");
    let mut f = fs::File::create(fasta).unwrap();
    writeln!(f, ">chr1").unwrap();
    writeln!(f, "ACGTCGACGT").unwrap();
}

/// Coverage file format: chr, start, end, pct_meth, count_meth, count_unmeth (1-based).
/// Covering position 5 (the C in CpG): 80% methylated, 4 meth reads, 1 unmeth read.
fn write_coverage(path: &PathBuf) {
    let mut f = fs::File::create(path).unwrap();
    writeln!(f, "chr1\t5\t5\t80.0\t4\t1").unwrap();
}

// ─── --version ───────────────────────────────────────────────────────────────

#[test]
fn test_version() {
    let out = Command::new(bin()).arg("--version").output().expect("binary not found");
    let combined = format!(
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(
        combined.contains("0.25.1") || combined.contains("coverage2cytosine") || combined.contains("Bismark"),
        "expected version string, got: {combined}"
    );
}

// ─── Basic CpG report generation ─────────────────────────────────────────────

#[test]
fn test_basic_cpg_report() {
    let dir = tempfile::tempdir().unwrap();
    write_genome(dir.path());
    let cov = dir.path().join("sample.cov");
    write_coverage(&cov);
    let out_report = dir.path().join("output.CpG_report.txt");

    let output = Command::new(bin())
        .arg("--genome_folder").arg(dir.path())
        .arg("--output").arg(&out_report)
        .arg(&cov)
        .output()
        .expect("binary failed");

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    assert!(out_report.exists(), "CpG report not created at {:?}", out_report);

    let content = fs::read_to_string(&out_report).unwrap();
    // Must have some content with chr1
    assert!(!content.is_empty(), "output is empty");
    assert!(content.contains("chr1"), "chr1 not in output:\n{content}");
    // Position 5 should be in the report (the covered CpG C)
    assert!(content.contains("5"), "position 5 not in report:\n{content}");
    // Methylation counts: 4 meth, 1 unmeth
    assert!(content.contains("4"), "meth count 4 not in report:\n{content}");
    assert!(content.contains("1"), "unmeth count 1 not in report:\n{content}");
    // Context should be CG
    assert!(content.contains("CG"), "CG context not in report:\n{content}");
}

// ─── Output includes uncovered cytosines ─────────────────────────────────────

#[test]
fn test_uncovered_cytosines_in_report() {
    let dir = tempfile::tempdir().unwrap();
    write_genome(dir.path());
    let cov = dir.path().join("sample.cov");
    write_coverage(&cov);
    let out_report = dir.path().join("out.CpG_report.txt");

    let output = Command::new(bin())
        .arg("--genome_folder").arg(dir.path())
        .arg("--output").arg(&out_report)
        .arg(&cov)
        .output()
        .expect("binary failed");

    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));

    let content = fs::read_to_string(&out_report).unwrap();
    // Every CpG in the genome should appear, even with 0 coverage.
    // The genome ACGTCGACGT has one CpG (C at pos 5, G at pos 6).
    // The report should have exactly the positions covered (pos 5, 6).
    let cpg_lines: Vec<&str> = content.lines()
        .filter(|l| l.contains("CG") || l.split('\t').nth(5).map(|c| c == "CG").unwrap_or(false))
        .collect();
    assert!(!cpg_lines.is_empty(), "no CG lines in report:\n{content}");
}
