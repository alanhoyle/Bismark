// Integration tests for bismark2bedGraph.

use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_bismark2bedGraph"))
}

/// Write a minimal bismark_methylation_extractor output file (CpG context, no header).
fn write_extractor_output(path: &PathBuf, records: &[&str]) {
    let mut f = fs::File::create(path).unwrap();
    for rec in records {
        writeln!(f, "{rec}").unwrap();
    }
}

// ─── --version (via --help which doesn't enforce required args) ───────────────

#[test]
fn test_version() {
    let out = Command::new(bin()).arg("--help").output().expect("binary not found");
    let combined = format!(
        "{}{}",
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(
        combined.contains("bismark2bedGraph") || combined.contains("Bismark"),
        "expected help/version string, got: {combined}"
    );
}

// ─── basic bedGraph output ────────────────────────────────────────────────────

#[test]
fn test_basic_bedgraph_output() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("CpG_OT_sample.txt");
    write_extractor_output(&input, &[
        "read1\t+\tchr1\t100\tZ",
        "read2\t-\tchr1\t100\tz",
        "read3\t+\tchr1\t200\tZ",
    ]);

    // --output must be a bare filename; use --dir for the directory.
    let output = Command::new(bin())
        .arg("--output").arg("sample.bedGraph")
        .arg("--dir").arg(dir.path())
        .arg("--no_header")
        .arg(&input)
        .output()
        .expect("binary failed");

    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));
    assert!(dir.path().join("sample.bedGraph.gz").exists() || dir.path().join("sample.bedGraph").exists(),
        "bedGraph output not created in {:?}", dir.path());
}

// ─── compat flags accepted without error ─────────────────────────────────────

#[test]
fn test_counts_flag_accepted() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("CpG_OT_sample.txt");
    write_extractor_output(&input, &["read1\t+\tchr1\t100\tZ"]);

    let output = Command::new(bin())
        .arg("--output").arg("sample.bedGraph")
        .arg("--dir").arg(dir.path())
        .arg("--no_header")
        .arg("--counts")
        .arg(&input)
        .output()
        .expect("binary failed");

    assert!(
        output.status.success(),
        "--counts caused unexpected failure: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}
