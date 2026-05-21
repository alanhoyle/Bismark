// Integration tests for NOMe_filtering.
//
// The tool requires a genome folder (FASTA) and a YACHT-format input file.
// Tests use a minimal synthetic genome and input to verify flag acceptance.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_NOMe_filtering"))
}

/// Write a minimal genome FASTA to `dir/chr1.fa`.
/// Sequence "NNNNATACGTNNN" has a CpG at positions 8–9 (1-based) in ACG context.
fn write_genome(dir: &Path) {
    let fa = dir.join("chr1.fa");
    let mut f = fs::File::create(&fa).unwrap();
    // ACG at pos 4 (0-based), so genomic pos 5 (1-based) has tri_nt CGT → CpG,
    // upstream = ACG → valid NOMe CpG.
    writeln!(f, ">chr1").unwrap();
    writeln!(f, "NNNNATACGTNNNN").unwrap();
}

/// Write a minimal YACHT-format input file.
/// ReadID  state  chr  pos  context  start  end  strand
fn write_yacht(path: &Path, records: &[&str]) {
    let mut f = fs::File::create(path).unwrap();
    for rec in records {
        writeln!(f, "{rec}").unwrap();
    }
}

/// One read with a methylated CpG at position 8 (1-based, 'A' at pos 7, 'C' at 8, 'G' at 9).
/// The genome "NNNNATACGTNNNN" has: pos 8 = 'C', pos 9 = 'G' → CpG.
/// Upstream (pos 7,8,9) = 'A','C','G' = ACG → valid NOMe CpG call.
fn yacht_record() -> &'static str {
    "read1\t+\tchr1\t8\tZ\t1\t14\t+"
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
        combined.contains("0.25.1") || combined.contains("NOMe") || combined.contains("Bismark"),
        "expected version string, got: {combined}"
    );
}

// ─── basic run ────────────────────────────────────────────────────────────────

#[test]
fn test_basic_run() {
    let dir = tempfile::tempdir().unwrap();
    let genome_dir = dir.path().join("genome");
    fs::create_dir(&genome_dir).unwrap();
    write_genome(&genome_dir);

    // The tool strips the input extension and adds .manOwar.txt.gz — use a
    // plain name so the output path is predictable.
    let yacht = dir.path().join("sample.txt");
    write_yacht(&yacht, &[yacht_record()]);

    let output = Command::new(bin())
        .arg("--genome_folder").arg(&genome_dir)
        .arg("--dir").arg(dir.path())
        .arg(&yacht)
        .output()
        .expect("binary failed");

    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));

    let out_gz = dir.path().join("sample.manOwar.txt.gz");
    assert!(out_gz.exists(), "output .gz not created at {out_gz:?}");
}

// ─── Perl-compat flags accepted without error ─────────────────────────────────
//
// All flags below are accepted by the Perl NOMe_filtering but are effectively
// no-ops in that implementation (parsed but never wired into the processing).
// The Rust port accepts them for drop-in compatibility.

fn run_with_flag(flag: &str) -> bool {
    let dir = tempfile::tempdir().unwrap();
    let genome_dir = dir.path().join("genome");
    fs::create_dir(&genome_dir).unwrap();
    write_genome(&genome_dir);

    let yacht = dir.path().join("sample.txt");
    write_yacht(&yacht, &[yacht_record()]);

    Command::new(bin())
        .arg("--genome_folder").arg(&genome_dir)
        .arg("--dir").arg(dir.path())
        .arg(flag)
        .arg(&yacht)
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

#[test]
fn test_zero_based_flag_accepted() {
    assert!(run_with_flag("--zero_based"), "--zero_based caused unexpected failure");
}

#[test]
fn test_cx_context_flag_accepted() {
    assert!(run_with_flag("--CX"), "--CX caused unexpected failure");
}

#[test]
fn test_merge_cpgs_flag_accepted() {
    assert!(run_with_flag("--merge_CpGs"), "--merge_CpGs caused unexpected failure");
}

#[test]
fn test_gc_context_flag_accepted() {
    assert!(run_with_flag("--GC"), "--GC caused unexpected failure");
}

#[test]
fn test_gzip_flag_accepted() {
    assert!(run_with_flag("--gzip"), "--gzip caused unexpected failure");
}

#[test]
fn test_nome_seq_flag_accepted() {
    assert!(run_with_flag("--nome-seq"), "--nome-seq caused unexpected failure");
}

#[test]
fn test_parent_dir_flag_accepted() {
    let dir = tempfile::tempdir().unwrap();
    let genome_dir = dir.path().join("genome");
    fs::create_dir(&genome_dir).unwrap();
    write_genome(&genome_dir);

    let yacht = dir.path().join("sample.txt");
    write_yacht(&yacht, &[yacht_record()]);

    let output = Command::new(bin())
        .arg("--genome_folder").arg(&genome_dir)
        .arg("--dir").arg(dir.path())
        .arg("--parent_dir").arg(dir.path())
        .arg(&yacht)
        .output()
        .expect("binary failed");

    assert!(output.status.success(), "--parent_dir caused failure: {}", String::from_utf8_lossy(&output.stderr));
}
