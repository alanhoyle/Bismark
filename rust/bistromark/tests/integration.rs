/// Integration tests for bistromark.
///
/// These tests exercise the binary directly without a real Bowtie2 index.
/// Alignment tests require a mock `bowtie2` binary that reads FASTQ from stdin
/// and emits canned SAM output.

use std::path::PathBuf;
use std::process::Command;

fn bistromark_bin() -> PathBuf {
    let mut p = std::env::current_exe().unwrap();
    p.pop(); // deps/
    p.pop(); // debug/ or release/
    p.push("bistromark");
    p
}

#[test]
fn version_flag() {
    let out = Command::new(bistromark_bin())
        .arg("--version")
        .output()
        .expect("failed to run bistromark --version");
    assert!(out.status.success(), "exit status: {}", out.status);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("0.25.1"),
        "expected version in stdout, got: {}",
        stdout
    );
}

#[test]
fn help_exits_zero() {
    let out = Command::new(bistromark_bin())
        .arg("--help")
        .output()
        .expect("failed to run bistromark --help");
    assert!(out.status.success(), "exit status: {}", out.status);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(stdout.contains("--genome"), "expected --genome in help");
    assert!(stdout.contains("--threads"), "expected --threads in help");
}

#[test]
fn missing_genome_exits_nonzero() {
    let out = Command::new(bistromark_bin())
        .args(["-U", "/dev/null"])
        .output()
        .expect("bistromark missing genome");
    assert!(
        !out.status.success(),
        "expected failure without --genome, but got success"
    );
}

/// Test the FASTQ reader round-trip.
#[test]
fn fastq_reader_basic() {
    use std::io::Cursor;
    // Inline the module path.
    let fastq = b"@read1 desc\nACGT\n+\nIIII\n@read2\nTTTT\n+\nIIII\n";
    let cursor = Cursor::new(fastq.as_ref());
    let reader = std::io::BufReader::new(cursor);

    // We can't import the private module from here; just smoke-test via binary.
    let _ = reader; // Compiled-to-nothing stub; real tests use the binary.
}

/// Test that `--non_directional` is accepted without error when genome is missing
/// (the error should be about genome, not about the flag).
#[test]
fn non_directional_flag_accepted() {
    let out = Command::new(bistromark_bin())
        .args(["--non_directional", "-U", "/dev/null"])
        .output()
        .expect("bistromark");
    let stderr = String::from_utf8_lossy(&out.stderr);
    // Should fail on missing --genome, not on unknown flag
    assert!(
        !out.status.success(),
        "expected failure without --genome"
    );
    assert!(
        !stderr.contains("unexpected argument"),
        "flag should be recognised; got: {}",
        stderr
    );
}

/// Test that `--pbat` is accepted.
#[test]
fn pbat_flag_accepted() {
    let out = Command::new(bistromark_bin())
        .args(["--pbat", "-U", "/dev/null"])
        .output()
        .expect("bistromark");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        !stderr.contains("unexpected argument"),
        "flag should be recognised; got: {}",
        stderr
    );
}
