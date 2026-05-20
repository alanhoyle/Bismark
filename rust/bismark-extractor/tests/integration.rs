// Integration tests for bismark_methylation_extractor.
// Each test writes a minimal SAM file, runs the binary, and checks the output.

use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_bismark_methylation_extractor"))
}

/// Write a minimal SAM with a header and one or more alignment lines.
fn write_sam(path: &PathBuf, records: &[&str]) {
    let mut f = fs::File::create(path).unwrap();
    writeln!(f, "@HD\tVN:1.6\tSO:unsorted").unwrap();
    writeln!(f, "@SQ\tSN:chr1\tLN:10000").unwrap();
    writeln!(f, "@PG\tID:Bismark\tPN:Bismark\tVN:v0.25.1\tCL:bismark --genome /g -1 r1.fq").unwrap();
    for rec in records {
        writeln!(f, "{rec}").unwrap();
    }
}

/// Single-end OT read: XR=CT, XG=CT, FLAG=0 (forward), POS=100, 10M.
/// XM="ZzXxHh...." → Z at pos 100, z at 101, X at 102, x at 103, H at 104, h at 105.
fn ot_se_record() -> &'static str {
    "read1\t0\tchr1\t100\t255\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tXM:Z:ZzXxHh....\tXR:Z:CT\tXG:Z:CT"
}

/// Single-end OB read: XR=CT, XG=GA, FLAG=16 (reverse), POS=200, 10M.
/// XM="ZzXxHh...." is in SAM/+ orientation: XM[0]='Z' maps to leftmost pos = POS = 200.
/// Reverse formula: start = end_pos = 209, pos = start - reversed_i → Z at 200, z at 201.
fn ob_se_record() -> &'static str {
    "read2\t16\tchr1\t200\t255\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tXM:Z:ZzXxHh....\tXR:Z:CT\tXG:Z:GA"
}

fn pe_overlap_records() -> [&'static str; 2] {
    [
        "pair1/1\t99\tchr1\t100\t255\t10M\t=\t105\t15\tACGTACGTAC\tIIIIIIIIII\tXM:Z:ZzZzZzZzZz\tXR:Z:CT\tXG:Z:CT",
        "pair1/2\t147\tchr1\t105\t255\t10M\t=\t100\t-15\tACGTACGTAC\tIIIIIIIIII\tXM:Z:ZzZzZzZzZz\tXR:Z:CT\tXG:Z:GA",
    ]
}

fn read_context_file(dir: &std::path::Path, needle: &str) -> String {
    let path = fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .find(|p| p.file_name().unwrap().to_string_lossy().contains(needle))
        .unwrap_or_else(|| panic!("{needle} output not found in {dir:?}"));
    fs::read_to_string(path).unwrap()
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
        combined.contains("0.25.1") || combined.contains("Bismark"),
        "expected version string, got: {combined}"
    );
}

// ─── Single-end OT extraction (strand-specific mode) ─────────────────────────

#[test]
fn test_se_ot_extraction() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("se_ot.sam");
    write_sam(&sam, &[ot_se_record()]);

    // Strand-specific mode (default): creates CpG_OT_<stem>.txt, CHG_OT_<stem>.txt, etc.
    let output = Command::new(bin())
        .args(["--single", "--no_header"])
        .arg("--dir").arg(dir.path())
        .arg(&sam)
        .output()
        .expect("binary failed");

    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));

    let cpg_ot: Vec<_> = fs::read_dir(dir.path())
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| {
            let name = e.file_name();
            let name = name.to_string_lossy();
            name.contains("CpG_OT") && name.ends_with(".txt")
        })
        .collect();

    assert!(!cpg_ot.is_empty(), "CpG_OT output not found in {:?}", dir.path());

    let content = fs::read_to_string(cpg_ot[0].path()).unwrap();
    assert!(content.contains("chr1"), "chr1 not in CpG_OT output:\n{content}");
    // Z at position 100 (methylated CpG on OT strand)
    assert!(content.contains("100"), "position 100 not in CpG_OT output:\n{content}");
    // z at position 101 (unmethylated CpG on OT strand)
    assert!(content.contains("101"), "position 101 not in CpG_OT output:\n{content}");
}

// ─── Single-end OB extraction (reverse read, fixed position logic) ────────────

#[test]
fn test_se_ob_extraction() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("se_ob.sam");
    write_sam(&sam, &[ob_se_record()]);

    let output = Command::new(bin())
        .args(["--single", "--no_header"])
        .arg("--dir").arg(dir.path())
        .arg(&sam)
        .output()
        .expect("binary failed");

    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));

    let cpg_ob: Vec<_> = fs::read_dir(dir.path())
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| {
            let name = e.file_name();
            let name = name.to_string_lossy();
            name.contains("CpG_OB") && name.ends_with(".txt")
        })
        .collect();

    assert!(!cpg_ob.is_empty(), "CpG_OB output not found");

    let content = fs::read_to_string(cpg_ob[0].path()).unwrap();
    assert!(content.contains("chr1"), "chr1 not in OB output:\n{content}");

    // OB: XM in SAM/+ orientation; XM[0]='Z' at POS=200, XM[1]='z' at 201.
    // Reverse formula: end_pos=209, pos = end_pos - reversed_i gives Z@200, z@201.
    assert!(content.contains("200"), "Z at pos 200 not in CpG_OB output:\n{content}");
    assert!(content.contains("201"), "z at pos 201 not in CpG_OB output:\n{content}");
}

// ─── Splitting report is written with --report flag ──────────────────────────

#[test]
fn test_splitting_report_written() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("se.sam");
    write_sam(&sam, &[ot_se_record()]);

    let output = Command::new(bin())
        .args(["--single", "--no_header", "--report"])
        .arg("--dir").arg(dir.path())
        .arg(&sam)
        .output()
        .expect("binary failed");

    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));

    let report: Vec<_> = fs::read_dir(dir.path())
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name().to_string_lossy().contains("splitting_report"))
        .collect();

    assert!(!report.is_empty(), "splitting report not found in {:?}", dir.path());
    let content = fs::read_to_string(report[0].path()).unwrap();
    assert!(content.contains("Bismark"), "report header missing:\n{content}");
    assert!(content.contains("CpG"), "CpG line missing:\n{content}");
}

// ─── Comprehensive mode creates context files ─────────────────────────────────

#[test]
fn test_comprehensive_mode_creates_context_files() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("comp.sam");
    write_sam(&sam, &[ot_se_record()]);

    let output = Command::new(bin())
        .args(["--single", "--no_header", "--comprehensive"])
        .arg("--dir").arg(dir.path())
        .arg(&sam)
        .output()
        .expect("binary failed");

    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));

    let ctx_files: Vec<_> = fs::read_dir(dir.path())
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name().to_string_lossy().contains("_context_"))
        .collect();

    assert_eq!(ctx_files.len(), 3, "expected CpG/CHG/CHH context files, got {:?}", ctx_files);
}

// ─── Paired-end overlap handling ─────────────────────────────────────────────

#[test]
fn test_pe_no_overlap_is_default() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("pe.sam");
    let records = pe_overlap_records();
    write_sam(&sam, &records);

    let output = Command::new(bin())
        .args(["--paired", "--no_header", "--comprehensive"])
        .arg("--dir").arg(dir.path())
        .arg(&sam)
        .output()
        .expect("binary failed");

    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));

    let cpg = read_context_file(dir.path(), "CpG_context");
    assert_eq!(cpg.lines().count(), 15, "default PE mode should suppress overlapping R2 calls:\n{cpg}");
    assert!(cpg.contains("pair1/2\t-\tchr1\t110\tz"), "expected first non-overlapping R2 call:\n{cpg}");
    assert!(!cpg.contains("pair1/2\t+\tchr1\t109\tZ"), "overlapping R2 call should be suppressed:\n{cpg}");
}

#[test]
fn test_pe_include_overlap_keeps_all_calls() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("pe.sam");
    let records = pe_overlap_records();
    write_sam(&sam, &records);

    let output = Command::new(bin())
        .args(["--paired", "--include_overlap", "--no_header", "--comprehensive"])
        .arg("--dir").arg(dir.path())
        .arg(&sam)
        .output()
        .expect("binary failed");

    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));

    let cpg = read_context_file(dir.path(), "CpG_context");
    assert_eq!(cpg.lines().count(), 20, "--include_overlap should keep all R1 and R2 calls:\n{cpg}");
    assert!(cpg.contains("pair1/2\t+\tchr1\t109\tZ"), "overlapping R2 call should be present:\n{cpg}");
}
