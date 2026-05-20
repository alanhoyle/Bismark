use anyhow::{bail, Context, Result};
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};

/// Read SAM/BAM via `samtools view -h`.
pub struct BamReader {
    child: Child,
    reader: BufReader<ChildStdout>,
}

impl BamReader {
    /// Open a BAM/SAM file for reading. Extra args are passed before the path.
    pub fn open(samtools: &str, path: &Path, extra_args: &[&str]) -> Result<Self> {
        let mut cmd = Command::new(samtools);
        cmd.arg("view").arg("-h");
        for a in extra_args {
            cmd.arg(a);
        }
        cmd.arg(path);
        cmd.stdout(Stdio::piped()).stderr(Stdio::inherit());

        let mut child = cmd
            .spawn()
            .with_context(|| format!("failed to spawn samtools for {}", path.display()))?;
        let stdout = child.stdout.take().unwrap();
        Ok(BamReader {
            child,
            reader: BufReader::with_capacity(1 << 20, stdout),
        })
    }

    pub fn lines(&mut self) -> &mut BufReader<ChildStdout> {
        &mut self.reader
    }

    /// Wait for the subprocess to finish and check exit status.
    pub fn finish(mut self) -> Result<()> {
        let status = self.child.wait()?;
        if !status.success() {
            bail!("samtools view exited with status {}", status);
        }
        Ok(())
    }
}

/// Write SAM records to a BAM file via `samtools view -bS`.
pub struct BamWriter {
    child: Child,
    writer: BufWriter<ChildStdin>,
}

impl BamWriter {
    pub fn open(samtools: &str, path: &Path) -> Result<Self> {
        Self::open_with_threads(samtools, path, 1)
    }

    pub fn open_with_threads(samtools: &str, path: &Path, threads: u32) -> Result<Self> {
        let mut child = Command::new(samtools)
            .args(["view", "-bS", "--threads"])
            .arg(threads.to_string())
            .arg("-")
            .arg("-o")
            .arg(path)
            .stdin(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .with_context(|| format!("failed to spawn samtools for {}", path.display()))?;
        let stdin = child.stdin.take().unwrap();
        Ok(BamWriter {
            child,
            writer: BufWriter::with_capacity(1 << 20, stdin),
        })
    }

    pub fn write_line(&mut self, line: &[u8]) -> Result<()> {
        self.writer.write_all(line)?;
        self.writer.write_all(b"\n")?;
        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        drop(self.writer); // close stdin so samtools knows we're done
        let status = self.child.wait()?;
        if !status.success() {
            bail!("samtools view -bS exited with status {}", status);
        }
        Ok(())
    }
}

/// Detect samtools on PATH or at a user-specified path.
pub fn find_samtools(user_path: Option<&str>) -> Result<String> {
    if let Some(p) = user_path {
        // Accept either a directory or a full path to the binary
        let path = Path::new(p);
        let binary = if path.is_dir() {
            path.join("samtools")
        } else {
            path.to_path_buf()
        };
        if binary.exists() {
            return Ok(binary.to_string_lossy().into_owned());
        }
        bail!(
            "Could not find samtools at '{}'. Please respecify with --samtools_path.",
            binary.display()
        );
    }
    // Check PATH
    if which("samtools") {
        return Ok("samtools".to_string());
    }
    bail!("No samtools installation found. Add samtools to PATH or use --samtools_path.")
}

fn which(name: &str) -> bool {
    Command::new("which")
        .arg(name)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

/// Check whether a BAM file appears truncated (samtools quickcheck).
pub fn bam_is_truncated(samtools: &str, path: &Path) -> bool {
    !Command::new(samtools)
        .args(["quickcheck", "-q"])
        .arg(path)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

/// Determine SE or PE mode by reading the @PG header of a BAM file.
/// Returns `(single_end, paired_end)`.
pub fn determine_file_type(samtools: &str, path: &Path) -> Result<(bool, bool)> {
    let mut cmd = Command::new(samtools);
    cmd.args(["view", "-H"]).arg(path);
    cmd.stdout(Stdio::piped()).stderr(Stdio::null());

    let output = cmd.output().context("samtools view -H failed")?;
    for line in output.stdout.split(|&b| b == b'\n') {
        if line.starts_with(b"@PG") {
            let s = String::from_utf8_lossy(line);
            if s.contains("--paired-end") || s.contains("-1 ") {
                return Ok((false, true));
            }
            // single-end is the default if no PE flags found in @PG
        }
    }
    Ok((true, false))
}

/// Assert a BAM file is sorted by name (not coordinate), as required by the
/// methylation extractor for PE mode.
pub fn assert_name_sorted(samtools: &str, path: &Path) -> Result<()> {
    let output = Command::new(samtools)
        .args(["view", "-H"])
        .arg(path)
        .output()
        .context("samtools view -H failed")?;

    for line in output.stdout.split(|&b| b == b'\n') {
        if line.starts_with(b"@HD") {
            let s = String::from_utf8_lossy(line);
            if s.contains("SO:coordinate") {
                bail!(
                    "BAM file {} appears to be coordinate-sorted. \
                     Paired-end processing requires name-sorted input. \
                     Re-sort with: samtools sort -n",
                    path.display()
                );
            }
        }
    }
    Ok(())
}
