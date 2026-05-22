use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{bail, Result};
use flate2::read::MultiGzDecoder;

#[derive(Debug, Clone)]
pub struct FastqRecord {
    #[allow(dead_code)]
    /// Read name, without the leading `@`, trimmed at first whitespace.
    pub id: Vec<u8>,
    /// Raw read name including description (after `@`).
    pub id_full: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

impl FastqRecord {
    #[allow(dead_code)]
    /// Truncated name used for SAM QNAME (everything up to first space/tab).
    pub fn qname(&self) -> &[u8] {
        &self.id
    }
}

pub struct FastqReader<R: BufRead> {
    inner: R,
    line: Vec<u8>,
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        FastqReader {
            inner: reader,
            line: Vec::with_capacity(256),
        }
    }
}

impl<R: BufRead> Iterator for FastqReader<R> {
    type Item = Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line.clear();
            match self.inner.read_until(b'\n', &mut self.line) {
                Err(e) => return Some(Err(e.into())),
                Ok(0) => return None,
                Ok(_) => {}
            }
            let trimmed = trim_end(&self.line);
            if trimmed.is_empty() {
                continue;
            }
            if trimmed[0] != b'@' {
                continue;
            }
            let header = trimmed[1..].to_vec();
            let id = header
                .iter()
                .position(|&b| b == b' ' || b == b'\t')
                .map(|i| header[..i].to_vec())
                .unwrap_or_else(|| header.clone());

            let mut seq = Vec::new();
            self.line.clear();
            match self.inner.read_until(b'\n', &mut self.line) {
                Err(e) => return Some(Err(e.into())),
                Ok(0) => return Some(Err(anyhow::anyhow!("truncated FASTQ: missing seq line"))),
                Ok(_) => seq.extend_from_slice(trim_end(&self.line)),
            }

            // separator line (+)
            self.line.clear();
            match self.inner.read_until(b'\n', &mut self.line) {
                Err(e) => return Some(Err(e.into())),
                Ok(0) => return Some(Err(anyhow::anyhow!("truncated FASTQ: missing + line"))),
                Ok(_) => {}
            }

            let mut qual = Vec::new();
            self.line.clear();
            match self.inner.read_until(b'\n', &mut self.line) {
                Err(e) => return Some(Err(e.into())),
                Ok(0) => return Some(Err(anyhow::anyhow!("truncated FASTQ: missing qual line"))),
                Ok(_) => qual.extend_from_slice(trim_end(&self.line)),
            }

            return Some(Ok(FastqRecord {
                id,
                id_full: header,
                seq,
                qual,
            }));
        }
    }
}

fn trim_end(s: &[u8]) -> &[u8] {
    let mut end = s.len();
    while end > 0 && (s[end - 1] == b'\n' || s[end - 1] == b'\r') {
        end -= 1;
    }
    &s[..end]
}

/// Open a FASTQ file (plain or gzip-compressed) and return a boxed `BufRead`.
pub fn open_fastq(path: &Path) -> Result<Box<dyn BufRead>> {
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    if name.ends_with(".gz") {
        let file = std::fs::File::open(path)?;
        Ok(Box::new(BufReader::with_capacity(
            1 << 20,
            MultiGzDecoder::new(file),
        )))
    } else {
        let file = std::fs::File::open(path)?;
        Ok(Box::new(BufReader::with_capacity(1 << 20, file)))
    }
}

/// Read all records from a FASTQ file into a Vec.
pub fn read_all(path: &Path) -> Result<Vec<FastqRecord>> {
    let reader = open_fastq(path)?;
    let mut records = Vec::new();
    for rec in FastqReader::new(reader) {
        records.push(rec?);
    }
    Ok(records)
}

/// Read from a paired pair of paths, returning interleaved (r1, r2) pairs.
pub fn read_all_pe(r1_path: &Path, r2_path: &Path) -> Result<Vec<(FastqRecord, FastqRecord)>> {
    let r1 = open_fastq(r1_path)?;
    let r2 = open_fastq(r2_path)?;
    let mut iter1 = FastqReader::new(r1);
    let mut iter2 = FastqReader::new(r2);
    let mut pairs = Vec::new();
    loop {
        match (iter1.next(), iter2.next()) {
            (None, None) => break,
            (Some(a), Some(b)) => pairs.push((a?, b?)),
            _ => bail!("R1 and R2 FASTQ files have different numbers of records"),
        }
    }
    Ok(pairs)
}

/// Write FASTQ records to a writer (used by the aligner stdin writer thread).
pub fn write_record<W: std::io::Write>(w: &mut W, rec: &FastqRecord) -> std::io::Result<()> {
    w.write_all(b"@")?;
    w.write_all(&rec.id_full)?;
    w.write_all(b"\n")?;
    w.write_all(&rec.seq)?;
    w.write_all(b"\n+\n")?;
    w.write_all(&rec.qual)?;
    w.write_all(b"\n")
}

/// Read a FASTA file (plain or gzip) and return records as FastqRecord.
/// Quality scores are set to 'I' (Phred 40) for all bases, matching Perl bismark.
pub fn read_fasta(path: &Path) -> Result<Vec<FastqRecord>> {
    let reader = open_fastq(path)?;
    let mut records = Vec::new();
    let mut cur_id: Vec<u8> = Vec::new();
    let mut cur_id_full: Vec<u8> = Vec::new();
    let mut in_record = false;
    let mut cur_seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if in_record {
                let qual = vec![b'I'; cur_seq.len()];
                records.push(FastqRecord {
                    id: cur_id.clone(),
                    id_full: cur_id_full.clone(),
                    seq: std::mem::take(&mut cur_seq),
                    qual,
                });
            }
            let header = line[1..].as_bytes().to_vec();
            let id = header
                .iter()
                .position(|&b| b == b' ' || b == b'\t')
                .map(|i| header[..i].to_vec())
                .unwrap_or_else(|| header.clone());
            cur_id = id;
            cur_id_full = header;
            in_record = true;
        } else if in_record {
            cur_seq.extend_from_slice(line.as_bytes());
        }
    }
    if in_record {
        let qual = vec![b'I'; cur_seq.len()];
        records.push(FastqRecord {
            id: cur_id,
            id_full: cur_id_full,
            seq: cur_seq,
            qual,
        });
    }
    Ok(records)
}

/// Read paired FASTA files as FastqRecord pairs.
pub fn read_fasta_pe(r1_path: &Path, r2_path: &Path) -> Result<Vec<(FastqRecord, FastqRecord)>> {
    let r1 = read_fasta(r1_path)?;
    let r2 = read_fasta(r2_path)?;
    if r1.len() != r2.len() {
        bail!(
            "R1 and R2 FASTA files have different numbers of records ({} vs {})",
            r1.len(),
            r2.len()
        );
    }
    Ok(r1.into_iter().zip(r2).collect())
}
