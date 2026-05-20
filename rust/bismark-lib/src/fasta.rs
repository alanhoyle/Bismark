use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use indexmap::IndexMap;
use walkdir::WalkDir;

/// Reference genome: chromosome name → uppercase sequence bytes.
/// Uses `IndexMap` to preserve the FASTA file order, which matters for
/// cytosine traversal and SAM header construction.
pub struct Genome {
    pub sequences: IndexMap<String, Vec<u8>>,
}

impl Genome {
    /// Load all `*.fa`, `*.fa.gz`, `*.fasta`, `*.fasta.gz` files from
    /// `genome_folder`, replicating Bismark's `read_genome_into_memory()`.
    pub fn load(genome_folder: &Path) -> Result<Self> {
        let mut sequences: IndexMap<String, Vec<u8>> = IndexMap::new();

        let mut fasta_files: Vec<PathBuf> = WalkDir::new(genome_folder)
            .max_depth(1)
            .into_iter()
            .filter_map(|e| e.ok())
            .filter(|e| e.file_type().is_file())
            .map(|e| e.into_path())
            .filter(|p| is_fasta(p))
            .collect();

        if fasta_files.is_empty() {
            bail!(
                "No FASTA files found in genome folder: {}",
                genome_folder.display()
            );
        }

        fasta_files.sort();

        for path in &fasta_files {
            eprintln!("Loading genome from file: {}", path.display());
            load_fasta_file(path, &mut sequences)
                .with_context(|| format!("failed to load {}", path.display()))?;
        }

        eprintln!(
            "Genome loaded ({} sequences, {} total bp)",
            sequences.len(),
            sequences.values().map(|v| v.len()).sum::<usize>()
        );

        Ok(Genome { sequences })
    }

    /// Return a slice of the sequence at `chr`, positions `start..start+len`
    /// (0-based, half-open). Returns `None` if out of bounds.
    pub fn slice(&self, chr: &str, start: usize, len: usize) -> Option<&[u8]> {
        let seq = self.sequences.get(chr)?;
        let end = start + len;
        if end > seq.len() {
            None
        } else {
            Some(&seq[start..end])
        }
    }

    /// Get a single base (0-based). Returns `b'N'` if out of bounds.
    pub fn base(&self, chr: &str, pos: usize) -> u8 {
        self.sequences
            .get(chr)
            .and_then(|s| s.get(pos))
            .copied()
            .unwrap_or(b'N')
    }
}

fn is_fasta(p: &Path) -> bool {
    let name = p.file_name().unwrap_or_default().to_string_lossy();
    name.ends_with(".fa")
        || name.ends_with(".fa.gz")
        || name.ends_with(".fasta")
        || name.ends_with(".fasta.gz")
}

fn load_fasta_file(
    path: &Path,
    sequences: &mut IndexMap<String, Vec<u8>>,
) -> Result<()> {
    let file = std::fs::File::open(path)?;

    if path.extension().and_then(|e| e.to_str()) == Some("gz")
        || path.to_str().map(|s| s.ends_with(".fa.gz") || s.ends_with(".fasta.gz")).unwrap_or(false)
    {
        let decoder = MultiGzDecoder::new(file);
        parse_fasta(BufReader::new(decoder), sequences)
    } else {
        parse_fasta(BufReader::new(file), sequences)
    }
}

fn parse_fasta<R: BufRead>(
    reader: R,
    sequences: &mut IndexMap<String, Vec<u8>>,
) -> Result<()> {
    let mut current_name: Option<String> = None;
    let mut current_seq: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            // Flush previous sequence
            if let Some(name) = current_name.take() {
                if sequences.contains_key(&name) {
                    eprintln!(
                        "Warning: duplicate chromosome name '{}' — skipping duplicate",
                        name
                    );
                } else {
                    sequences.insert(name, current_seq.clone());
                }
                current_seq.clear();
            }
            // Extract first whitespace-delimited token as chromosome name
            let name = line[1..]
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            current_name = Some(name);
        } else if current_name.is_some() {
            // Append uppercase sequence bytes
            for b in line.bytes() {
                current_seq.push(b.to_ascii_uppercase());
            }
        }
    }

    // Flush last sequence
    if let Some(name) = current_name {
        if !sequences.contains_key(&name) {
            sequences.insert(name, current_seq);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_fasta_basic() {
        let fa = b">chr1 Some description\nACGTACGT\nACGT\n>chr2\nNNNN\n";
        let mut seqs = IndexMap::new();
        parse_fasta(Cursor::new(fa.as_ref()), &mut seqs).unwrap();
        assert_eq!(seqs["chr1"], b"ACGTACGTACGT");
        assert_eq!(seqs["chr2"], b"NNNN");
        // Insertion order preserved
        assert_eq!(seqs.get_index(0).unwrap().0, "chr1");
        assert_eq!(seqs.get_index(1).unwrap().0, "chr2");
    }

    #[test]
    fn test_genome_slice() {
        let fa = b">chr1\nACGTACGTACGT\n";
        let mut seqs = IndexMap::new();
        parse_fasta(std::io::Cursor::new(fa.as_ref()), &mut seqs).unwrap();
        let g = Genome { sequences: seqs };
        assert_eq!(g.slice("chr1", 0, 4), Some(b"ACGT".as_ref()));
        assert_eq!(g.slice("chr1", 8, 4), Some(b"ACGT".as_ref()));
        assert_eq!(g.slice("chr1", 10, 5), None); // out of bounds
        assert_eq!(g.base("chr1", 0), b'A');
        assert_eq!(g.base("chr2", 0), b'N'); // missing chr
    }
}
