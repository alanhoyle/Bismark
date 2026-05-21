# Bismark — Rust port

This directory contains a Rust reimplementation of the Bismark downstream
processing tools. The Rust binaries are drop-in replacements for their Perl
counterparts: they accept the same command-line flags and produce byte-identical
output.

The main Bismark aligner is not yet ported; alignment still uses the Perl
`bismark` script backed by Bowtie2, HISAT2, or minimap2.

## Tools

| Binary                          | Replaces                        | Purpose                                                                |
| ------------------------------- | ------------------------------- | ---------------------------------------------------------------------- |
| `bismark_methylation_extractor` | `bismark_methylation_extractor` | Extract per-cytosine methylation calls from Bismark BAM/SAM/CRAM files |
| `deduplicate_bismark`           | `deduplicate_bismark`           | Remove PCR duplicates from Bismark alignments                          |
| `bismark2bedGraph`              | `bismark2bedGraph`              | Convert methylation call files to bedGraph and coverage format         |
| `coverage2cytosine`             | `coverage2cytosine`             | Generate genome-wide cytosine methylation reports                      |
| `filter_non_conversion`         | `filter_non_conversion`         | Filter reads with high non-bisulfite-conversion rates                  |
| `bismark_genome_preparation`    | `bismark_genome_preparation`    | Prepare bisulfite-converted genome indices                             |
| `bam2nuc`                       | `bam2nuc`                       | Calculate nucleotide coverage from BAM files                           |
| `bismark2report`                | `bismark2report`                | Generate per-sample HTML reports                                       |
| `bismark2summary`               | `bismark2summary`               | Generate multi-sample summary report                                   |
| `methylation_consistency`       | `methylation_consistency`       | Assess read-level methylation consistency                              |
| `NOMe_filtering`                | `NOMe_filtering`                | Filter NOMe-seq cytosine reports                                       |

## Building

```bash
cd rust
cargo build --release --workspace
```

Binaries are written to `rust/target/release/`. Samtools must be available in
`PATH` (or passed via `--samtools_path`) at runtime.

**Minimum Rust version:** current stable (2021 edition)

## Using the Rust tools alongside the Perl aligner

`rust/bismark` is a wrapper script that prepends `rust/target/release/` to
`PATH` and then invokes the Perl `bismark` aligner. Running it means every
downstream tool spawned during or after alignment — `bismark_methylation_extractor`,
`deduplicate_bismark`, `bismark2bedGraph`, etc. — automatically resolves to the
Rust binary.

Use it in place of the bare `bismark` command:

```bash
./rust/bismark --genome /path/to/genome -1 r1.fq.gz -2 r2.fq.gz
```

Or add `rust/` to `PATH` for a session-wide effect:

```bash
export PATH="/path/to/Bismark/rust:$PATH"
bismark --genome /path/to/genome -1 r1.fq.gz -2 r2.fq.gz
```

## Testing

### Unit tests

```bash
cargo test --workspace
```

### Differential tests

The differential test suite runs both the Perl and Rust versions on the same
input and diffs every output file byte-for-byte:

```bash
./tests/differential_test.sh
```

Pass `--keep` to retain temporary output directories on failure, or
`--test-files` to run a full alignment first and then compare downstream tools.

### Performance benchmark

```bash
./tests/performance.sh
```

## Project layout

```
rust/
├── Cargo.toml                    # workspace manifest
├── bismark-lib/                  # shared library (BAM I/O, FASTA loading, …)
├── bismark-extractor/            # bismark_methylation_extractor
├── bismark-dedup/                # deduplicate_bismark
├── bismark-bedgraph/             # bismark2bedGraph
├── bismark-coverage2cytosine/    # coverage2cytosine
├── bismark-filter-non-conversion/
├── bismark-genome-prep/
├── bismark-bam2nuc/
├── bismark-report/
├── bismark-summary/
├── bismark-methylation-consistency/
├── bismark-nome-filtering/
└── tests/                        # integration / differential test scripts
```

## CLI compatibility

The Rust tools accept the same flags as the Perl originals. The notes below
cover flags where behaviour differs from the Perl version.

### `deduplicate_bismark`

- `--bam` — accepted for compatibility; BAM is already the default output
  format (Perl required this flag explicitly to get BAM output).
- `--sam` — writes plain SAM output; not available in the Perl version.
- `--representative` — was removed upstream; the Rust tool exits immediately
  with an error, matching Perl behaviour since that version.

### `coverage2cytosine`

- `--parent_dir` — accepted for compatibility; prints a deprecation notice and
  has no effect (the Rust version does not use `chdir` for path resolution).

### `bismark2bedGraph`

- `--buffer_size <size>` — enables a native Rust external merge sort rather
  than shelling out to the system `sort` utility. Sorted runs are flushed to
  temp files and merged with a k-way min-heap; no external process is spawned.
  Parse rules: `G`/`g` = GiB, `M`/`m` = MiB, `K`/`k` = KiB; bare numbers
  are bytes.
- `--ample_memory` — in-memory sort is the default when `--buffer_size` is
  not set, so this flag is a no-op.
- `--gazillion`/`--scaffolds` — no-op; the native sort handles any number of
  scaffolds without special casing.

Run any binary with `--help` for the full flag reference.
