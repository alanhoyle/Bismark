# Bismark — Rust port

This directory contains a Rust reimplementation of the Bismark downstream
processing tools, based on **[Bismark v0.25.1](https://github.com/FelixKrueger/Bismark)**
by Felix Krueger (Altos Bioinformatics). The Rust binaries are drop-in
replacements for their Perl counterparts: they accept the same command-line
flags and produce byte-identical output.

This port follows the [rewrites.bio](https://rewrites.bio/) principles:
credit the original authors, emulate outputs exactly, and be transparent
about AI assistance and validation.

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

## Attribution

Bismark was created by Felix Krueger and Simon Andrews at the Babraham Institute
and is now maintained at Altos Bioinformatics.

- **Upstream repository:** <https://github.com/FelixKrueger/Bismark>
- **Citation:** Krueger F & Andrews SR (2011). Bismark: A flexible aligner and
  methylation caller for Bismark-Seq applications. _Bioinformatics_ 27(11):1571–2.
  <https://doi.org/10.1093/bioinformatics/btr167>

This Rust port is a derivative work. Please cite the original Bismark paper when
using these tools in published research.

## AI assistance & validation

This port was written with [Claude Code](https://claude.ai/code) (Anthropic) and
[OpenAI Codex](https://openai.com/codex), with Codex used in test development.

**Correctness validation:**

- A differential test suite (`tests/differential_test.sh`) runs both the Perl
  and Rust implementations on identical inputs and diffs every output file
  byte-for-byte.
- End-to-end validation was performed against a fork of the
  [nf-core/methylseq](https://github.com/nf-core/methylseq) pipeline, comparing
  Perl and Rust outputs at every stage.

**Known gaps:**

- The main `bismark` aligner is not ported; only downstream tools are covered.
- Validation used paired-end Bowtie2 alignments. HISAT2, minimap2, and single-end
  modes have lighter test coverage.
- NOMe-seq and SLAM-seq code paths have not been validated against real data.
- `bismark2summary` sorts samples lexicographically; the Perl version emits them
  in filesystem glob order (non-deterministic). Row order in the summary report
  may differ from Perl when processing multiple samples.

## Building

```bash
cd rust
cargo build --release --workspace
```

Binaries are written to `rust/target/release/`. Samtools must be available in
`PATH` (or passed via `--samtools_path`) at runtime.

**Minimum Rust version:** current stable (2021 edition)

## Docker

A `Dockerfile` at the repository root builds a two-stage image: the Rust
downstream tools are compiled in a `rust:slim-bookworm` builder stage; the
runtime stage is `ubuntu:24.04` with Perl, samtools, Bowtie2, HISAT2, and
minimap2 pre-installed alongside the Rust binaries.

**Single-platform local build** (run from the repository root):

```bash
docker build -t bismark .
```

**Multi-platform build and push** (requires `docker buildx`):

```bash
docker buildx build --platform linux/amd64,linux/arm64 \
    -t bismark --push .
```

**Run** (mount your working directory as `/data`):

```bash
docker run --rm -v "$PWD:/data" -w /data bismark \
    bismark --genome /data/genome -1 r1.fq.gz -2 r2.fq.gz

docker run --rm -v "$PWD:/data" -w /data bismark \
    bismark_methylation_extractor --paired --comprehensive sample.bam
```

All Bismark tools — both the Perl aligner and the Rust downstream binaries —
are on `PATH` inside the container at `/bismark`.

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

Pass `--keep` to retain temporary output directories on failure, `--test-files`
to run a full alignment first and compare downstream tools, or supply your own
data directly:

```bash
./tests/differential_test.sh --fasta genome.fa.gz \
    --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz
```

### Performance benchmark

```bash
./tests/performance.sh
```

Supply custom data the same way:

```bash
./tests/performance.sh --fasta genome.fa.gz \
    --fastq1 R1.fastq.gz --fastq2 R2.fastq.gz
```

## Project layout

```text
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
