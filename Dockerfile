###############################################################################
# Bismark — Rust-enabled build
#
# The Perl bismark aligner handles alignment via Bowtie2 (included).
# All downstream tools (extraction, dedup, bedGraph, coverage report, …) run
# as native Rust binaries built in stage 1.
#
# Multi-platform build (linux/amd64 + linux/arm64):
#   docker buildx build --platform linux/amd64,linux/arm64 \
#       -t bismark --push .
#
# Single-platform local build:
#   docker build -t bismark .
#
# Run (mount your working directory as /data):
#   docker run --rm -v "$PWD:/data" -w /data bismark \
#       bismark --genome /data/genome -1 r1.fq.gz -2 r2.fq.gz
#
#   docker run --rm -v "$PWD:/data" -w /data bismark \
#       bismark_methylation_extractor --paired --comprehensive sample.bam
#
# To use HISAT2 or minimap2 instead of Bowtie2, build a derived image that
# installs those aligners and pass --hisat2 / --minimap2 to bismark.
###############################################################################

# ── Stage 1: compile Rust downstream tools ────────────────────────────────────
# No --platform override: each target arch builds its own binaries natively
# under QEMU emulation when cross-building, avoiding cross-linker complexity.
FROM rust:slim-bookworm AS builder

RUN apt-get update \
 && apt-get upgrade -y \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY rust/ .
# bismark-report and bismark-summary embed plotly assets via include_str! with
# paths that go three directories above their source files (../../../plotly/).
# With WORKDIR=/build that resolves to /plotly/, so copy the assets there.
COPY plotly/ /plotly/
RUN cargo build --release --workspace

RUN mkdir /out \
 && cp target/release/bismark_methylation_extractor \
       target/release/deduplicate_bismark \
       target/release/bismark2bedGraph \
       target/release/coverage2cytosine \
       target/release/bismark_genome_preparation \
       target/release/filter_non_conversion \
       target/release/bam2nuc \
       target/release/bismark2report \
       target/release/bismark2summary \
       target/release/methylation_consistency \
       target/release/NOMe_filtering \
       /out/

# ── Stage 2: runtime ──────────────────────────────────────────────────────────
FROM ubuntu:24.04

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        perl \
        samtools \
        bowtie2 \
 && rm -rf /var/lib/apt/lists/*

# All tools live in /bismark:
#   bismark          — Perl aligner (invokes bowtie2, writes BAM + tags)
#   everything else  — Rust binaries (drop-in replacements for Perl originals)
WORKDIR /bismark

COPY bismark .
COPY --from=builder /out/ ./

ENV PATH="/bismark:$PATH"

WORKDIR /data
CMD ["bismark", "--help"]
