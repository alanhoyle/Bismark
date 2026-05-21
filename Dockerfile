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
# TARGETPLATFORM is injected by buildx; the Rust toolchain cross-compiles
# to the correct architecture automatically.
FROM --platform=$BUILDPLATFORM rust:1.75-slim-bookworm AS builder

ARG TARGETPLATFORM
ARG TARGETARCH
# Map Docker arch names to Rust target triples
RUN case "$TARGETARCH" in \
        amd64) echo x86_64-unknown-linux-gnu   > /target ;; \
        arm64) echo aarch64-unknown-linux-gnu  > /target ;; \
        *) echo "Unsupported arch: $TARGETARCH" >&2; exit 1 ;; \
    esac \
 && rustup target add "$(cat /target)"

WORKDIR /build
COPY rust/ .
RUN cargo build --release --workspace --target "$(cat /target)"

RUN mkdir /out \
 && cp target/"$(cat /target)"/release/bismark_methylation_extractor \
       target/"$(cat /target)"/release/deduplicate_bismark \
       target/"$(cat /target)"/release/bismark2bedGraph \
       target/"$(cat /target)"/release/coverage2cytosine \
       target/"$(cat /target)"/release/bismark_genome_preparation \
       target/"$(cat /target)"/release/filter_non_conversion \
       target/"$(cat /target)"/release/bam2nuc \
       target/"$(cat /target)"/release/bismark2report \
       target/"$(cat /target)"/release/bismark2summary \
       target/"$(cat /target)"/release/methylation_consistency \
       target/"$(cat /target)"/release/NOMe_filtering \
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
