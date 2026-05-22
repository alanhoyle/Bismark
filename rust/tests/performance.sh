#!/usr/bin/env bash
# Performance comparison: run Perl and Rust versions of Bismark downstream
# tools on synthetic inputs and report wall-clock timing summaries.
#
# Usage:
#   ./rust/tests/performance.sh [--runs N] [--records N] [--threads N]
#                                [--mem] [--keep] [--test-files]
#
# Requirements:
#   - samtools in PATH
#   - Perl scripts at ../  (relative to rust/)
#   - Rust binaries built: cargo build --release --workspace
#   - bowtie2/bowtie2-build in PATH when using --test-files
#
# Options:
#   --runs N      Number of benchmark repetitions [default: 3]
#   --records N   Synthetic input size [default: 50000]
#   --threads N   Thread count forwarded to tools that support it [default: 1]
#   --mem         Track peak RSS memory usage via /usr/bin/time.
#                 Prefers gtime (brew install gnu-time) on macOS for clean
#                 stderr separation; falls back to /usr/bin/time -l otherwise.
#   --keep        Keep temporary output directories on failure for inspection
#   --test-files  Use test_files/ data instead of synthetic inputs
#   --fasta FILE  Reference genome FASTA (replaces test_files/NC_010473.fa.gz);
#                 implies --test-files
#   --fastq1 FILE R1 FASTQ (replaces test_files/test_R1.fastq.gz)
#   --fastq2 FILE R2 FASTQ (replaces test_files/test_R2.fastq.gz)
#
# Notes:
#   This is a lightweight benchmark harness, not a statistical benchmark suite.
#   Run it on an otherwise quiet machine for the most useful numbers.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RUST_BIN="$SCRIPT_DIR/../target/release"
PERL_BIN="$REPO_ROOT"
TEST_FILES="$REPO_ROOT/test_files"

RUNS=3
RECORDS=50000
THREADS=1
KEEP=0
MEASURE_MEM=0
USE_TEST_FILES=0
CUSTOM_FASTA=""
CUSTOM_FASTQ1=""
CUSTOM_FASTQ2=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --runs)
            RUNS="${2:-}"
            shift 2
            ;;
        --records)
            RECORDS="${2:-}"
            shift 2
            ;;
        --threads)
            THREADS="${2:-}"
            shift 2
            ;;
        --mem)
            MEASURE_MEM=1
            shift
            ;;
        --keep)
            KEEP=1
            shift
            ;;
        --test-files)
            USE_TEST_FILES=1
            shift
            ;;
        --fasta)
            CUSTOM_FASTA="${2:-}"; USE_TEST_FILES=1
            shift 2
            ;;
        --fastq1)
            CUSTOM_FASTQ1="${2:-}"; USE_TEST_FILES=1
            shift 2
            ;;
        --fastq2)
            CUSTOM_FASTQ2="${2:-}"; USE_TEST_FILES=1
            shift 2
            ;;
        -h|--help)
            sed -n '1,31p' "$0"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 2
            ;;
    esac
done

die() { echo "FATAL: $*" >&2; exit 1; }

# ─── Memory measurement setup ─────────────────────────────────────────────────
# TIME_CMD / TIME_ARGS: the command used to wrap benchmarks when --mem is set.
# TIME_MODE: "gnu" (gtime or Linux time -v) or "darwin" (/usr/bin/time -l).
# TIME_SUPPORTS_OUTFILE: 1 when the -o flag is available to write time output
#   to a separate file (keeps command stderr clean); 0 on macOS fallback.
TIME_CMD=""
TIME_ARGS=""
TIME_MODE=""
TIME_SUPPORTS_OUTFILE=0

setup_time_cmd() {
    if [[ "$MEASURE_MEM" -eq 0 ]]; then
        return
    fi
    local platform
    platform="$(uname -s)"
    if command -v gtime >/dev/null 2>&1; then
        TIME_CMD="gtime"
        TIME_ARGS="-v"
        TIME_MODE="gnu"
        TIME_SUPPORTS_OUTFILE=1
    elif [[ "$platform" == "Linux" ]]; then
        TIME_CMD="/usr/bin/time"
        TIME_ARGS="-v"
        TIME_MODE="gnu"
        TIME_SUPPORTS_OUTFILE=1
    elif [[ "$platform" == "Darwin" ]]; then
        TIME_CMD="/usr/bin/time"
        TIME_ARGS="-l"
        TIME_MODE="darwin"
        TIME_SUPPORTS_OUTFILE=0
        echo "Note: gtime not found; falling back to /usr/bin/time -l." \
             "Command stderr and timing output will be mixed in .timemem files." \
             "(Install: brew install gnu-time)" >&2
    else
        die "--mem: cannot find a supported time command (try: brew install gnu-time)"
    fi
}

# Extract peak RSS in bytes from a time output file.
extract_rss() {
    local timelog="$1"
    [[ -f "$timelog" ]] || { echo 0; return; }
    if [[ "$TIME_MODE" == "darwin" ]]; then
        # /usr/bin/time -l line: "  32751616  maximum resident set size"
        grep "maximum resident set size" "$timelog" \
            | awk '{print $1}' | head -1 || echo 0
    else
        # gtime -v / time -v line: "Maximum resident set size (kbytes): 31984"
        local kb
        kb="$(grep "Maximum resident set size" "$timelog" \
              | awk '{print $NF}' | head -1)" || true
        echo "$(( ${kb:-0} * 1024 ))"
    fi
}

# ─── Timing helpers ───────────────────────────────────────────────────────────

now_seconds() {
    perl -MTime::HiRes=time -e 'printf "%.6f\n", time'
}

elapsed_seconds() {
    local start="$1" end="$2"
    awk -v s="$start" -v e="$end" 'BEGIN { printf "%.6f", e - s }'
}

# Globals set by time_command; read by append_result in the same shell.
# time_command must NOT be called inside $() — that creates a subshell and the
# globals would be invisible to the caller.
LAST_LOG_PREFIX=""
LAST_ELAPSED="0"

time_command() {
    LAST_LOG_PREFIX="$1"; shift
    local start end
    start="$(now_seconds)"
    if [[ "$MEASURE_MEM" -eq 1 ]]; then
        local timemem="${LAST_LOG_PREFIX}.timemem"
        if [[ "$TIME_SUPPORTS_OUTFILE" -eq 1 ]]; then
            # Write time stats to a separate file; command stderr goes to .err
            "$TIME_CMD" $TIME_ARGS -o "$timemem" \
                "$@" >"${LAST_LOG_PREFIX}.out" 2>"${LAST_LOG_PREFIX}.err"
        else
            # macOS fallback: command stderr and time stats both go to .timemem
            "$TIME_CMD" $TIME_ARGS \
                "$@" >"${LAST_LOG_PREFIX}.out" 2>"$timemem"
        fi
    else
        "$@" >"${LAST_LOG_PREFIX}.out" 2>"${LAST_LOG_PREFIX}.err"
    fi
    end="$(now_seconds)"
    LAST_ELAPSED="$(elapsed_seconds "$start" "$end")"
}

# ─── Statistics helpers ───────────────────────────────────────────────────────

ratio() {
    local perl_time="$1" rust_time="$2"
    awk -v p="$perl_time" -v r="$rust_time" 'BEGIN {
        if (r == 0) { printf "inf" }
        else { printf "%.2fx", p / r }
    }'
}

mem_ratio() {
    local perl_mem="$1" rust_mem="$2"
    awk -v p="$perl_mem" -v r="$rust_mem" 'BEGIN {
        if (r == 0) { printf "inf" }
        else { printf "%.2fx", p / r }
    }'
}

mean_csv() {
    local csv="$1"
    awk -F, '$3 ~ /^[0-9.]+$/ { sum += $3; n++ } \
             END { if (n) printf "%.6f", sum / n; else printf "0.000000" }' "$csv"
}

min_csv() {
    local csv="$1"
    awk -F, '$3 ~ /^[0-9.]+$/ { if (n == 0 || $3 < min) min = $3; n++ } \
             END { if (n) printf "%.6f", min; else printf "0.000000" }' "$csv"
}

# Mean peak RSS across runs (column 5), in bytes.
mean_rss_csv() {
    local csv="$1"
    awk -F, '$5 ~ /^[0-9]+$/ { sum += $5; n++ } \
             END { if (n) printf "%.0f", sum / n; else printf "0" }' "$csv"
}

# Max peak RSS across runs (column 5), in bytes.
max_rss_csv() {
    local csv="$1"
    awk -F, '$5 ~ /^[0-9]+$/ { if (n == 0 || $5 > max) max = $5; n++ } \
             END { if (n) printf "%.0f", max; else printf "0" }' "$csv"
}

format_bytes() {
    local bytes="$1"
    awk -v b="$bytes" 'BEGIN {
        if      (b >= 1073741824) printf "%.1f GiB", b / 1073741824
        else if (b >= 1048576)    printf "%.1f MiB", b / 1048576
        else if (b >= 1024)       printf "%.1f KiB", b / 1024
        else                      printf "%d B",     b
    }'
}

append_result() {
    local case_name="$1" impl="$2" run="$3" seconds="$4"
    local rss=0
    if [[ "$MEASURE_MEM" -eq 1 ]]; then
        rss="$(extract_rss "${LAST_LOG_PREFIX}.timemem")"
        rss="${rss:-0}"
    fi
    printf "%s,%s,%s,%s,%s\n" "$case_name" "$impl" "$seconds" "$run" "$rss" >> "$RESULTS"
}

normalise_dir() {
    local dir="$1"
    [[ "$dir" == */ ]] && printf "%s" "$dir" || printf "%s/" "$dir"
}

make_workdir() {
    mktemp -d
}

cleanup() {
    local dir="$1"
    if [[ "$KEEP" -eq 0 ]]; then
        rm -rf "$dir"
    else
        echo "Kept working directory: $dir"
    fi
}

# ─── Prereq check ─────────────────────────────────────────────────────────────

check_prereq() {
    [[ "$RUNS"    =~ ^[0-9]+$ && "$RUNS"    -gt 0 ]] || die "--runs must be a positive integer"
    [[ "$RECORDS" =~ ^[0-9]+$ && "$RECORDS" -gt 0 ]] || die "--records must be a positive integer"
    [[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -gt 0 ]] || die "--threads must be a positive integer"
    command -v samtools >/dev/null 2>&1 || die "samtools not in PATH"
    [[ -f "$RUST_BIN/bismark_methylation_extractor" ]] \
        || die "Rust binaries not built - run: cargo build --release --workspace"
    if [[ "$USE_TEST_FILES" -eq 1 ]]; then
        local fa="${CUSTOM_FASTA:-$TEST_FILES/NC_010473.fa.gz}"
        local fq1="${CUSTOM_FASTQ1:-$TEST_FILES/test_R1.fastq.gz}"
        local fq2="${CUSTOM_FASTQ2:-$TEST_FILES/test_R2.fastq.gz}"
        [[ -f "$fa"  ]] || die "FASTA not found: $fa"
        [[ -f "$fq1" ]] || die "FASTQ R1 not found: $fq1"
        [[ -f "$fq2" ]] || die "FASTQ R2 not found: $fq2"
        command -v bowtie2 >/dev/null 2>&1 \
            || die "bowtie2 not in PATH (required for --test-files)"
        command -v bowtie2-build >/dev/null 2>&1 \
            || die "bowtie2-build not in PATH (required for --test-files)"
    fi
}

# ─── Input generation ─────────────────────────────────────────────────────────

make_fake_aligner_dir() {
    local wd="$1"
    local bin_dir="$wd/fake_aligner"
    mkdir -p "$bin_dir"
    cat > "$bin_dir/bowtie2-build" << 'EOF'
#!/usr/bin/env bash
exit 0
EOF
    chmod +x "$bin_dir/bowtie2-build"
    echo "$bin_dir"
}

make_large_sam() {
    local out="$1" count="$2"
    perl -e '
        my ($out, $count) = @ARGV;
        open my $fh, ">", $out or die "open $out: $!";
        print $fh "\@HD\tVN:1.6\tSO:unsorted\n";
        print $fh "\@SQ\tSN:chr1\tLN:100000000\n";
        print $fh "\@PG\tID:Bismark\tPN:Bismark\tVN:v0.25.1\tCL:bismark --genome /g -1 r1.fq\n";
        my @xr = qw(CT CT GA GA);
        my @xg = qw(CT GA CT GA);
        for my $i (1..$count) {
            my $j = ($i - 1) % 4;
            my $flag = ($j == 1 || $j == 2) ? 16 : 0;
            my $pos = 100 + ($i * 17);
            print $fh join("\t",
                "read$i", $flag, "chr1", $pos, 255, "20M", "*", 0, 0,
                "ACGTACGTACGTACGTACGT", "IIIIIIIIIIIIIIIIIIII",
                "XM:Z:ZzXxHh..ZzXxHh..ZzXx", "XR:Z:$xr[$j]", "XG:Z:$xg[$j]"
            ), "\n";
        }
    ' "$out" "$count"
}

make_large_methylation() {
    local out="$1" count="$2"
    perl -e '
        my ($out, $count) = @ARGV;
        open my $fh, ">", $out or die "open $out: $!";
        for my $i (1..$count) {
            my $pos = 100 + $i;
            my $state = ($i % 2) ? "Z" : "z";
            my $strand = ($i % 2) ? "+" : "-";
            print $fh "read$i\t$strand\tchr1\t$pos\t$state\n";
        }
    ' "$out" "$count"
}

make_large_genome_and_coverage() {
    local genome_dir="$1" cov="$2" count="$3"
    mkdir -p "$genome_dir"
    perl -e '
        my ($genome_dir, $cov, $count) = @ARGV;
        my $len = $count + 1000;
        open my $fa, ">", "$genome_dir/chr1.fa" or die "open fasta: $!";
        print $fa ">chr1\n";
        my $chunk = "ACGT" x 250;
        my $written = 0;
        while ($written < $len) {
            my $take = $len - $written > length($chunk) ? length($chunk) : $len - $written;
            print $fa substr($chunk, 0, $take), "\n";
            $written += $take;
        }
        close $fa;

        open my $cfh, ">", $cov or die "open coverage: $!";
        for (my $pos = 2; $pos <= $count; $pos += 4) {
            print $cfh "chr1\t$pos\t$pos\t50.0\t1\t1\n";
        }
    ' "$genome_dir" "$cov" "$count"
}

prepare_inputs() {
    local wd="$1"
    make_large_sam "$wd/large.sam" "$RECORDS"
    samtools view -bS "$wd/large.sam" > "$wd/large.bam"
    samtools index "$wd/large.bam"
    make_large_methylation "$wd/CpG_OT_large.txt" "$RECORDS"
    make_large_genome_and_coverage "$wd/genome" "$wd/large.cov" "$RECORDS"
}

prepare_test_files_inputs() {
    local wd="$1"
    local fa="${CUSTOM_FASTA:-$TEST_FILES/NC_010473.fa.gz}"
    local fq1="${CUSTOM_FASTQ1:-$TEST_FILES/test_R1.fastq.gz}"
    local fq2="${CUSTOM_FASTQ2:-$TEST_FILES/test_R2.fastq.gz}"
    local genome_dir="$wd/test_files"
    mkdir -p "$genome_dir"
    cp "$fa"  "$genome_dir/"
    cp "$fq1" "$genome_dir/"
    cp "$fq2" "$genome_dir/"
    local fq1_base; fq1_base=$(basename "$fq1")
    local fq2_base; fq2_base=$(basename "$fq2")

    echo "Preparing genome..."
    (cd "$wd" && perl "$PERL_BIN/bismark_genome_preparation" "$genome_dir" \
        >/dev/null 2>"$wd/logs/genome_preparation.err")

    echo "Aligning paired-end FASTQs with Perl Bismark..."
    (cd "$wd" && perl "$PERL_BIN/bismark" \
        --genome "$genome_dir" \
        -1 "$genome_dir/$fq1_base" \
        -2 "$genome_dir/$fq2_base" \
        >/dev/null 2>"$wd/logs/bismark_align.err")

    local stem; stem=$(basename "$fq1_base" .gz); stem="${stem%.fastq}"; stem="${stem%.fq}"
    [[ -f "$wd/${stem}_bismark_bt2_pe.bam" ]] || die "Expected Bismark BAM not found"
}

# ─── Benchmark functions ───────────────────────────────────────────────────────

bench_genome_prep() {
    local wd="$1" run="$2"
    local fake_aligner="$3"
    local perl_genome="$wd/run${run}/genome_prep/perl_genome"
    local rust_genome="$wd/run${run}/genome_prep/rust_genome"
    mkdir -p "$perl_genome" "$rust_genome"
    cp "$wd/genome/chr1.fa" "$perl_genome/chr1.fa"
    cp "$wd/genome/chr1.fa" "$rust_genome/chr1.fa"

    local perl_cmd=(perl "$PERL_BIN/bismark_genome_preparation" --path_to_aligner "$fake_aligner")
    local rust_cmd=("$RUST_BIN/bismark_genome_preparation" --path_to_aligner "$fake_aligner")
    if [[ "$THREADS" -gt 1 ]]; then
        perl_cmd+=(--parallel "$THREADS")
        rust_cmd+=(--parallel "$THREADS")
    fi

    time_command "$wd/logs/genome_prep_perl_$run" "${perl_cmd[@]}" "$perl_genome"
    append_result "bismark_genome_preparation" "perl" "$run" "$LAST_ELAPSED"

    time_command "$wd/logs/genome_prep_rust_$run" "${rust_cmd[@]}" "$rust_genome"
    append_result "bismark_genome_preparation" "rust" "$run" "$LAST_ELAPSED"
}

bench_test_files_genome_prep() {
    local wd="$1" run="$2"
    local fake_aligner="$3"
    local perl_genome="$wd/run${run}/genome_prep/perl_genome"
    local rust_genome="$wd/run${run}/genome_prep/rust_genome"
    mkdir -p "$perl_genome" "$rust_genome"
    local fa="${CUSTOM_FASTA:-$TEST_FILES/NC_010473.fa.gz}"
    cp "$fa" "$perl_genome/"
    cp "$fa" "$rust_genome/"

    local perl_cmd=(perl "$PERL_BIN/bismark_genome_preparation" --path_to_aligner "$fake_aligner")
    local rust_cmd=("$RUST_BIN/bismark_genome_preparation" --path_to_aligner "$fake_aligner")
    if [[ "$THREADS" -gt 1 ]]; then
        perl_cmd+=(--parallel "$THREADS")
        rust_cmd+=(--parallel "$THREADS")
    fi

    time_command "$wd/logs/test_files_genome_prep_perl_$run" "${perl_cmd[@]}" "$perl_genome"
    append_result "test_files/bismark_genome_preparation" "perl" "$run" "$LAST_ELAPSED"

    time_command "$wd/logs/test_files_genome_prep_rust_$run" "${rust_cmd[@]}" "$rust_genome"
    append_result "test_files/bismark_genome_preparation" "rust" "$run" "$LAST_ELAPSED"
}

bench_extractor() {
    local wd="$1" run="$2"
    local perl_dir="$wd/run${run}/extractor/perl"
    local rust_dir="$wd/run${run}/extractor/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    time_command "$wd/logs/extractor_perl_$run" \
        perl "$PERL_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off --comprehensive --parallel "$THREADS" \
        --output "$perl_dir" "$wd/large.sam"
    append_result "bismark_methylation_extractor" "perl" "$run" "$LAST_ELAPSED"

    time_command "$wd/logs/extractor_rust_$run" \
        "$RUST_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off --comprehensive --parallel "$THREADS" \
        --dir "$rust_dir" "$wd/large.sam"
    append_result "bismark_methylation_extractor" "rust" "$run" "$LAST_ELAPSED"
}

bench_test_files_extractor() {
    local wd="$1" run="$2"
    local perl_dir="$wd/run${run}/extractor/perl"
    local rust_dir="$wd/run${run}/extractor/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    time_command "$wd/logs/test_files_extractor_perl_$run" \
        perl "$PERL_BIN/bismark_methylation_extractor" \
        --paired --no_header --mbias_off --comprehensive --parallel "$THREADS" \
        --output "$perl_dir" "$wd/test_R1_bismark_bt2_pe.bam"
    append_result "test_files/bismark_methylation_extractor" "perl" "$run" "$LAST_ELAPSED"

    time_command "$wd/logs/test_files_extractor_rust_$run" \
        "$RUST_BIN/bismark_methylation_extractor" \
        --paired --no_header --mbias_off --comprehensive --parallel "$THREADS" \
        --dir "$rust_dir" "$wd/test_R1_bismark_bt2_pe.bam"
    append_result "test_files/bismark_methylation_extractor" "rust" "$run" "$LAST_ELAPSED"
}

bench_dedup() {
    local wd="$1" run="$2"
    local perl_dir="$wd/run${run}/dedup/perl"
    local rust_dir="$wd/run${run}/dedup/rust"
    mkdir -p "$perl_dir" "$rust_dir"
    cp "$wd/large.bam"     "$perl_dir/large.bam"
    cp "$wd/large.bam.bai" "$perl_dir/large.bam.bai"
    cp "$wd/large.bam"     "$rust_dir/large.bam"
    cp "$wd/large.bam.bai" "$rust_dir/large.bam.bai"

    time_command "$wd/logs/dedup_perl_$run" \
        perl "$PERL_BIN/deduplicate_bismark" \
        --single --parallel "$THREADS" --output_dir "$perl_dir" "$perl_dir/large.bam"
    append_result "deduplicate_bismark" "perl" "$run" "$LAST_ELAPSED"

    time_command "$wd/logs/dedup_rust_$run" \
        "$RUST_BIN/deduplicate_bismark" \
        --single --parallel "$THREADS" --output_dir "$rust_dir" "$rust_dir/large.bam"
    append_result "deduplicate_bismark" "rust" "$run" "$LAST_ELAPSED"
}

bench_test_files_dedup() {
    local wd="$1" run="$2"
    local perl_dir="$wd/run${run}/dedup/perl"
    local rust_dir="$wd/run${run}/dedup/rust"
    mkdir -p "$perl_dir" "$rust_dir"
    cp "$wd/test_R1_bismark_bt2_pe.bam" "$perl_dir/test.bam"
    cp "$wd/test_R1_bismark_bt2_pe.bam" "$rust_dir/test.bam"

    time_command "$wd/logs/test_files_dedup_perl_$run" \
        perl "$PERL_BIN/deduplicate_bismark" \
        --paired --parallel "$THREADS" --output_dir "$perl_dir" "$perl_dir/test.bam"
    append_result "test_files/deduplicate_bismark" "perl" "$run" "$LAST_ELAPSED"

    time_command "$wd/logs/test_files_dedup_rust_$run" \
        "$RUST_BIN/deduplicate_bismark" \
        --paired --parallel "$THREADS" --output_dir "$rust_dir" "$rust_dir/test.bam"
    append_result "test_files/deduplicate_bismark" "rust" "$run" "$LAST_ELAPSED"
}

bench_bedgraph() {
    local wd="$1" run="$2"
    local perl_dir="$wd/run${run}/bedgraph/perl"
    local rust_dir="$wd/run${run}/bedgraph/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    time_command "$wd/logs/bedgraph_perl_$run" \
        perl "$PERL_BIN/bismark2bedGraph" \
        --output large.bedGraph --no_header --dir "$perl_dir" "$wd/CpG_OT_large.txt"
    append_result "bismark2bedGraph" "perl" "$run" "$LAST_ELAPSED"

    time_command "$wd/logs/bedgraph_rust_$run" \
        "$RUST_BIN/bismark2bedGraph" \
        --output large.bedGraph --no_header --dir "$rust_dir" "$wd/CpG_OT_large.txt"
    append_result "bismark2bedGraph" "rust" "$run" "$LAST_ELAPSED"
}

bench_test_files_bedgraph() {
    local wd="$1" run="$2"
    local extractor_dir="$wd/run${run}/extractor/rust"
    local cpg_file
    cpg_file="$(find "$extractor_dir" -maxdepth 1 -type f -name 'CpG_context_*.txt' | sort | head -1)"
    [[ -n "$cpg_file" ]] || die "No test_files CpG extractor output found for bedGraph benchmark"

    local perl_dir="$wd/run${run}/bedgraph/perl"
    local rust_dir="$wd/run${run}/bedgraph/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    time_command "$wd/logs/test_files_bedgraph_perl_$run" \
        perl "$PERL_BIN/bismark2bedGraph" \
        --output test_files.bedGraph --no_header --dir "$perl_dir" "$cpg_file"
    append_result "test_files/bismark2bedGraph" "perl" "$run" "$LAST_ELAPSED"

    time_command "$wd/logs/test_files_bedgraph_rust_$run" \
        "$RUST_BIN/bismark2bedGraph" \
        --output test_files.bedGraph --no_header --dir "$rust_dir" "$cpg_file"
    append_result "test_files/bismark2bedGraph" "rust" "$run" "$LAST_ELAPSED"
}

bench_coverage2cytosine() {
    local wd="$1" run="$2"
    local perl_dir="$wd/run${run}/coverage2cytosine/perl"
    local rust_dir="$wd/run${run}/coverage2cytosine/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    time_command "$wd/logs/coverage2cytosine_perl_$run" \
        perl "$PERL_BIN/coverage2cytosine" \
        --genome_folder "$wd/genome" \
        --output "$perl_dir/large.CpG_report.txt" "$wd/large.cov"
    append_result "coverage2cytosine" "perl" "$run" "$LAST_ELAPSED"

    time_command "$wd/logs/coverage2cytosine_rust_$run" \
        "$RUST_BIN/coverage2cytosine" \
        --genome_folder "$wd/genome" \
        --output "$rust_dir/large.CpG_report.txt" "$wd/large.cov"
    append_result "coverage2cytosine" "rust" "$run" "$LAST_ELAPSED"
}

bench_test_files_coverage2cytosine() {
    local wd="$1" run="$2"
    local perl_dir="$wd/run${run}/coverage2cytosine/perl"
    local rust_dir="$wd/run${run}/coverage2cytosine/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    local cov="$wd/run${run}/bedgraph/rust/test_files.bismark.cov.gz"
    [[ -f "$cov" ]] || die "No test_files coverage file found for coverage2cytosine benchmark"

    time_command "$wd/logs/test_files_coverage2cytosine_perl_$run" \
        perl "$PERL_BIN/coverage2cytosine" \
        --genome_folder "$wd/test_files" \
        --output "$perl_dir/test_files.CpG_report.txt" "$cov"
    append_result "test_files/coverage2cytosine" "perl" "$run" "$LAST_ELAPSED"

    time_command "$wd/logs/test_files_coverage2cytosine_rust_$run" \
        "$RUST_BIN/coverage2cytosine" \
        --genome_folder "$wd/test_files" \
        --output "$rust_dir/test_files.CpG_report.txt" "$cov"
    append_result "test_files/coverage2cytosine" "rust" "$run" "$LAST_ELAPSED"
}

# ─── Summary ──────────────────────────────────────────────────────────────────

print_summary() {
    local cases
    cases="$(awk -F, 'NR > 1 { seen[$1] = 1 } END { for (c in seen) print c }' \
             "$RESULTS" | sort)"
    local case_width
    case_width="$(awk -F, 'NR > 1 { if (length($1) > max) max = length($1) } \
                            END { print (max > 4 ? max : 4) }' "$RESULTS")"

    echo ""
    if [[ "$MEASURE_MEM" -eq 1 ]]; then
        echo "Performance + memory summary (wall-clock seconds; peak RSS averaged across runs)"
        printf "%-*s  %10s  %10s  %10s  %12s  %12s  %9s\n" \
            "$case_width" "case" \
            "perl avg" "rust avg" "speedup" \
            "perl RAM" "rust RAM" "RAM ratio"
        printf "%-*s  %10s  %10s  %10s  %12s  %12s  %9s\n" \
            "$case_width" "----" \
            "--------" "--------" "-------" \
            "--------" "--------" "---------"
    else
        echo "Performance summary (wall-clock seconds; lower is better)"
        printf "%-*s  %10s  %10s  %10s  %10s  %10s\n" \
            "$case_width" "case" \
            "perl avg" "rust avg" "speedup" "perl min" "rust min"
        printf "%-*s  %10s  %10s  %10s  %10s  %10s\n" \
            "$case_width" "----" \
            "--------" "--------" "-------" "--------" "--------"
    fi

    while IFS= read -r case_name; do
        [[ -n "$case_name" ]] || continue
        local perl_csv rust_csv
        perl_csv="$(mktemp)"
        rust_csv="$(mktemp)"
        awk -F, -v c="$case_name" '$1 == c && $2 == "perl" { print $0 }' \
            "$RESULTS" > "$perl_csv"
        awk -F, -v c="$case_name" '$1 == c && $2 == "rust" { print $0 }' \
            "$RESULTS" > "$rust_csv"

        local perl_avg rust_avg
        perl_avg="$(mean_csv "$perl_csv")"
        rust_avg="$(mean_csv "$rust_csv")"

        if [[ "$MEASURE_MEM" -eq 1 ]]; then
            local perl_rss rust_rss
            perl_rss="$(mean_rss_csv "$perl_csv")"
            rust_rss="$(mean_rss_csv "$rust_csv")"
            printf "%-*s  %10.3f  %10.3f  %10s  %12s  %12s  %9s\n" \
                "$case_width" "$case_name" \
                "$perl_avg" "$rust_avg" \
                "$(ratio "$perl_avg" "$rust_avg")" \
                "$(format_bytes "$perl_rss")" \
                "$(format_bytes "$rust_rss")" \
                "$(mem_ratio "$perl_rss" "$rust_rss")"
        else
            local perl_min rust_min
            perl_min="$(min_csv "$perl_csv")"
            rust_min="$(min_csv "$rust_csv")"
            printf "%-*s  %10.3f  %10.3f  %10s  %10.3f  %10.3f\n" \
                "$case_width" "$case_name" \
                "$perl_avg" "$rust_avg" \
                "$(ratio "$perl_avg" "$rust_avg")" \
                "$perl_min" "$rust_min"
        fi
        rm -f "$perl_csv" "$rust_csv"
    done <<< "$cases"

    echo ""
    if [[ "$KEEP" -eq 1 ]]; then
        echo "Raw results: $RESULTS"
    else
        echo "Raw results are kept only with --keep."
    fi
    if [[ "$MEASURE_MEM" -eq 1 ]]; then
        echo "RAM figures are mean peak RSS across $RUNS run(s)."
        if [[ "$TIME_SUPPORTS_OUTFILE" -eq 0 ]]; then
            echo "Note: command stderr was merged with timing output (gtime not available)."
        fi
    fi
}

# ─── Main ─────────────────────────────────────────────────────────────────────

check_prereq
setup_time_cmd

WD="$(make_workdir)"
RESULTS="$WD/results.csv"
mkdir -p "$WD/logs"
printf "case,implementation,seconds,run,rss_bytes\n" > "$RESULTS"
FAKE_ALIGNER="$(make_fake_aligner_dir "$WD")"

if [[ "$USE_TEST_FILES" -eq 1 ]]; then
    echo "Preparing test_files inputs..."
    prepare_test_files_inputs "$WD"
else
    echo "Preparing synthetic inputs ($RECORDS records)..."
    prepare_inputs "$WD"
fi

MEM_MSG=""
[[ "$MEASURE_MEM" -eq 1 ]] && MEM_MSG=" with RAM tracking (${TIME_CMD} ${TIME_ARGS})"
echo "Running $RUNS benchmark run(s) with $THREADS thread(s) where supported${MEM_MSG}..."

for run in $(seq 1 "$RUNS"); do
    echo ""
    echo "Run $run/$RUNS"
    if [[ "$USE_TEST_FILES" -eq 1 ]]; then
        bench_test_files_genome_prep "$WD" "$run" "$FAKE_ALIGNER"
        bench_test_files_extractor   "$WD" "$run"
        bench_test_files_dedup       "$WD" "$run"
        bench_test_files_bedgraph    "$WD" "$run"
        bench_test_files_coverage2cytosine "$WD" "$run"
    else
        bench_genome_prep      "$WD" "$run" "$FAKE_ALIGNER"
        bench_extractor        "$WD" "$run"
        bench_dedup            "$WD" "$run"
        bench_bedgraph         "$WD" "$run"
        bench_coverage2cytosine "$WD" "$run"
    fi
done

print_summary
cleanup "$WD"
