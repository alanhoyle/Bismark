#!/usr/bin/env bash
# Performance comparison: run Perl and Rust versions of Bismark downstream
# tools on synthetic inputs and report wall-clock timing summaries.
#
# Usage:
#   ./rust/tests/performance.sh [--runs N] [--records N] [--threads N] [--keep] [--test-files]
#
# Requirements:
#   - samtools in PATH
#   - Perl scripts at ../  (relative to rust/)
#   - Rust binaries built: cargo build --release --workspace
#   - bowtie2/bowtie2-build in PATH when using --test-files
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
USE_TEST_FILES=0

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
        --keep)
            KEEP=1
            shift
            ;;
        --test-files)
            USE_TEST_FILES=1
            shift
            ;;
        -h|--help)
            sed -n '1,18p' "$0"
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 2
            ;;
    esac
done

die() { echo "FATAL: $*" >&2; exit 1; }

check_prereq() {
    [[ "$RUNS" =~ ^[0-9]+$ && "$RUNS" -gt 0 ]] || die "--runs must be a positive integer"
    [[ "$RECORDS" =~ ^[0-9]+$ && "$RECORDS" -gt 0 ]] || die "--records must be a positive integer"
    [[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -gt 0 ]] || die "--threads must be a positive integer"
    command -v samtools >/dev/null 2>&1 || die "samtools not in PATH"
    [[ -f "$RUST_BIN/bismark_methylation_extractor" ]] \
        || die "Rust binaries not built - run: cargo build --release --workspace"
    if [[ "$USE_TEST_FILES" -eq 1 ]]; then
        [[ -f "$TEST_FILES/NC_010473.fa.gz" ]] || die "NC_010473.fa.gz not found"
        [[ -f "$TEST_FILES/test_R1.fastq.gz" ]] || die "test_R1.fastq.gz not found"
        [[ -f "$TEST_FILES/test_R2.fastq.gz" ]] || die "test_R2.fastq.gz not found"
        command -v bowtie2 >/dev/null 2>&1 || die "bowtie2 not in PATH (required for --test-files)"
        command -v bowtie2-build >/dev/null 2>&1 || die "bowtie2-build not in PATH (required for --test-files)"
    fi
}

now_seconds() {
    perl -MTime::HiRes=time -e 'printf "%.6f\n", time'
}

elapsed_seconds() {
    local start="$1" end="$2"
    awk -v s="$start" -v e="$end" 'BEGIN { printf "%.6f", e - s }'
}

time_command() {
    local log_prefix="$1"
    shift
    local start end
    start="$(now_seconds)"
    "$@" >"${log_prefix}.out" 2>"${log_prefix}.err"
    end="$(now_seconds)"
    elapsed_seconds "$start" "$end"
}

ratio() {
    local perl_time="$1" rust_time="$2"
    awk -v p="$perl_time" -v r="$rust_time" 'BEGIN {
        if (r == 0) {
            printf "inf"
        } else {
            printf "%.2fx", p / r
        }
    }'
}

mean_csv() {
    local csv="$1"
    awk -F, '$3 ~ /^[0-9.]+$/ { sum += $3; n++ } END { if (n) printf "%.6f", sum / n; else printf "0.000000" }' "$csv"
}

min_csv() {
    local csv="$1"
    awk -F, '$3 ~ /^[0-9.]+$/ { if (n == 0 || $3 < min) min = $3; n++ } END { if (n) printf "%.6f", min; else printf "0.000000" }' "$csv"
}

append_result() {
    local case_name="$1" impl="$2" run="$3" seconds="$4"
    printf "%s,%s,%s,%s\n" "$case_name" "$impl" "$seconds" "$run" >> "$RESULTS"
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
    local genome_dir="$wd/test_files"
    mkdir -p "$genome_dir"
    cp "$TEST_FILES/NC_010473.fa.gz" "$genome_dir/"
    cp "$TEST_FILES/test_R1.fastq.gz" "$genome_dir/"
    cp "$TEST_FILES/test_R2.fastq.gz" "$genome_dir/"

    echo "Preparing copied test_files genome..."
    (cd "$wd" && perl "$PERL_BIN/bismark_genome_preparation" "$genome_dir" >/dev/null 2>"$wd/logs/genome_preparation.err")

    echo "Aligning test_files paired-end FASTQs with Perl Bismark..."
    (cd "$wd" && perl "$PERL_BIN/bismark" \
        --genome "$genome_dir" \
        -1 "$genome_dir/test_R1.fastq.gz" \
        -2 "$genome_dir/test_R2.fastq.gz" \
        >/dev/null 2>"$wd/logs/bismark_align.err")

    [[ -f "$wd/test_R1_bismark_bt2_pe.bam" ]] || die "Expected Bismark BAM not found"
}

bench_genome_prep() {
    local wd="$1" run="$2"
    local fake_aligner="$3"
    local perl_genome="$wd/run${run}/genome_prep/perl_genome"
    local rust_genome="$wd/run${run}/genome_prep/rust_genome"
    mkdir -p "$perl_genome" "$rust_genome"
    cp "$wd/genome/chr1.fa" "$perl_genome/chr1.fa"
    cp "$wd/genome/chr1.fa" "$rust_genome/chr1.fa"

    local t
    local perl_cmd=(perl "$PERL_BIN/bismark_genome_preparation" --path_to_aligner "$fake_aligner")
    local rust_cmd=("$RUST_BIN/bismark_genome_preparation" --path_to_aligner "$fake_aligner")
    if [[ "$THREADS" -gt 1 ]]; then
        perl_cmd+=(--parallel "$THREADS")
        rust_cmd+=(--parallel "$THREADS")
    fi

    t="$(time_command "$wd/logs/genome_prep_perl_$run" \
        "${perl_cmd[@]}" "$perl_genome")"
    append_result "bismark_genome_preparation" "perl" "$run" "$t"

    t="$(time_command "$wd/logs/genome_prep_rust_$run" \
        "${rust_cmd[@]}" "$rust_genome")"
    append_result "bismark_genome_preparation" "rust" "$run" "$t"
}

bench_test_files_genome_prep() {
    local wd="$1" run="$2"
    local fake_aligner="$3"
    local perl_genome="$wd/run${run}/genome_prep/perl_genome"
    local rust_genome="$wd/run${run}/genome_prep/rust_genome"
    mkdir -p "$perl_genome" "$rust_genome"
    cp "$TEST_FILES/NC_010473.fa.gz" "$perl_genome/"
    cp "$TEST_FILES/NC_010473.fa.gz" "$rust_genome/"

    local t
    local perl_cmd=(perl "$PERL_BIN/bismark_genome_preparation" --path_to_aligner "$fake_aligner")
    local rust_cmd=("$RUST_BIN/bismark_genome_preparation" --path_to_aligner "$fake_aligner")
    if [[ "$THREADS" -gt 1 ]]; then
        perl_cmd+=(--parallel "$THREADS")
        rust_cmd+=(--parallel "$THREADS")
    fi

    t="$(time_command "$wd/logs/test_files_genome_prep_perl_$run" \
        "${perl_cmd[@]}" "$perl_genome")"
    append_result "test_files/bismark_genome_preparation" "perl" "$run" "$t"

    t="$(time_command "$wd/logs/test_files_genome_prep_rust_$run" \
        "${rust_cmd[@]}" "$rust_genome")"
    append_result "test_files/bismark_genome_preparation" "rust" "$run" "$t"
}

bench_extractor() {
    local wd="$1" run="$2"
    local perl_dir="$wd/run${run}/extractor/perl"
    local rust_dir="$wd/run${run}/extractor/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    local t
    t="$(time_command "$wd/logs/extractor_perl_$run" \
        perl "$PERL_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off --comprehensive --parallel "$THREADS" \
        --output "$perl_dir" "$wd/large.sam")"
    append_result "bismark_methylation_extractor" "perl" "$run" "$t"

    t="$(time_command "$wd/logs/extractor_rust_$run" \
        "$RUST_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off --comprehensive --parallel "$THREADS" \
        --dir "$rust_dir" "$wd/large.sam")"
    append_result "bismark_methylation_extractor" "rust" "$run" "$t"
}

bench_test_files_extractor() {
    local wd="$1" run="$2"
    local perl_dir="$wd/run${run}/extractor/perl"
    local rust_dir="$wd/run${run}/extractor/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    local t
    t="$(time_command "$wd/logs/test_files_extractor_perl_$run" \
        perl "$PERL_BIN/bismark_methylation_extractor" \
        --paired --no_header --mbias_off --comprehensive --parallel "$THREADS" \
        --output "$perl_dir" "$wd/test_R1_bismark_bt2_pe.bam")"
    append_result "test_files/bismark_methylation_extractor" "perl" "$run" "$t"

    t="$(time_command "$wd/logs/test_files_extractor_rust_$run" \
        "$RUST_BIN/bismark_methylation_extractor" \
        --paired --no_header --mbias_off --comprehensive --parallel "$THREADS" \
        --dir "$rust_dir" "$wd/test_R1_bismark_bt2_pe.bam")"
    append_result "test_files/bismark_methylation_extractor" "rust" "$run" "$t"
}

bench_dedup() {
    local wd="$1" run="$2"
    local perl_dir="$wd/run${run}/dedup/perl"
    local rust_dir="$wd/run${run}/dedup/rust"
    mkdir -p "$perl_dir" "$rust_dir"
    cp "$wd/large.bam" "$perl_dir/large.bam"
    cp "$wd/large.bam.bai" "$perl_dir/large.bam.bai"
    cp "$wd/large.bam" "$rust_dir/large.bam"
    cp "$wd/large.bam.bai" "$rust_dir/large.bam.bai"

    local t
    t="$(time_command "$wd/logs/dedup_perl_$run" \
        perl "$PERL_BIN/deduplicate_bismark" \
        --single --output_dir "$perl_dir" "$perl_dir/large.bam")"
    append_result "deduplicate_bismark" "perl" "$run" "$t"

    t="$(time_command "$wd/logs/dedup_rust_$run" \
        "$RUST_BIN/deduplicate_bismark" \
        --single --output_dir "$rust_dir" "$rust_dir/large.bam")"
    append_result "deduplicate_bismark" "rust" "$run" "$t"
}

bench_test_files_dedup() {
    local wd="$1" run="$2"
    local perl_dir="$wd/run${run}/dedup/perl"
    local rust_dir="$wd/run${run}/dedup/rust"
    mkdir -p "$perl_dir" "$rust_dir"
    cp "$wd/test_R1_bismark_bt2_pe.bam" "$perl_dir/test.bam"
    cp "$wd/test_R1_bismark_bt2_pe.bam" "$rust_dir/test.bam"

    local t
    t="$(time_command "$wd/logs/test_files_dedup_perl_$run" \
        perl "$PERL_BIN/deduplicate_bismark" \
        --paired --output_dir "$perl_dir" "$perl_dir/test.bam")"
    append_result "test_files/deduplicate_bismark" "perl" "$run" "$t"

    t="$(time_command "$wd/logs/test_files_dedup_rust_$run" \
        "$RUST_BIN/deduplicate_bismark" \
        --paired --output_dir "$rust_dir" "$rust_dir/test.bam")"
    append_result "test_files/deduplicate_bismark" "rust" "$run" "$t"
}

bench_bedgraph() {
    local wd="$1" run="$2"
    local perl_dir="$wd/run${run}/bedgraph/perl"
    local rust_dir="$wd/run${run}/bedgraph/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    local t
    t="$(time_command "$wd/logs/bedgraph_perl_$run" \
        perl "$PERL_BIN/bismark2bedGraph" \
        --output large.bedGraph --no_header --dir "$perl_dir" "$wd/CpG_OT_large.txt")"
    append_result "bismark2bedGraph" "perl" "$run" "$t"

    t="$(time_command "$wd/logs/bedgraph_rust_$run" \
        "$RUST_BIN/bismark2bedGraph" \
        --output large.bedGraph --no_header --dir "$rust_dir" "$wd/CpG_OT_large.txt")"
    append_result "bismark2bedGraph" "rust" "$run" "$t"
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

    local t
    t="$(time_command "$wd/logs/test_files_bedgraph_perl_$run" \
        perl "$PERL_BIN/bismark2bedGraph" \
        --output test_files.bedGraph --no_header --dir "$perl_dir" "$cpg_file")"
    append_result "test_files/bismark2bedGraph" "perl" "$run" "$t"

    t="$(time_command "$wd/logs/test_files_bedgraph_rust_$run" \
        "$RUST_BIN/bismark2bedGraph" \
        --output test_files.bedGraph --no_header --dir "$rust_dir" "$cpg_file")"
    append_result "test_files/bismark2bedGraph" "rust" "$run" "$t"
}

bench_coverage2cytosine() {
    local wd="$1" run="$2"
    local perl_dir="$wd/run${run}/coverage2cytosine/perl"
    local rust_dir="$wd/run${run}/coverage2cytosine/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    local t
    t="$(time_command "$wd/logs/coverage2cytosine_perl_$run" \
        perl "$PERL_BIN/coverage2cytosine" \
        --genome_folder "$wd/genome" --output "$perl_dir/large.CpG_report.txt" "$wd/large.cov")"
    append_result "coverage2cytosine" "perl" "$run" "$t"

    t="$(time_command "$wd/logs/coverage2cytosine_rust_$run" \
        "$RUST_BIN/coverage2cytosine" \
        --genome_folder "$wd/genome" --output "$rust_dir/large.CpG_report.txt" "$wd/large.cov")"
    append_result "coverage2cytosine" "rust" "$run" "$t"
}

bench_test_files_coverage2cytosine() {
    local wd="$1" run="$2"
    local perl_dir="$wd/run${run}/coverage2cytosine/perl"
    local rust_dir="$wd/run${run}/coverage2cytosine/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    local cov="$wd/run${run}/bedgraph/rust/test_files.bismark.cov.gz"
    [[ -f "$cov" ]] || die "No test_files coverage file found for coverage2cytosine benchmark"

    local t
    t="$(time_command "$wd/logs/test_files_coverage2cytosine_perl_$run" \
        perl "$PERL_BIN/coverage2cytosine" \
        --genome_folder "$wd/test_files" --output "$perl_dir/test_files.CpG_report.txt" "$cov")"
    append_result "test_files/coverage2cytosine" "perl" "$run" "$t"

    t="$(time_command "$wd/logs/test_files_coverage2cytosine_rust_$run" \
        "$RUST_BIN/coverage2cytosine" \
        --genome_folder "$wd/test_files" --output "$rust_dir/test_files.CpG_report.txt" "$cov")"
    append_result "test_files/coverage2cytosine" "rust" "$run" "$t"
}

print_summary() {
    local cases
    cases="$(awk -F, 'NR > 1 { seen[$1] = 1 } END { for (c in seen) print c }' "$RESULTS" | sort)"
    local case_width
    case_width="$(awk -F, 'NR > 1 { if (length($1) > max) max = length($1) } END { print (max > 4 ? max : 4) }' "$RESULTS")"

    echo ""
    echo "Performance summary (seconds; lower is better)"
    printf "%-*s  %10s  %10s  %10s  %10s  %10s\n" "$case_width" "case" "perl avg" "rust avg" "speedup" "perl min" "rust min"
    printf "%-*s  %10s  %10s  %10s  %10s  %10s\n" "$case_width" "----" "--------" "--------" "-------" "--------" "--------"

    while IFS= read -r case_name; do
        [[ -n "$case_name" ]] || continue
        local perl_csv rust_csv perl_avg rust_avg perl_min rust_min
        perl_csv="$(mktemp)"
        rust_csv="$(mktemp)"
        awk -F, -v c="$case_name" '$1 == c && $2 == "perl" { print $0 }' "$RESULTS" > "$perl_csv"
        awk -F, -v c="$case_name" '$1 == c && $2 == "rust" { print $0 }' "$RESULTS" > "$rust_csv"
        perl_avg="$(mean_csv "$perl_csv")"
        rust_avg="$(mean_csv "$rust_csv")"
        perl_min="$(min_csv "$perl_csv")"
        rust_min="$(min_csv "$rust_csv")"
        printf "%-*s  %10.3f  %10.3f  %10s  %10.3f  %10.3f\n" \
            "$case_width" "$case_name" "$perl_avg" "$rust_avg" "$(ratio "$perl_avg" "$rust_avg")" "$perl_min" "$rust_min"
        rm -f "$perl_csv" "$rust_csv"
    done <<< "$cases"

    echo ""
    if [[ "$KEEP" -eq 1 ]]; then
        echo "Raw timings: $RESULTS"
    else
        echo "Raw timings are kept only with --keep."
    fi
}

check_prereq

WD="$(make_workdir)"
RESULTS="$WD/results.csv"
mkdir -p "$WD/logs"
printf "case,implementation,seconds,run\n" > "$RESULTS"
FAKE_ALIGNER="$(make_fake_aligner_dir "$WD")"

if [[ "$USE_TEST_FILES" -eq 1 ]]; then
    echo "Preparing test_files inputs..."
    prepare_test_files_inputs "$WD"
else
    echo "Preparing synthetic inputs ($RECORDS records)..."
    prepare_inputs "$WD"
fi

echo "Running $RUNS benchmark run(s) with $THREADS thread(s) where supported..."
for run in $(seq 1 "$RUNS"); do
    echo ""
    echo "Run $run/$RUNS"
    if [[ "$USE_TEST_FILES" -eq 1 ]]; then
        bench_test_files_genome_prep "$WD" "$run" "$FAKE_ALIGNER"
        bench_test_files_extractor "$WD" "$run"
        bench_test_files_dedup "$WD" "$run"
        bench_test_files_bedgraph "$WD" "$run"
        bench_test_files_coverage2cytosine "$WD" "$run"
    else
        bench_genome_prep "$WD" "$run" "$FAKE_ALIGNER"
        bench_extractor "$WD" "$run"
        bench_dedup "$WD" "$run"
        bench_bedgraph "$WD" "$run"
        bench_coverage2cytosine "$WD" "$run"
    fi
done

print_summary
cleanup "$WD"
