#!/usr/bin/env bash
# Differential test: run Perl and Rust versions of Bismark downstream tools
# on the same input, then diff every output file byte-for-byte.
#
# Usage:
#   ./rust/tests/differential_test.sh [--keep] [--test-files]
#
# Requirements:
#   - samtools in PATH
#   - Perl scripts at ../  (relative to rust/)
#   - Rust binaries built: cargo build --release --workspace
#   - Test data at ../test_files/
#
# Options:
#   --keep   Keep temporary output directories on failure for inspection
#   --test-files
#            Use test_files/test_R1.fastq.gz and test_R2.fastq.gz by running
#            Perl Bismark first, then compare downstream Perl/Rust tools.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
RUST_BIN="$SCRIPT_DIR/../target/release"
PERL_BIN="$REPO_ROOT"
TEST_FILES="$REPO_ROOT/test_files"

KEEP=0
USE_TEST_FILES=0
for arg in "$@"; do
    case "$arg" in
        --keep) KEEP=1 ;;
        --test-files) USE_TEST_FILES=1 ;;
        -h|--help)
            sed -n '1,16p' "$0"
            exit 0
            ;;
        *)
            echo "Unknown option: $arg" >&2
            exit 2
            ;;
    esac
done

PASS=0; FAIL=0

# ─── Helpers ─────────────────────────────────────────────────────────────────

die() { echo "FATAL: $*" >&2; exit 1; }

check_prereq() {
    [[ -d "$TEST_FILES" ]] || die "test_files/ not found at $TEST_FILES"
    [[ -f "$TEST_FILES/NC_010473.fa.gz" ]] || die "NC_010473.fa.gz not found"
    [[ -f "$TEST_FILES/test_R1.fastq.gz" ]] || die "test_R1.fastq.gz not found"
    [[ -f "$TEST_FILES/test_R2.fastq.gz" ]] || die "test_R2.fastq.gz not found"
    command -v samtools >/dev/null 2>&1 || die "samtools not in PATH"
    [[ -f "$RUST_BIN/bismark_methylation_extractor" ]] \
        || die "Rust binaries not built — run: cargo build --release --workspace"
    if [[ "$USE_TEST_FILES" -eq 1 ]]; then
        command -v bowtie2 >/dev/null 2>&1 || die "bowtie2 not in PATH (required for --test-files)"
        command -v bowtie2-build >/dev/null 2>&1 || die "bowtie2-build not in PATH (required for --test-files)"
    fi
}

run_diff() {
    local name="$1" perl_file="$2" rust_file="$3"
    if diff -q "$perl_file" "$rust_file" >/dev/null 2>&1; then
        echo "  PASS  $name"
        PASS=$(( PASS + 1 ))
    else
        echo "  FAIL  $name"
        diff --unified=3 "$perl_file" "$rust_file" | head -30 || true
        FAIL=$(( FAIL + 1 ))
    fi
}

run_diff_no_paths() {
    # Compare two files after stripping lines that contain absolute file paths
    # (dedup reports embed the input BAM path, which differs by temp dir)
    local name="$1" perl_file="$2" rust_file="$3"
    local strip='Total number of alignments analysed in'
    if diff -q <(grep -v "$strip" "$perl_file") <(grep -v "$strip" "$rust_file") >/dev/null 2>&1; then
        echo "  PASS  $name"
        PASS=$(( PASS + 1 ))
    else
        echo "  FAIL  $name"
        diff --unified=3 <(grep -v "$strip" "$perl_file") <(grep -v "$strip" "$rust_file") | head -30 || true
        FAIL=$(( FAIL + 1 ))
    fi
}

run_diff_gz() {
    # Compare two gzip-compressed files by decompressing on the fly (gzip -dc for macOS compat)
    local name="$1" perl_file="$2" rust_file="$3"
    if diff -q <(gzip -dc "$perl_file") <(gzip -dc "$rust_file") >/dev/null 2>&1; then
        echo "  PASS  $name"
        PASS=$(( PASS + 1 ))
    else
        echo "  FAIL  $name"
        diff --unified=3 <(gzip -dc "$perl_file") <(gzip -dc "$rust_file") | head -30 || true
        FAIL=$(( FAIL + 1 ))
    fi
}

run_diff_sorted() {
    # For BAM outputs: sort SAM records before diffing (order may differ)
    local name="$1" perl_bam="$2" rust_bam="$3"
    local p_sorted r_sorted
    p_sorted=$(mktemp)
    r_sorted=$(mktemp)
    samtools view "$perl_bam" | sort > "$p_sorted"
    samtools view "$rust_bam" | sort > "$r_sorted"
    if diff -q "$p_sorted" "$r_sorted" >/dev/null 2>&1; then
        echo "  PASS  $name (sorted BAM)"
        PASS=$(( PASS + 1 ))
    else
        echo "  FAIL  $name (sorted BAM)"
        diff --unified=3 "$p_sorted" "$r_sorted" | head -30 || true
        FAIL=$(( FAIL + 1 ))
    fi
    rm -f "$p_sorted" "$r_sorted"
}

make_workdir() {
    local d; d=$(mktemp -d)
    echo "$d"
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

prepare_test_files_alignment() {
    local wd="$1"
    local genome_dir="$wd/test_files"
    mkdir -p "$genome_dir"
    cp "$TEST_FILES/NC_010473.fa.gz" "$genome_dir/"
    cp "$TEST_FILES/test_R1.fastq.gz" "$genome_dir/"
    cp "$TEST_FILES/test_R2.fastq.gz" "$genome_dir/"

    echo "  Preparing copied test_files genome..." >&2
    (cd "$wd" && perl "$PERL_BIN/bismark_genome_preparation" "$genome_dir" >/dev/null 2>"$wd/genome_preparation.err")

    echo "  Aligning test_files paired-end FASTQs with Perl Bismark..." >&2
    (cd "$wd" && perl "$PERL_BIN/bismark" \
        --genome "$genome_dir" \
        -1 "$genome_dir/test_R1.fastq.gz" \
        -2 "$genome_dir/test_R2.fastq.gz" \
        >/dev/null 2>"$wd/bismark_align.err")

    local bam="$wd/test_R1_bismark_bt2_pe.bam"
    [[ -f "$bam" ]] || die "Expected Bismark BAM not found at $bam"
    echo "$bam"
}

make_tiny_genome() {
    local dir="$1"
    mkdir -p "$dir"
    cat > "$dir/test.fa" << 'EOF'
>chr1 description
ACGTNacgtn
>chr2
CCCCGGGGAAAATTTTNNNN
EOF
}

cleanup() {
    local dir="$1"
    if [[ $KEEP -eq 0 || $FAIL -eq 0 ]]; then
        rm -rf "$dir"
    else
        echo "  (kept $dir for inspection)"
    fi
}

# ─── Build a small SAM file from FASTQ test data if a BAM is needed ──────────
# This section requires a Bismark-aligned BAM; use a pre-canned tiny SAM instead.

make_tiny_sam() {
    local out="$1"
    cat > "$out" << 'SAMEOF'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:10000
@PG	ID:Bismark	PN:Bismark	VN:v0.25.1	CL:bismark --genome /g -1 r1.fq
r1	0	chr1	100	255	20M	*	0	0	ACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIII	XM:Z:ZzXxHh..ZzXxHh..ZzXx	XR:Z:CT	XG:Z:CT
r2	0	chr1	200	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XM:Z:ZzXxHh....	XR:Z:CT	XG:Z:CT	NM:i:0
r3	16	chr1	300	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XM:Z:ZzXxHh....	XR:Z:CT	XG:Z:GA	NM:i:0
r4	16	chr1	400	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XM:Z:ZzXxHh....	XR:Z:GA	XG:Z:CT	NM:i:0
r5	0	chr1	500	255	10M	*	0	0	ACGTACGTAC	IIIIIIIIII	XM:Z:ZzXxHh....	XR:Z:GA	XG:Z:GA	NM:i:0
SAMEOF
}

make_pe_overlap_sam() {
    local out="$1"
    cat > "$out" << 'SAMEOF'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:10000
@PG	ID:Bismark	PN:Bismark	VN:v0.25.1	CL:bismark --genome /g -1 r1.fq -2 r2.fq
pair1/1	99	chr1	100	255	10M	=	105	15	ACGTACGTAC	IIIIIIIIII	XM:Z:ZzZzZzZzZz	XR:Z:CT	XG:Z:CT
pair1/2	147	chr1	105	255	10M	=	100	-15	ACGTACGTAC	IIIIIIIIII	XM:Z:ZzZzZzZzZz	XR:Z:CT	XG:Z:GA
SAMEOF
}

find_one() {
    local dir="$1" pattern="$2"
    find "$dir" -maxdepth 1 -type f -name "$pattern" | sort | head -1
}

compare_context_outputs() {
    local label="$1" perl_dir="$2" rust_dir="$3" gz="${4:-0}"
    local ctx strand perl_f rust_f
    for ctx in CpG CHG CHH; do
        perl_f=$(find_one "$perl_dir" "${ctx}_context_*.txt$([[ "$gz" == 1 ]] && echo .gz)")
        rust_f=$(find_one "$rust_dir" "${ctx}_context_*.txt$([[ "$gz" == 1 ]] && echo .gz)")
        if [[ -z "$perl_f" && -z "$rust_f" ]]; then
            continue
        elif [[ -n "$perl_f" && -n "$rust_f" ]]; then
            if [[ "$gz" == 1 ]]; then
                run_diff_gz "$label/$ctx" "$perl_f" "$rust_f"
            else
                run_diff "$label/$ctx" "$perl_f" "$rust_f"
            fi
        else
            echo "  FAIL  $label/$ctx (missing file: perl=$perl_f rust=$rust_f)"
            FAIL=$(( FAIL + 1 ))
        fi
    done
    for strand in OT CTOT CTOB OB; do
        for ctx in CpG CHG CHH; do
            perl_f=$(find_one "$perl_dir" "${ctx}_${strand}_*.txt")
            rust_f=$(find_one "$rust_dir" "${ctx}_${strand}_*.txt")
            if [[ -n "$perl_f" || -n "$rust_f" ]]; then
                if [[ -n "$perl_f" && -n "$rust_f" ]]; then
                    run_diff "$label/${ctx}_${strand}" "$perl_f" "$rust_f"
                else
                    echo "  FAIL  $label/${ctx}_${strand} (missing file: perl=$perl_f rust=$rust_f)"
                    FAIL=$(( FAIL + 1 ))
                fi
            fi
        done
    done
}

# ─── Test: bismark_methylation_extractor ─────────────────────────────────────

test_extractor() {
    echo ""
    echo "=== bismark_methylation_extractor ==="

    local wd; wd=$(make_workdir)
    local perl_dir="$wd/perl" rust_dir="$wd/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    local sam="$wd/test.sam"
    make_tiny_sam "$sam"

    # Perl
    perl "$PERL_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off --comprehensive \
        --output "$perl_dir" "$sam" 2>/dev/null

    # Rust
    "$RUST_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off --comprehensive \
        --dir "$rust_dir" "$sam" 2>/dev/null

    for ctx in CpG CHG CHH; do
        local perl_f rust_f
        perl_f=$(ls "$perl_dir/${ctx}_context_"*.txt 2>/dev/null | head -1)
        rust_f=$(ls "$rust_dir/${ctx}_context_"*.txt 2>/dev/null | head -1)
        if [[ -z "$perl_f" || -z "$rust_f" ]]; then
            echo "  SKIP  ${ctx}_context (file missing: perl=$perl_f rust=$rust_f)"
            continue
        fi
        run_diff "${ctx}_context" "$perl_f" "$rust_f"
    done

    cleanup "$wd"
}

test_extractor_strand_specific() {
    echo ""
    echo "=== bismark_methylation_extractor strand-specific ==="

    local wd; wd=$(make_workdir)
    local perl_dir="$wd/perl" rust_dir="$wd/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    local sam="$wd/test.sam"
    make_tiny_sam "$sam"

    perl "$PERL_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off \
        --output "$perl_dir" "$sam" 2>/dev/null

    "$RUST_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off \
        --dir "$rust_dir" "$sam" 2>/dev/null

    compare_context_outputs "strand_specific" "$perl_dir" "$rust_dir"
    cleanup "$wd"
}

test_extractor_modes() {
    echo ""
    echo "=== bismark_methylation_extractor modes ==="

    local wd; wd=$(make_workdir)
    local sam="$wd/test.sam"
    make_tiny_sam "$sam"

    local perl_dir="$wd/perl_merge" rust_dir="$wd/rust_merge"
    mkdir -p "$perl_dir" "$rust_dir"
    perl "$PERL_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off --comprehensive --merge_non_CpG \
        --output "$perl_dir" "$sam" 2>/dev/null
    "$RUST_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off --comprehensive --merge_non_CpG \
        --dir "$rust_dir" "$sam" 2>/dev/null
    run_diff "extractor/merge_CpG" "$(find_one "$perl_dir" "CpG_context_*.txt")" "$(find_one "$rust_dir" "CpG_context_*.txt")"
    run_diff "extractor/merge_Non_CpG" "$(find_one "$perl_dir" "Non_CpG_context_*.txt")" "$(find_one "$rust_dir" "Non_CpG_context_*.txt")"

    perl_dir="$wd/perl_gzip"; rust_dir="$wd/rust_gzip"
    mkdir -p "$perl_dir" "$rust_dir"
    perl "$PERL_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off --comprehensive --gzip \
        --output "$perl_dir" "$sam" 2>/dev/null
    "$RUST_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off --comprehensive --gzip \
        --dir "$rust_dir" "$sam" 2>/dev/null
    compare_context_outputs "gzip_comprehensive" "$perl_dir" "$rust_dir" 1

    perl_dir="$wd/perl_yacht"; rust_dir="$wd/rust_yacht"
    mkdir -p "$perl_dir" "$rust_dir"
    perl "$PERL_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off --yacht \
        --output "$perl_dir" "$sam" 2>/dev/null
    "$RUST_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off --yacht \
        --dir "$rust_dir" "$sam" 2>/dev/null
    run_diff "extractor/yacht" "$(find_one "$perl_dir" "any_C_context_*.txt")" "$(find_one "$rust_dir" "any_C_context_*.txt")"

    perl_dir="$wd/perl_ignore"; rust_dir="$wd/rust_ignore"
    mkdir -p "$perl_dir" "$rust_dir"
    perl "$PERL_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off --comprehensive --ignore 2 --ignore_3prime 1 \
        --output "$perl_dir" "$sam" 2>/dev/null
    "$RUST_BIN/bismark_methylation_extractor" \
        --single --no_header --mbias_off --comprehensive --ignore 2 --ignore_3prime 1 \
        --dir "$rust_dir" "$sam" 2>/dev/null
    compare_context_outputs "ignore" "$perl_dir" "$rust_dir"

    cleanup "$wd"
}

test_extractor_paired_overlap() {
    echo ""
    echo "=== bismark_methylation_extractor paired overlap ==="

    local wd; wd=$(make_workdir)
    local sam="$wd/pe.sam"
    make_pe_overlap_sam "$sam"

    local perl_dir="$wd/perl_default" rust_dir="$wd/rust_default"
    mkdir -p "$perl_dir" "$rust_dir"
    perl "$PERL_BIN/bismark_methylation_extractor" \
        --paired --no_header --mbias_off --comprehensive \
        --output "$perl_dir" "$sam" 2>/dev/null
    "$RUST_BIN/bismark_methylation_extractor" \
        --paired --no_header --mbias_off --comprehensive \
        --dir "$rust_dir" "$sam" 2>/dev/null
    run_diff "extractor/paired_default_no_overlap" "$(find_one "$perl_dir" "CpG_context_*.txt")" "$(find_one "$rust_dir" "CpG_context_*.txt")"

    perl_dir="$wd/perl_include" rust_dir="$wd/rust_include"
    mkdir -p "$perl_dir" "$rust_dir"
    perl "$PERL_BIN/bismark_methylation_extractor" \
        --paired --include_overlap --no_header --mbias_off --comprehensive \
        --output "$perl_dir" "$sam" 2>/dev/null
    "$RUST_BIN/bismark_methylation_extractor" \
        --paired --include_overlap --no_header --mbias_off --comprehensive \
        --dir "$rust_dir" "$sam" 2>/dev/null
    run_diff "extractor/paired_include_overlap" "$(find_one "$perl_dir" "CpG_context_*.txt")" "$(find_one "$rust_dir" "CpG_context_*.txt")"

    cleanup "$wd"
}

# ─── Test: deduplicate_bismark ────────────────────────────────────────────────

test_dedup() {
    echo ""
    echo "=== deduplicate_bismark ==="

    local wd; wd=$(make_workdir)
    local sam="$wd/test.sam"
    make_tiny_sam "$sam"

    # Convert to BAM
    local bam="$wd/test.bam"
    samtools view -bS "$sam" > "$bam"
    samtools index "$bam"

    local perl_dir="$wd/perl" rust_dir="$wd/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    cp "$bam" "$bam.bai" "$perl_dir/" 2>/dev/null || cp "$bam" "$perl_dir/"
    cp "$bam" "$bam.bai" "$rust_dir/" 2>/dev/null || cp "$bam" "$rust_dir/"

    # Perl
    perl "$PERL_BIN/deduplicate_bismark" \
        --single --output_dir "$perl_dir" "$perl_dir/test.bam" 2>/dev/null

    # Rust
    "$RUST_BIN/deduplicate_bismark" \
        --single --output_dir "$rust_dir" "$rust_dir/test.bam" 2>/dev/null

    local perl_rep rust_rep
    perl_rep=$(ls "$perl_dir/"*deduplication_report* 2>/dev/null | head -1)
    rust_rep=$(ls "$rust_dir/"*deduplication_report* 2>/dev/null | head -1)
    if [[ -n "$perl_rep" && -n "$rust_rep" ]]; then
        run_diff_no_paths "dedup_report" "$perl_rep" "$rust_rep"
    else
        echo "  SKIP  dedup_report (Perl or Rust report not found)"
    fi

    local perl_bam rust_bam
    perl_bam=$(ls "$perl_dir/"*deduplicated.bam 2>/dev/null | head -1)
    rust_bam=$(ls "$rust_dir/"*deduplicated.bam 2>/dev/null | head -1)
    if [[ -n "$perl_bam" && -n "$rust_bam" ]]; then
        run_diff_sorted "dedup_bam" "$perl_bam" "$rust_bam"
    else
        echo "  SKIP  dedup_bam (file not found)"
    fi

    cleanup "$wd"
}

# ─── Test: bismark2bedGraph ───────────────────────────────────────────────────

test_bedgraph() {
    echo ""
    echo "=== bismark2bedGraph ==="

    local wd; wd=$(make_workdir)

    # Create a minimal CpG_OT methylation extractor output
    local meth_file="$wd/CpG_OT_test.txt"
    cat > "$meth_file" << 'EOF'
r1	+	chr1	100	Z
r1	-	chr1	101	z
r2	+	chr1	200	Z
r2	-	chr1	201	z
r2	+	chr1	205	Z
EOF

    local perl_dir="$wd/perl" rust_dir="$wd/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    # Perl
    perl "$PERL_BIN/bismark2bedGraph" \
        --output test.bedGraph \
        --no_header \
        --dir "$perl_dir" \
        "$meth_file" 2>/dev/null

    # Rust
    "$RUST_BIN/bismark2bedGraph" \
        --output test.bedGraph \
        --no_header \
        --dir "$rust_dir" \
        "$meth_file" 2>/dev/null

    # Both Perl and Rust write gzip-compressed output files
    for f in test.bedGraph.gz test.bismark.cov.gz; do
        local perl_f="$perl_dir/$f" rust_f="$rust_dir/$f"
        if [[ -f "$perl_f" && -f "$rust_f" ]]; then
            run_diff_gz "bedGraph/$f" "$perl_f" "$rust_f"
        else
            echo "  SKIP  bedGraph/$f (file missing: perl=$(test -f "$perl_f" && echo y || echo n) rust=$(test -f "$rust_f" && echo y || echo n))"
        fi
    done

    cleanup "$wd"
}

# ─── Test: coverage2cytosine ─────────────────────────────────────────────────

test_coverage2cytosine() {
    echo ""
    echo "=== coverage2cytosine ==="

    local wd; wd=$(make_workdir)

    # Tiny genome
    local genome_dir="$wd/genome"
    mkdir -p "$genome_dir"
    printf ">chr1\nACGTCGACGTCGACGT\n" > "$genome_dir/chr1.fa"

    # Coverage input (bismark2bedGraph .bismark.cov format)
    local cov="$wd/sample.cov"
    printf "chr1\t5\t5\t80.0\t4\t1\nchr1\t11\t11\t50.0\t2\t2\n" > "$cov"

    local perl_dir="$wd/perl" rust_dir="$wd/rust"
    mkdir -p "$perl_dir" "$rust_dir"

    # Perl
    perl "$PERL_BIN/coverage2cytosine" \
        --genome_folder "$genome_dir" \
        --output "$perl_dir/sample.CpG_report.txt" \
        "$cov" 2>/dev/null

    # Rust
    "$RUST_BIN/coverage2cytosine" \
        --genome_folder "$genome_dir" \
        --output "$rust_dir/sample.CpG_report.txt" \
        "$cov" 2>/dev/null

    local perl_f="$perl_dir/sample.CpG_report.txt"
    local rust_f="$rust_dir/sample.CpG_report.txt"
    if [[ -f "$perl_f" && -f "$rust_f" ]]; then
        run_diff "CpG_report" "$perl_f" "$rust_f"
    else
        echo "  SKIP  CpG_report (Perl=$perl_f exists=$(test -f $perl_f && echo y || echo n); Rust=$(test -f $rust_f && echo y || echo n))"
    fi

    cleanup "$wd"
}

test_genome_preparation() {
    echo ""
    echo "=== bismark_genome_preparation ==="

    local wd; wd=$(make_workdir)
    local fake_aligner; fake_aligner=$(make_fake_aligner_dir "$wd")
    local perl_genome="$wd/perl_genome" rust_genome="$wd/rust_genome"

    make_tiny_genome "$perl_genome"
    mkdir -p "$rust_genome"
    cp "$perl_genome/test.fa" "$rust_genome/test.fa"

    perl "$PERL_BIN/bismark_genome_preparation" \
        --path_to_aligner "$fake_aligner" "$perl_genome" >/dev/null 2>"$wd/perl_genome_prep.err"

    "$RUST_BIN/bismark_genome_preparation" \
        --path_to_aligner "$fake_aligner" "$rust_genome" >/dev/null 2>"$wd/rust_genome_prep.err"

    run_diff "genome_prep/CT_conversion" \
        "$perl_genome/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa" \
        "$rust_genome/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"

    run_diff "genome_prep/GA_conversion" \
        "$perl_genome/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa" \
        "$rust_genome/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"

    cleanup "$wd"
}

test_genome_preparation_test_files() {
    echo ""
    echo "=== bismark_genome_preparation test_files ==="

    local wd; wd=$(make_workdir)
    local fake_aligner; fake_aligner=$(make_fake_aligner_dir "$wd")
    local perl_genome="$wd/perl_genome" rust_genome="$wd/rust_genome"
    mkdir -p "$perl_genome" "$rust_genome"
    cp "$TEST_FILES/NC_010473.fa.gz" "$perl_genome/"
    cp "$TEST_FILES/NC_010473.fa.gz" "$rust_genome/"

    perl "$PERL_BIN/bismark_genome_preparation" \
        --path_to_aligner "$fake_aligner" "$perl_genome" >/dev/null 2>"$wd/perl_genome_prep.err"

    "$RUST_BIN/bismark_genome_preparation" \
        --path_to_aligner "$fake_aligner" "$rust_genome" >/dev/null 2>"$wd/rust_genome_prep.err"

    run_diff "test_files/genome_prep_CT" \
        "$perl_genome/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa" \
        "$rust_genome/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"

    run_diff "test_files/genome_prep_GA" \
        "$perl_genome/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa" \
        "$rust_genome/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"

    cleanup "$wd"
}

test_test_files_inputs() {
    echo ""
    echo "=== test_files FASTQ-derived downstream parity ==="

    local wd; wd=$(make_workdir)
    local bam
    bam=$(prepare_test_files_alignment "$wd")

    local perl_dir="$wd/extractor_perl" rust_dir="$wd/extractor_rust"
    mkdir -p "$perl_dir" "$rust_dir"
    perl "$PERL_BIN/bismark_methylation_extractor" \
        --paired --no_header --mbias_off --comprehensive \
        --output "$perl_dir" "$bam" 2>/dev/null
    "$RUST_BIN/bismark_methylation_extractor" \
        --paired --no_header --mbias_off --comprehensive \
        --dir "$rust_dir" "$bam" 2>/dev/null
    compare_context_outputs "test_files/extractor" "$perl_dir" "$rust_dir"

    perl_dir="$wd/dedup_perl"; rust_dir="$wd/dedup_rust"
    mkdir -p "$perl_dir" "$rust_dir"
    cp "$bam" "$perl_dir/test.bam"
    cp "$bam" "$rust_dir/test.bam"
    perl "$PERL_BIN/deduplicate_bismark" \
        --paired --output_dir "$perl_dir" "$perl_dir/test.bam" 2>/dev/null
    "$RUST_BIN/deduplicate_bismark" \
        --paired --output_dir "$rust_dir" "$rust_dir/test.bam" 2>/dev/null
    local perl_rep rust_rep perl_bam rust_bam
    perl_rep=$(find_one "$perl_dir" "*deduplication_report*")
    rust_rep=$(find_one "$rust_dir" "*deduplication_report*")
    run_diff_no_paths "test_files/dedup_report" "$perl_rep" "$rust_rep"
    perl_bam=$(find_one "$perl_dir" "*deduplicated.bam")
    rust_bam=$(find_one "$rust_dir" "*deduplicated.bam")
    run_diff_sorted "test_files/dedup_bam" "$perl_bam" "$rust_bam"

    perl_dir="$wd/bedgraph_perl"; rust_dir="$wd/bedgraph_rust"
    mkdir -p "$perl_dir" "$rust_dir"
    perl "$PERL_BIN/bismark2bedGraph" \
        --output test_files.bedGraph --no_header --dir "$perl_dir" \
        "$(find_one "$wd/extractor_perl" "CpG_context_*.txt")" 2>/dev/null
    "$RUST_BIN/bismark2bedGraph" \
        --output test_files.bedGraph --no_header --dir "$rust_dir" \
        "$(find_one "$wd/extractor_rust" "CpG_context_*.txt")" 2>/dev/null
    run_diff_gz "test_files/bedGraph" "$perl_dir/test_files.bedGraph.gz" "$rust_dir/test_files.bedGraph.gz"
    run_diff_gz "test_files/coverage" "$perl_dir/test_files.bismark.cov.gz" "$rust_dir/test_files.bismark.cov.gz"

    perl_dir="$wd/cytosine_perl"; rust_dir="$wd/cytosine_rust"
    mkdir -p "$perl_dir" "$rust_dir"
    perl "$PERL_BIN/coverage2cytosine" \
        --genome_folder "$wd/test_files" \
        --output "$perl_dir/test_files.CpG_report.txt" \
        "$perl_dir/../bedgraph_perl/test_files.bismark.cov.gz" 2>/dev/null
    "$RUST_BIN/coverage2cytosine" \
        --genome_folder "$wd/test_files" \
        --output "$rust_dir/test_files.CpG_report.txt" \
        "$rust_dir/../bedgraph_rust/test_files.bismark.cov.gz" 2>/dev/null
    run_diff "test_files/CpG_report" "$perl_dir/test_files.CpG_report.txt" "$rust_dir/test_files.CpG_report.txt"

    cleanup "$wd"
}

# ─── Main ─────────────────────────────────────────────────────────────────────

check_prereq

if [[ "$USE_TEST_FILES" -eq 1 ]]; then
    test_genome_preparation_test_files
    test_test_files_inputs
else
    test_genome_preparation
    test_extractor
    test_extractor_strand_specific
    test_extractor_modes
    test_extractor_paired_overlap
    test_dedup
    test_bedgraph
    test_coverage2cytosine
fi

echo ""
echo "═══════════════════════════════════════"
echo " Results: $PASS passed, $FAIL failed"
echo "═══════════════════════════════════════"

[[ $FAIL -eq 0 ]] && exit 0 || exit 1
