#!/usr/bin/env bash
set -euo pipefail

# Run lofreq against synthetic two-chromosome test data and compare outputs
# against committed golden files.
#
# Usage:
#   ./run_tests.sh                      # writes to tests/output (golden-file generation)
#   ./run_tests.sh tests/docker_output  # used by CI to write Docker outputs for comparison
#
# Generating / updating golden files:
#   Run this script from within the container and commit tests/output/.
#   Example:
#     docker run --rm \
#       -v /path/to/lofreq/2.1.5-post:/scripts \
#       <image> bash /scripts/run_tests.sh
#
# Test data:
#   tests/input/ref.fa      - two 200 bp synthetic chromosomes (ACGT×50, GCTA×50)
#   tests/input/reads.sam   - 50 reads per chromosome; 15/50 carry a 30% AF SNV
#                             at position 125: chr1 A>T, chr2 G>A
#
# What is tested:
#   1. lofreq call (single-threaded) correctly calls both SNVs
#   2. lofreq call-parallel splits by chromosome, merges with bcftools concat
#      (the VCF-merging fix from PR #109), and produces identical output

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${1:-$SCRIPT_DIR/tests/output}"
mkdir -p "$OUTDIR"

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# lofreq call emits AF=0.300000 while bcftools (used by call-parallel) normalises
# to AF=0.3. Strip trailing zeros so both paths produce identical output.
normalize_af() {
    sed -E 's/AF=([0-9]+\.[0-9]*[1-9])0+/AF=\1/g; s/AF=([0-9]+)\.0+/AF=\1/g'
}

# ── Build BAM ────────────────────────────────────────────────────────────────
samtools sort -O BAM -o "$TMPDIR/reads.bam" "$SCRIPT_DIR/tests/input/reads.sam"
samtools index "$TMPDIR/reads.bam"
lofreq indelqual --dindel \
    -f "$SCRIPT_DIR/tests/input/ref.fa" \
    -o "$TMPDIR/reads_iq.bam" \
    "$TMPDIR/reads.bam"
samtools index "$TMPDIR/reads_iq.bam"

# ── Single-threaded call ─────────────────────────────────────────────────────
echo "Testing single-threaded lofreq call..."
lofreq call \
    -f "$SCRIPT_DIR/tests/input/ref.fa" \
    -o "$TMPDIR/single.vcf" \
    "$TMPDIR/reads_iq.bam"
grep -v '^##' "$TMPDIR/single.vcf" | normalize_af > "$OUTDIR/single.vcf"

# ── Parallel call ─────────────────────────────────────────────────────────────
echo "Testing parallel lofreq call (bcftools concat path from PR #109)..."
lofreq call-parallel --pp-threads 2 \
    -f "$SCRIPT_DIR/tests/input/ref.fa" \
    -o "$TMPDIR/parallel.vcf.gz" \
    "$TMPDIR/reads_iq.bam"
zgrep -v '^##' "$TMPDIR/parallel.vcf.gz" | normalize_af > "$OUTDIR/parallel.vcf"

# ── Internal consistency check ───────────────────────────────────────────────
diff "$OUTDIR/single.vcf" "$OUTDIR/parallel.vcf" \
    && echo "PASS: single and parallel outputs agree" \
    || { echo "FAIL: single and parallel outputs differ"; exit 1; }
