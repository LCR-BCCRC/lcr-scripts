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
#   3. lofreq somatic --continue succeeds when normal preprocessing files are
#      already present in the output directory (regression test for the
#      morinlab fork's --continue mode fix)

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

# ── lofreq somatic --continue with pre-existing normal files ──────────────────
# Regression test for the morinlab fork's --continue mode fix: the CSB5
# upstream tried to parse test counts from the normal log file (where every
# line has a "stderr: " prefix), causing ValueError and a silent exit.
# The morinlab fork skips normal processing entirely in --continue mode.
echo "Testing lofreq somatic --continue with pre-existing normal files..."

NORMAL_DIR="$TMPDIR/normal_preprocessing"
mkdir -p "$NORMAL_DIR"

# Run normal-only preprocessing (same BAM for both -t and -n, as the
# pipeline does for unmatched normals via _lofreq_preprocess_normal)
lofreq somatic --normal_only --threads 2 \
    -t "$TMPDIR/reads_iq.bam" \
    -n "$TMPDIR/reads_iq.bam" \
    -f "$SCRIPT_DIR/tests/input/ref.fa" \
    -o "$NORMAL_DIR/" \
    || { echo "FAIL: lofreq somatic --normal_only"; exit 1; }

# lofreq filter does not auto-index its output; index any stringent VCFs
# that are missing a TBI so that lofreq call's -S argument can seek into them
for vcf in "$NORMAL_DIR"/normal_stringent.*.vcf.gz; do
    [ -f "${vcf}.tbi" ] || bcftools index -t "$vcf"
done

# Populate the --continue output directory with the normal files already in
# place, mirroring what _lofreq_link_to_preprocessed does via symlinks
CONTINUE_DIR="$TMPDIR/somatic_continue"
mkdir -p "$CONTINUE_DIR"
for f in normal_relaxed.vcf.gz normal_relaxed.vcf.gz.tbi normal_relaxed.log \
          normal_stringent.snvs.vcf.gz normal_stringent.snvs.vcf.gz.tbi \
          normal_stringent.indels.vcf.gz normal_stringent.indels.vcf.gz.tbi; do
    ln -s "$NORMAL_DIR/$f" "$CONTINUE_DIR/$f"
done

# Run --continue; the morinlab fork must skip normal calling and proceed
# directly to tumour calling without attempting to parse the normal log file
lofreq somatic --continue --threads 2 \
    -t "$TMPDIR/reads_iq.bam" \
    -n "$TMPDIR/reads_iq.bam" \
    -f "$SCRIPT_DIR/tests/input/ref.fa" \
    -o "$CONTINUE_DIR/" \
    || { echo "FAIL: lofreq somatic --continue with pre-existing normal files"; exit 1; }

for f in somatic_final.snvs.vcf.gz somatic_final.indels.vcf.gz; do
    [ -f "$CONTINUE_DIR/$f" ] \
        || { echo "FAIL: missing $f after --continue run"; exit 1; }
done
echo "PASS: lofreq somatic --continue handled pre-existing normal files correctly"
