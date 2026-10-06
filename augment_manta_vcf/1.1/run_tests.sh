#!/usr/bin/env bash
# Run unit tests for augment_manta_vcf.py to check for unexpected changes.
#
# Usage: run from the directory containing this script:
#   ./run_tests.sh [output_dir]
#
# output_dir defaults to tests/output (used for generating the golden baseline).
# Pass a different directory to compare against the committed baseline without
# overwriting it (used by CI to validate the Docker container):
#   ./run_tests.sh tests/docker_output
#
# Generating the golden baseline:
#   Activate the conda environment for this script, then run:
#     ./run_tests.sh
#   Inspect the outputs, then commit tests/output/ as the reference.
#
# On subsequent local runs, verify nothing changed with:
#   git diff tests/output/
#
# Note: ##cmdline and ##regions_bed header lines embed absolute paths and
# will always differ across environments. The CI comparison filters them out.

set -euo pipefail

SCRIPT="./augment_manta_vcf.py"
IN="tests/input"
OUT="${1:-tests/output}"
mkdir -p "$OUT"

echo "Test 1: somaticSV — augment fields and rename sample IDs"
python $SCRIPT \
    --tumour_id TUMOUR_SAMPLE \
    --normal_id NORMAL_SAMPLE \
    ${IN}/somaticSV.vcf \
    ${OUT}/somaticSV.augmented.vcf

echo "Test 2: tumorSV — augment fields and rename sample ID"
python $SCRIPT \
    --tumour_id TUMOUR_SAMPLE \
    ${IN}/tumorSV.vcf \
    ${OUT}/tumorSV.augmented.vcf

echo "Test 3: somaticSV with BED regions — verify REGIONS INFO field"
python $SCRIPT \
    --tumour_id TUMOUR_SAMPLE \
    --normal_id NORMAL_SAMPLE \
    ${IN}/somaticSV.vcf \
    ${OUT}/somaticSV.with_regions.augmented.vcf \
    --bed_regions ${IN}/regions.bed

echo "All tests completed successfully."
echo "Verify outputs are unchanged: git diff tests/output/"
