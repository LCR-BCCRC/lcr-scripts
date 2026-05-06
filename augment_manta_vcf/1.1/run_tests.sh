#!/usr/bin/env bash
# Run unit tests for augment_manta_vcf.py to check for unexpected changes.
# Usage: run from the directory containing this script:
#   ./run_tests.sh
#
# After the first run, commit tests/output/ to establish the baseline.
# On subsequent runs, verify nothing changed with:
#   git diff tests/output/
#
# Variable header lines (##cmdline, ##regions_bed) will always differ
# across environments — compare only the variant records if doing a
# strict diff:
#   grep -v "^##cmdline\|^##regions_bed" tests/output/file.vcf | diff - <expected>

set -euo pipefail

SCRIPT="./augment_manta_vcf.py"
IN="tests/input"
OUT="tests/output"

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
    --bed_regions ${IN}/regions.bed \
    ${IN}/somaticSV.vcf \
    ${OUT}/somaticSV.with_regions.augmented.vcf

echo "All tests completed successfully."
echo "Verify outputs are unchanged: git diff tests/output/"
