#!/usr/bin/env bash
set -euo pipefail

# Smoke test: verify vcf2maf.pl, samtools, and VEP are installed and in PATH.
# Usage: ./run_tests.sh [output_dir]
# output_dir is accepted for CI compatibility but no golden files are produced.

OUTDIR="${1:-tests/output}"
mkdir -p "$OUTDIR"

echo "Checking samtools..."
samtools --version | head -1

echo "Checking vcf2maf.pl..."
# vcf2maf.pl exits non-zero without required args but prints version/usage
vcf2maf.pl 2>&1 | head -5 || true
which vcf2maf.pl

echo "Checking VEP..."
vep --help 2>&1 | head -5 || true
which vep

echo "PASS: vcf2maf.pl, samtools, and VEP are all in PATH"
