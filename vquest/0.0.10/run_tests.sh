#!/usr/bin/env bash
set -euo pipefail

# Run vquest with test sequences and compare output against committed golden files.
#
# Usage:
#   ./run_tests.sh                      # writes to tests/output (local golden-file generation)
#   ./run_tests.sh tests/docker_output  # used by CI to write Docker outputs for comparison
#
# Generating / updating golden files:
#   Run ./run_tests.sh from within the conda environment and commit tests/output/.
#   Note: IMGT occasionally updates their reference database, which can change
#   vquest_airr.tsv results. Re-run and re-commit tests/output/ when that happens.
#
# Note on moleculeType:
#   IMGT added a moleculeType field to their API. This is passed via the YAML config
#   file (tests/input/*.yml) rather than as a CLI flag, since vquest 0.0.10 does not
#   expose it as a command-line option but does forward unknown YAML keys to the API.

OUTDIR="${1:-tests/output}"
mkdir -p "$OUTDIR"

echo "Testing human IGH analysis..."
vquest \
    --outdir "$OUTDIR" \
    --fileSequences tests/input/human_igh.fasta \
    tests/input/human_igh.yml
