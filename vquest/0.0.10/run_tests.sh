#!/usr/bin/env bash
set -euo pipefail

# Smoke test: verify vquest is installed and produces valid AIRR output.
# Does not compare against golden files because the IMGT database is updated
# externally and results change without any change to this tool or its
# installation.

OUTDIR="${1:-tests/output}"
mkdir -p "$OUTDIR"

echo "Running vquest on test sequences..."
vquest \
    --outdir "$OUTDIR" \
    --fileSequences tests/input/human_igh.fasta \
    tests/input/human_igh.yml

AIRR="$OUTDIR/vquest_airr.tsv"
if [ ! -f "$AIRR" ]; then
    echo "FAIL: vquest_airr.tsv not produced"
    exit 1
fi

if ! head -1 "$AIRR" | grep -q "sequence_id"; then
    echo "FAIL: vquest_airr.tsv missing expected AIRR header"
    exit 1
fi

DATA_ROWS=$(awk 'NR>1 && NF' "$AIRR" | wc -l)
if [ "$DATA_ROWS" -lt 2 ]; then
    echo "FAIL: expected at least 2 data rows in vquest_airr.tsv, got $DATA_ROWS"
    exit 1
fi

echo "PASS: vquest produced valid AIRR output with $DATA_ROWS sequences"
