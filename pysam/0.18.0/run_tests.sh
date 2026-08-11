#!/usr/bin/env bash
set -euo pipefail

# Smoke test: verify pysam, pandas, and oncopipe are installed with expected versions.
# Usage: ./run_tests.sh [output_dir]
# output_dir is accepted for CI compatibility but no golden files are produced.

OUTDIR="${1:-tests/output}"
mkdir -p "$OUTDIR"

echo "Testing Python package imports..."
python3 - <<'EOF'
import pysam
import pandas
import oncopipe

parts = tuple(int(x) for x in pysam.__version__.split(".")[:2])
assert parts >= (0, 18), f"Expected pysam >=0.18.0, got {pysam.__version__}"
print(f"pysam {pysam.__version__}: OK")
print(f"pandas {pandas.__version__}: OK")

# Verify op.absolute_symlink is available (used by augment_ssm.py)
assert hasattr(oncopipe, 'absolute_symlink'), \
    "oncopipe.absolute_symlink not found"
print(f"oncopipe (absolute_symlink): OK")
EOF

echo "PASS: pysam, pandas, and oncopipe are all importable"
