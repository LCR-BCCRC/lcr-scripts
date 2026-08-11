#!/usr/bin/env bash
set -euo pipefail

# Smoke test: verify GAMBLR and its key dependencies load correctly in R.
# Does not test database connectivity (requires live credentials).
# Usage: ./run_tests.sh [output_dir]
# output_dir is accepted for CI compatibility but no golden files are produced.

OUTDIR="${1:-tests/output}"
mkdir -p "$OUTDIR"

echo "Testing R GAMBLR package..."
Rscript - <<'EOF'
suppressPackageStartupMessages({
    library(GAMBLR)
    library(data.table)
    library(tidyverse)
    library(glue)
})

# Verify functions called by deblacklist_ssm.R are exported
stopifnot(
    "fread_maf not exported from GAMBLR" =
        exists("fread_maf", where = "package:GAMBLR", mode = "function"),
    "annotate_ssm_blacklist not exported from GAMBLR" =
        exists("annotate_ssm_blacklist", where = "package:GAMBLR", mode = "function")
)

cat(sprintf("GAMBLR loaded (%d exported symbols)\n", length(ls("package:GAMBLR"))))
cat("data.table: OK\n")
cat("tidyverse: OK\n")
cat("glue: OK\n")
EOF

echo "PASS: GAMBLR environment is functional"
