#!/usr/bin/env bash
# Smoke tests for the battenberg-base image: every tool and R package Battenberg needs is present,
# and Battenberg itself is not (test images install it on top).
# Usage: ./run_tests.sh

PASS=true

check() {
    local desc="$1"
    local cmd="$2"
    local output
    output=$(eval "$cmd" 2>&1)
    if [ $? -eq 0 ]; then
        echo "PASS: $desc"
    else
        echo "FAIL: $desc"
        echo "$output" | sed 's/^/  /'
        PASS=false
    fi
}

check "alleleCounter is available"   "command -v alleleCounter"
check "impute2 is available"         "command -v impute2"
check "Rscript is available"         "command -v Rscript"
check "R is 4.5.3"                   "Rscript --vanilla -e 'stopifnot(getRversion() == \"4.5.3\")'"
for pkg in ASCAT optparse remotes RColorBrewer ggplot2 readr gtools gridExtra doParallel foreach parallel splines; do
    check "R package $pkg loads"     "Rscript --vanilla -e 'suppressMessages(library($pkg))'"
done
check "Battenberg is not installed"  "Rscript --vanilla -e 'stopifnot(!requireNamespace(\"Battenberg\", quietly = TRUE))'"

$PASS
