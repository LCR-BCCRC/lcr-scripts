#!/usr/bin/env bash
# Smoke tests for the battenberg container.
# Verifies that all required tools and R packages are present and loadable.
# Usage: ./run_tests.sh

PASS=true

check() {
    local desc="$1"
    local cmd="$2"
    if eval "$cmd" > /dev/null 2>&1; then
        echo "PASS: $desc"
    else
        echo "FAIL: $desc"
        PASS=false
    fi
}

check "alleleCounter is available"       "command -v alleleCounter"
check "impute2 is available"             "command -v impute2"
check "Rscript is available"             "command -v Rscript"
check "R package Battenberg loads"       "Rscript --vanilla -q -e 'library(Battenberg)'"
check "R package ASCAT loads"            "Rscript --vanilla -q -e 'library(ASCAT)'"
check "R package optparse loads"         "Rscript --vanilla -q -e 'library(optparse)'"

$PASS
