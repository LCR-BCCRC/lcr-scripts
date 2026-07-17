#!/usr/bin/env bash
# Smoke tests for the mfr container.
# Verifies that all required tools and R packages are present and loadable.
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

check "Rscript is available"        "command -v Rscript"
check "R package readr loads"       "Rscript --vanilla -e 'library(readr)'"
check "R package dplyr loads"       "Rscript --vanilla -e 'library(dplyr)'"
check "R package cluster loads"     "Rscript --vanilla -e 'library(cluster)'"
check "R package ggplot2 loads"     "Rscript --vanilla -e 'library(ggplot2)'"
check "python is available"         "command -v python"
check "tabix is available"          "command -v tabix"

$PASS
