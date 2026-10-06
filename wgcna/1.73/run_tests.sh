#!/usr/bin/env bash
# Smoke tests for the WGCNA container.
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

check "Rscript is available"             "command -v Rscript"
check "R package dplyr loads"            "Rscript --vanilla -e 'library(dplyr)'"
check "R package tidyr loads"            "Rscript --vanilla -e 'library(tidyr)'"
check "R package tibble loads"           "Rscript --vanilla -e 'library(tibble)'"
check "R package readr loads"            "Rscript --vanilla -e 'library(readr)'"
check "R package ggplot2 loads"          "Rscript --vanilla -e 'library(ggplot2)'"
check "R package matrixStats loads"      "Rscript --vanilla -e 'library(matrixStats)'"
check "R package DESeq2 loads"           "Rscript --vanilla -e 'library(DESeq2)'"
check "R package limma loads"            "Rscript --vanilla -e 'library(limma)'"
check "R package BiocParallel loads"     "Rscript --vanilla -e 'library(BiocParallel)'"
check "R package WGCNA loads"            "Rscript --vanilla -e 'library(WGCNA)'"
check "R package GAMBLR.helpers loads"   "Rscript --vanilla -e 'library(GAMBLR.helpers)'"

$PASS
