#!/usr/bin/env bash
# Smoke tests for the MutationTimeR container.
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

check "Rscript is available"                                "command -v Rscript"
check "R package tidyverse loads"                           "Rscript --vanilla -e 'library(tidyverse)'"
check "R package data.table loads"                          "Rscript --vanilla -e 'library(data.table)'"
check "R package MutationTimeR loads"                       "Rscript --vanilla -e 'library(MutationTimeR)'"
check "R package BSgenome loads"                            "Rscript --vanilla -e 'library(BSgenome)'"
check "R package BSgenome.Hsapiens.UCSC.hg19 loads"         "Rscript --vanilla -e 'library(BSgenome.Hsapiens.UCSC.hg19)'"
check "R package BSgenome.Hsapiens.UCSC.hg38 loads"         "Rscript --vanilla -e 'library(BSgenome.Hsapiens.UCSC.hg38)'"
check "R package GAMBLR.helpers loads"                      "Rscript --vanilla -e 'library(GAMBLR.helpers)'"
check "R package ggpubr loads"                              "Rscript --vanilla -e 'library(ggpubr)'"

## TEST

$PASS