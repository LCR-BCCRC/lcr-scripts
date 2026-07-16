#!/usr/bin/env bash
# Smoke tests for the DriveMuTR container.
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

check "Rscript is available"                        "command -v Rscript"
check "R package dplyr loads"                        "Rscript --vanilla -e 'library(dplyr)'"
check "R package tidyr loads"                        "Rscript --vanilla -e 'library(tidyr)'"
check "R package tibble loads"                       "Rscript --vanilla -e 'library(tibble)'"
check "R package purrr loads"                        "Rscript --vanilla -e 'library(purrr)'"
check "R package readr loads"                        "Rscript --vanilla -e 'library(readr)'"
check "R package stringr loads"                      "Rscript --vanilla -e 'library(stringr)'"
check "R package ggplot2 loads"                      "Rscript --vanilla -e 'library(ggplot2)'"
check "R package rlang loads"                        "Rscript --vanilla -e 'library(rlang)'"
check "R package broom loads"                        "Rscript --vanilla -e 'library(broom)'"
check "R package ggbeeswarm loads"                   "Rscript --vanilla -e 'library(ggbeeswarm)'"
check "R package iml loads"                          "Rscript --vanilla -e 'library(iml)'"
check "R package pre loads"                          "Rscript --vanilla -e 'library(pre)'"
check "R package BiocParallel loads"                 "Rscript --vanilla -e 'library(BiocParallel)'"
check "R package Biostrings loads"                   "Rscript --vanilla -e 'library(Biostrings)'"
check "R package GenomicRanges loads"                "Rscript --vanilla -e 'library(GenomicRanges)'"
check "R package IRanges loads"                      "Rscript --vanilla -e 'library(IRanges)'"
check "R package rtracklayer loads"                  "Rscript --vanilla -e 'library(rtracklayer)'"
check "R package MotifDb loads"                      "Rscript --vanilla -e 'library(MotifDb)'"
check "R package motifbreakR loads"                  "Rscript --vanilla -e 'library(motifbreakR)'"
check "R package BSgenome.Hsapiens.UCSC.hg19 loads"  "Rscript --vanilla -e 'library(BSgenome.Hsapiens.UCSC.hg19)'"
check "R package GAMBLR.data loads"                  "Rscript --vanilla -e 'library(GAMBLR.data)'"

$PASS