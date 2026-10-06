#!/usr/bin/env bash
# Smoke tests for the sniffles container.
# Verifies that all required tools and python packages are present and loadable.
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

check "sniffles is available"            "command -v sniffles"
check "sniffles reports version 2.7.5"   "sniffles --version | grep -q 'Version 2.7.5'"
check "python package pysam loads"       "python -c 'import pysam'"

$PASS
