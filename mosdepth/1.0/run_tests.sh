#!/usr/bin/env bash
# Smoke tests for the mosdepth container.
# Verifies that all required tools are present and loadable.
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

check "mosdepth is available"            "command -v mosdepth"
check "mosdepth reports version 0.3.13"  "mosdepth --version | grep -q '0.3.13'"

$PASS
