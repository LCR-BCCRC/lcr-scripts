#!/usr/bin/env bash

# This script runs some unit tests across the different file formats supported by cnv2igv.py
# to help reveal unexpected changes across versions

# Usage:
# In the same directory where the script resides, simply run 
# ./run_tests.sh > tests/output/run_tests.sh.log
# 
# IMPORTANT!
# First ensure completion of the script without errors. 
# After that, you should check that all the files in tests/output/ have fresh time stamps, indicating they were overwritten by the run
# Next, check for changes to their contents using git diff
# git diff tests/output/*
# The files should all be identical. 

in_dir="tests/input"
out_dir="tests/output"

echo "testing sequenza mode"
./cnv2igv.py  --sample "Sequenza_SAMPLE1" --mode sequenza $in_dir/sequenza_segments.txt > $out_dir/sequenza_segments.seg

echo "testing battenberg mode"
./cnv2igv.py --mode battenberg --sample SP116712 $in_dir/SP116712_subclones.txt > $out_dir/SP116712_subclones.seg
