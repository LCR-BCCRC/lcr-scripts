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
./cnv2igv.py  --sample "Sequenza_SAMPLE1" --mode sequenza --preserve_log_ratio $in_dir/sequenza_segments.txt > $out_dir/sequenza_segments.preserved.seg


echo "testing purecn mode"
./cnv2igv.py  --mode purecn $in_dir/DLBCL10538T_dnacopy.seg > $out_dir/DLBCL10538T_dnacopy.seg
./cnv2igv.py  --mode purecn --preserve_log_ratio $in_dir/DLBCL10538T_dnacopy.seg > $out_dir/DLBCL10538T_dnacopy.preserved.seg

echo "testing sample_id override"
./cnv2igv.py  --mode purecn $in_dir/DLBCL-RICOVER_258-Tumor_dnacopy.seg > $out_dir/DLBCL-RICOVER_258-Tumor_dnacopy.seg
./cnv2igv.py  --mode purecn --sample DLBCL-RICOVER_258-Tumor $in_dir/DLBCL-RICOVER_258-Tumor_dnacopy.seg > $out_dir/DLBCL-RICOVER_258-Tumor_dnacopy_explicit_sample.seg


echo "testing battenberg mode"
./cnv2igv.py --mode battenberg $in_dir/SP116712_subclones.txt > $out_dir/SP116712_subclones.seg
./cnv2igv.py --preserve_log_ratio --mode battenberg $in_dir/SP116712_subclones.txt > $out_dir/SP116712_subclones.preserved.seg

