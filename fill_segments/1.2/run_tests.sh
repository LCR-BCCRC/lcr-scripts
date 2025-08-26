#!/usr/bin/env bash
# This script runs some unit tests across the different file formats supported by fill_segments.sh
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

# conda activate /projects/rmorin_scratch/conda_environments/d4151b5f (bedtools 2.29.2)

hg38_chrArms="src/chromArm.hg38.bed"
grch37_chrArms="src/chromArm.grch37.bed"

hg38_blacklist="src/blacklisted.hg38.bed"
grch37_blacklist="src/blacklisted.grch37.bed"

in_dir="tests/input"
out_dir="tests/output"

echo "testing seg mode with purecn input"
./fill_segments.sh $hg38_chrArms $in_dir/DLBCL10538T_dnacopy.seg $hg38_blacklist $out_dir/DLBCL10538T_dnacopy.filled.seg DLBCL10538T SEG

echo "testing seg mode with purecn preserved input"
./fill_segments.sh $hg38_chrArms $in_dir/DLBCL10538T_dnacopy.preserved.seg $hg38_blacklist $out_dir/DLBCL10538T_dnacopy.preserved.filled.seg DLBCL10538T SEG

echo "testing seg mode with battenberg input"
./fill_segments.sh $hg38_chrArms $in_dir/SP116712_subclones.seg $hg38_blacklist $out_dir/SP116712_subclones.filled.seg SP116712 SEG

echo "testing seg mode with battenberg preserved input"
./fill_segments.sh $hg38_chrArms $in_dir/SP116712_subclones.preserved.seg $hg38_blacklist $out_dir/SP116712_subclones.preserved.filled.seg SP116712 SEG

echo "testing subclones mode"
./fill_segments.sh $hg38_chrArms $in_dir/SP116712_subclones.txt $hg38_blacklist $out_dir/SP116712_subclones.filled.txt SP116712 subclones

echo "testing sequenza mode"
./fill_segments.sh $hg38_chrArms $in_dir/sequenza_segments.txt $hg38_blacklist $out_dir/sequenza_segments.filled.txt Sequenza_SAMPLE1 sequenza
