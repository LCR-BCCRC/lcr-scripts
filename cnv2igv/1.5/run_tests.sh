#!/usr/bin/env bash

in_dir="tests/input"
out_dir="tests/output"

echo "testing sequenza mode"
./cnv2igv.py  --sample "Sequenza_SAMPLE1" --mode sequenza $in_dir/sequenza_segments.txt > $out_dir/sequenza_segments.seg
./cnv2igv.py  --sample "Sequenza_SAMPLE1" --mode sequenza --preserve_log_ratio $in_dir/sequenza_segments.txt > $out_dir/sequenza_segments.preserved.seg


echo "testing purecn mode"
./cnv2igv.py  --mode purecn $in_dir/DLBCL10538T_dnacopy.seg > $out_dir/DLBCL10538T_dnacopy.seg
./cnv2igv.py  --mode purecn --preserve_log_ratio $in_dir/DLBCL10538T_dnacopy.seg > $out_dir/DLBCL10538T_dnacopy.preserved.seg

echo "testing battenberg mode"
./cnv2igv.py --preserve_log_ratio --mode battenberg $in_dir/SP116712_subclones.txt > $out_dir/SP116712_subclones.seg
./cnv2igv.py --preserve_log_ratio --mode battenberg $in_dir/SP116712_subclones.txt > $out_dir/SP116712_subclones.preserved.seg