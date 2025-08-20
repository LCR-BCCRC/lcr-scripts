#!/usr/bin/env bash
# This script runs some unit tests across the different file formats supported by liftover.sh
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

hg38_chain="hg38ToHg19.over.chain.gz"
grch37_chain="hg19ToHg38.over.chain.gz"
if [ ! -f $hg38_chain ]; then
    echo "downloading hg38ToHg19.over.chain.gz"
    curl -o hg38ToHg19.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz 
fi

if [ ! -f $grch37_chain ]; then
    echo "downloading hg19ToHg38.over.chain.gz"
    curl -o hg19ToHg38.over.chain.gz https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
fi

report_change () {
    echo "$1: $2"
    echo "$3: $4"
    echo "$(echo "100*($2 - $4)/$2"|bc)% change"
}


# test SEG mode lifting both directions
#../liftover.sh SEG input/test.grch37.seg output/test.grch37.to_hg38.seg $grch37_chain NO
#../liftover.sh SEG input/test.hg38.seg output/test.hg38.to_grch37.seg $hg38_chain NO
IN="tests/input/test.grch37.withheader.seg"
OUT="tests/output/test.grch37.withheader.to_hg38.seg"
./liftover.sh SEG $IN $OUT $grch37_chain YES 0.95

IN_SIZE=`cat $IN | perl -ane '$size = $F[3]-$F[2];$total+=$size;END{print "$total\n";}'`
OUT_SIZE=`cat $OUT | perl -ane '$size = $F[3]-$F[2];$total+=$size;END{print "$total\n";}'`

report_change $IN $IN_SIZE $OUT $OUT_SIZE

IN="tests/input/test.hg38.withheader.seg"
OUT="tests/output/test.hg38.withheader.to_grch37.seg"

./liftover.sh SEG $IN $OUT $hg38_chain YES 0.95
IN_SIZE=`cat $IN | perl -ane '$size = $F[3]-$F[2];$total+=$size;END{print "$total\n";}'`
OUT_SIZE=`cat $OUT | perl -ane '$size = $F[3]-$F[2];$total+=$size;END{print "$total\n";}'`

report_change $IN $IN_SIZE $OUT $OUT_SIZE
IN="tests/input/test.grch37.bed"
OUT="tests/output/test.grch37.to_hg38.bed"
./liftover.sh BED $IN $OUT $grch37_chain NO 0.95

IN_SIZE=`cat $IN | perl -ane '$size = $F[2]-$F[1];$total+=$size;END{print "$total\n";}'`
OUT_SIZE=`cat $OUT | perl -ane '$size = $F[2]-$F[1];$total+=$size;END{print "$total\n";}'`
report_change $IN $IN_SIZE $OUT $OUT_SIZE

./liftover.sh BED tests/output/test.grch37.to_hg38.bed tests/output/test.grch37.roundtrip.bed $hg38_chain NO 0.95
