#!/usr/bin/env bash
here="tests"
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

# test SEG mode lifting both directions
../liftover.sh SEG input/test.grch37.seg output/test.grch37.to_hg38.seg $grch37_chain NO
../liftover.sh SEG input/test.hg38.seg output/test.hg38.to_grch37.seg $hg38_chain NO
../liftover.sh SEG input/test.grch37.withheader.seg output/test.grch37.withheader.to_hg38.seg $grch37_chain YES
../liftover.sh SEG input/test.hg38.withheader.seg output/test.hg38.withheader.to_grch37.seg $hg38_chain YES
