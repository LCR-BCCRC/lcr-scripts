#!/usr/bin/env bash

set -euf -o pipefail

INPUT="$1"
OUTPUT="$2"

EBV_CHROM=${EBV_CHROM:-chrEBV}

echo -e "sample	total_reads	total_EBV	fraction_EBV" >> $OUTPUT

while IFS=$'\t' read -r SAMPLE BAM;do


    TOTAL_COV=$(samtools idxstats $BAM | grep -v "*" | datamash sum 3)
    EBV_COV=$(samtools view -c ${BAM} ${EBV_CHROM})
    FRACTION=$(echo -e $EBV_COV/$TOTAL_COV | bc -l)
    
    echo -e "${SAMPLE}	${TOTAL_COV}	${EBV_COV}	${FRACTION}" >> $OUTPUT 
done < $INPUT


