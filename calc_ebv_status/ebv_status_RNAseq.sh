#!/usr/bin/env bash

set -euf -o pipefail

INPUT="$1"
OUTPUT="$2"

EBV_CHROM=${EBV_CHROM:-chrEBV}
echo -e "sample	total_EBV	EBER1	EBER2	total_EBER" >> $OUTPUT

while IFS=$'\t' read -r SAMPLE BAM;do

    EBV_COV=$(samtools view -c ${BAM} ${EBV_CHROM})
    EBER1_COV=$(samtools view -c ${BAM} ${EBV_CHROM}:6629-6795)
    EBER2_COV=$(samtools view -c ${BAM} ${EBV_CHROM}:6956-7128)
    TOTAL=$((EBER1_COV + EBER2_COV))

    echo -e "${SAMPLE}	${EBV_COV}	${EBER1_COV}	${EBER2_COV}	${TOTAL}" >> $OUTPUT 
done < $INPUT

