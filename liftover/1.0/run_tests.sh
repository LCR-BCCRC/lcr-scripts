#!/usr/bin/env bash
# This script runs some unit tests across the different file formats supported by liftover.sh
# to help reveal unexpected changes across versions
# Usage:
# In the same directory where the script resides, simply run
# ./run_tests.sh > tests/output/run_tests.sh.log
# Or pass an output directory as the first argument (used by CI):
# ./run_tests.sh tests/docker_output
#
# IMPORTANT!
# First ensure completion of the script without errors.
# After that, you should check that all the files in tests/output/ have fresh time stamps, indicating they were overwritten by the run
# Next, check for changes to their contents using git diff
# git diff tests/output/*
# The files should all be identical.

OUTDIR="${1:-tests/output}"
mkdir -p "$OUTDIR"

hg38_chain="data/hg38ToHg19.over.chain.gz"
grch37_chain="data/hg19ToHg38.over.chain.gz"
for chain in "$hg38_chain" "$grch37_chain"; do
    if [ ! -f "$chain" ]; then
        echo "missing chain file: $chain" >&2
        exit 1
    fi
done

report_change () {
    echo "$1: $2"
    echo "$3: $4"
    echo "$(echo "100*($2 - $4)/$2"|bc)% change"
}


# test SEG mode lifting both directions (no-header mode)
IN="tests/input/test.grch37.seg"
OUT="$OUTDIR/test.grch37.to_hg38.seg"
./liftover.sh SEG $IN $OUT $grch37_chain NO 0.95
IN="tests/input/test.hg38.seg"
OUT="$OUTDIR/test.hg38.to_grch37.seg"
./liftover.sh SEG $IN $OUT $hg38_chain NO 0.95


IN="tests/input/test.grch37.withheader.seg"
OUT="$OUTDIR/test.grch37.withheader.to_hg38.seg"
./liftover.sh SEG $IN $OUT $grch37_chain YES 0.95

IN_SIZE=`cat $IN | perl -ane '$size = $F[3]-$F[2];$total+=$size;END{print "$total\n";}'`
OUT_SIZE=`cat $OUT | perl -ane '$size = $F[3]-$F[2];$total+=$size;END{print "$total\n";}'`

report_change $IN $IN_SIZE $OUT $OUT_SIZE

IN="tests/input/test.hg38.withheader.seg"
OUT="$OUTDIR/test.hg38.withheader.to_grch37.seg"
./liftover.sh SEG $IN $OUT $hg38_chain YES 0.95

IN_SIZE=`cat $IN | perl -ane '$size = $F[3]-$F[2];$total+=$size;END{print "$total\n";}'`
OUT_SIZE=`cat $OUT | perl -ane '$size = $F[3]-$F[2];$total+=$size;END{print "$total\n";}'`

report_change $IN $IN_SIZE $OUT $OUT_SIZE

IN="tests/input/cnvkit_hg38.seg"
OUT="$OUTDIR/cnvkit_hg38.to_grch37.seg"
./liftover.sh SEG $IN $OUT $hg38_chain YES 0.95

IN_SIZE=`cat $IN | perl -ane '$size = $F[3]-$F[2];$total+=$size;END{print "$total\n";}'`
OUT_SIZE=`cat $OUT | perl -ane '$size = $F[3]-$F[2];$total+=$size;END{print "$total\n";}'`

report_change $IN $IN_SIZE $OUT $OUT_SIZE

IN="tests/input/test.grch37.bed"
OUT="$OUTDIR/test.grch37.to_hg38.bed"
./liftover.sh BED $IN $OUT $grch37_chain NO 0.95

IN_SIZE=`cat $IN | perl -ane '$size = $F[2]-$F[1];$total+=$size;END{print "$total\n";}'`
OUT_SIZE=`cat $OUT | perl -ane '$size = $F[2]-$F[1];$total+=$size;END{print "$total\n";}'`

report_change $IN $IN_SIZE $OUT $OUT_SIZE

#reverse the process
./liftover.sh BED $OUT $OUTDIR/test.grch37.roundtrip.bed $hg38_chain NO 0.95

IN="tests/input/battenberg_hg38_with_header.bed"
OUT="$OUTDIR/battenberg_hg38_with_header.to_grch37.bed"
./liftover.sh BED $IN $OUT $hg38_chain YES 0.95
