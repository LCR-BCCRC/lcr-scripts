#!/usr/bin/env bash

# This script will use bedtools subtract to leverage chromosome coordinates and seg file and fill all missing segments.
# Important! The chromosome prefix will be normalize that of the input arm bed file.
# This script is intended to be used within snakemake and expects to have environment with bedtools available.

# Usage:
# bash fill_segments.sh <input chromosome arm bed file> <input seg file> <input blacklisted bed file> <output seg file> <sample_id>
# Eexample:
# bash fill_segments.sh src/chromArm.hg19.bed TCRBOA7-T-WEX-test--matched.igv.seg TCRBOA7-T-WEX-test--matched.igv.filled.seg TCRBOA7-T-WEX
# 

set -e

# Read variables to store the arguments from command line
ARM_BED_PATH="$1"
SEG_PATH="$2"
BLACKLIST_PATH="$3"
RESULTS_PATH="$4"
THIS_SAMPLE_ID="$5"


# bedtools subtract doesn't like headers - so we need to strip it and store separately in temp file
cat $SEG_PATH | head -1 > $RESULTS_PATH.header

# Rearrange columns in the seg file to follow bed convention of chr-start-end
cat $SEG_PATH | grep -v "start" | perl -pale 'BEGIN { $"="\t"; } $_ = "@F[1..$#F,0]"' > $RESULTS_PATH.headerless.bed

# If there is inconsistent chr prefixing between SEG file and bed file of chromosome coordinates,
# normalize it to match between the two files.
ARM_IS_PREFIXED=$(head -1 $ARM_BED_PATH | cut -f 1)
SEG_IS_PREFIXED=$(head -1 $RESULTS_PATH.headerless.bed | cut -f 1)

# if the arm coordinates file is prefixed, but the SEG file is not
if [[ "$ARM_IS_PREFIXED" == *"chr"* && ! "$SEG_IS_PREFIXED" == *"chr"* ]]; then
    echo "Normalizing chr prefix in SEG file to match that of provided arm file..."
    # add chr prefix
    cat $RESULTS_PATH.headerless.bed | perl -lane '@a=split;$a[0] =~ s/^/chr/;print join "\t", @a;' > $RESULTS_PATH.normalized.bed
    # replace the headerless fiel with this normalized one
    cat $RESULTS_PATH.normalized.bed > $RESULTS_PATH.headerless.bed
    rm $RESULTS_PATH.normalized.bed
    echo "Done normalizing! Continuing filling segments..."
# if the arm coordinates file is not prefixed, but the SEG file does have chromosome names with chr
elif [[ ! "$ARM_IS_PREFIXED" == *"chr"* && "$SEG_IS_PREFIXED" == *"chr"* ]]; then
    echo "Normalizing chr prefix in SEG file to match that of provided arm file..."
    # strip chr prefix
    cat $RESULTS_PATH.headerless.bed | perl -lane '@a=split;$a[0] =~ s/chr//;print join "\t", @a;' > $RESULTS_PATH.normalized.bed
    # replace the headerless fiel with this normalized one
    cat $RESULTS_PATH.normalized.bed > $RESULTS_PATH.headerless.bed
    rm $RESULTS_PATH.normalized.bed
    echo "Done normalizing! Continuing filling segments..."
fi


# Make variable available to use within perl
export THIS_SAMPLE_ID

# Run bedtools substract to find regions that are missing from seg file.
# The perl pipe part will assign the neutral segment log.ratio and no LOH for the missing segments.
# For the missing segments, +1 is added to start position and 1 is substracted from the end position
# to ensure they do not overlap with original segments of seg file by 1 position.
# 0.00 are explicitly used to easily identify segments supplemented by this script
bedtools subtract -a $ARM_BED_PATH -b $RESULTS_PATH.headerless.bed | perl -lane '@a=split;$a[1] = ++$a[1];$a[2] = --$a[2]; $a[3]="0.00"; $a[4]="0.00"; $a[5]=$ENV{THIS_SAMPLE_ID}; print join "\t", @a;' > $RESULTS_PATH.temp

# Merge the initial seg file with the missing segments and rearrange columns to match the style of seg files
cat $RESULTS_PATH.headerless.bed $RESULTS_PATH.temp | sort -k1,1 -k2,2n -V  > $RESULTS_PATH.merged.seg

# Now, remove blacklisted regions (centromeres and p arm telomeres) from this file.
bedtools subtract -a $RESULTS_PATH.merged.seg -b $BLACKLIST_PATH | perl -pale 'BEGIN { $"="\t"; } $_ = "@F[$#F,0..$#F-1]"'  > $RESULTS_PATH.deblacklisted.seg

# Return back the header
cat $RESULTS_PATH.header $RESULTS_PATH.deblacklisted.seg > $RESULTS_PATH

# Cleanup
rm $RESULTS_PATH.temp
rm $RESULTS_PATH.header
rm $RESULTS_PATH.headerless.bed
rm $RESULTS_PATH.merged.seg
rm $RESULTS_PATH.deblacklisted.seg
