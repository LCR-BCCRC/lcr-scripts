#!/usr/bin/env bash

# This script will use bedtools subtract to leverage chromosome coordinates and seg file and fill all missing segments.
# Important! The chromosome prefix will be normalize that of the input arm bed file.
# This script is intended to be used within snakemake and expects to have environment with bedtools available.
# The parameter 'mode' is passed as argument on the command line and will determine whether it is `SEG` files to fill
# or `subclones` output filess from CNV caller battenberg.

# Usage:
# bash fill_segments.sh <input chromosome arm bed file> <input seg file> <input blacklisted bed file> <output seg file> <sample_id> <mode>
# Eexample:
# bash fill_segments.sh src/chromArm.hg19.bed TCRBOA7-T-WEX-test--matched.igv.seg src/blacklisted.hg19.bed TCRBOA7-T-WEX-test--matched.igv.filled.seg TCRBOA7-T-WEX SEG
#

set -e

# Read variables to store the arguments from command line
ARM_BED_PATH="$1"
INPUT_PATH="$2"
BLACKLIST_PATH="$3"
RESULTS_PATH="$4"
THIS_SAMPLE_ID="$5"
MODE="$6"


# bedtools subtract doesn't like headers - so we need to strip it and store separately in temp file
cat $INPUT_PATH | head -1 > $RESULTS_PATH.header

if [[ "$MODE" == *"SEG"* ]]; then
    # Rearrange columns in the seg file to follow bed convention of chr-start-end and strip header
    cat $INPUT_PATH | grep -v "start" | perl -pale 'BEGIN { $"="\t"; } $_ = "@F[1..$#F,0]"' > $RESULTS_PATH.headerless.bed
elif [[ "$MODE" == *"subclones"* ]] || [[ "$MODE" == *"sequenza"* ]] ; then
    # Strip header
    cat $INPUT_PATH | grep -v "start" > $RESULTS_PATH.headerless.bed
else
    echo "You specified mode $MODE, which is not supported. Please provide seg file or battenberg subclones."
    exit 1 # terminate and indicate error
fi

# Count number of columns  - to be used to fill the subclones file
NUM_COLUMNS=$(cat $RESULTS_PATH.headerless.bed | perl -lane 'print $#F; exit')

# If there is inconsistent chr prefixing between SEG file and bed file of chromosome coordinates,
# normalize it to match between the two files.
ARM_IS_PREFIXED=$(head -1 $ARM_BED_PATH | cut -f 1)
INPUT_IS_PREFIXED=$(head -1 $RESULTS_PATH.headerless.bed | cut -f 1)

# if the arm coordinates file is prefixed, but the SEG file is not
if [[ "$ARM_IS_PREFIXED" == *"chr"* && ! "$INPUT_IS_PREFIXED" == *"chr"* ]]; then
    echo "Normalizing chr prefix in SEG file to match that of provided arm file..."
    # add chr prefix
    cat $RESULTS_PATH.headerless.bed | perl -lane '@a=split;$a[0] =~ s/^/chr/;print join "\t", @a;' > $RESULTS_PATH.normalized.bed
    # replace the headerless fiel with this normalized one
    cat $RESULTS_PATH.normalized.bed > $RESULTS_PATH.headerless.bed
    rm $RESULTS_PATH.normalized.bed
    echo "Done normalizing! Continuing filling segments..."
# if the arm coordinates file is not prefixed, but the SEG file does have chromosome names with chr
elif [[ ! "$ARM_IS_PREFIXED" == *"chr"* && "$INPUT_IS_PREFIXED" == *"chr"* ]]; then
    echo "Normalizing chr prefix in SEG file to match that of provided arm file..."
    # strip chr prefix
    cat $RESULTS_PATH.headerless.bed | perl -lane '@a=split;$a[0] =~ s/chr//;print join "\t", @a;' > $RESULTS_PATH.normalized.bed
    # replace the headerless fiel with this normalized one
    cat $RESULTS_PATH.normalized.bed > $RESULTS_PATH.headerless.bed
    rm $RESULTS_PATH.normalized.bed
    echo "Done normalizing! Continuing filling segments..."
fi


# Run bedtools substract to find regions that are missing from seg file.
# The perl pipe part will assign the neutral segment log.ratio and no LOH for the missing segments.
# For the missing segments, +1 is added to start position and 1 is substracted from the end position
# to ensure they do not overlap with original segments of seg file by 1 position.
# Filled segments are given normal copy number values and are tracked by the dummy_segment column.

if [[ "$MODE" == *"SEG"* ]]; then
    # Make variable available to use within perl
    export THIS_SAMPLE_ID
    bedtools subtract -a $ARM_BED_PATH -b $RESULTS_PATH.headerless.bed | perl -lane '@a=split;$a[1] = ++$a[1];$a[2] = --$a[2]; $a[3]="99"; $a[4]="99"; $a[5]=$ENV{THIS_SAMPLE_ID}; print join "\t", @a;' > $RESULTS_PATH.temp

elif [[ "$MODE" == *"subclones"* ]]; then
    # Make variable available to use within perl
    export NUM_COLUMNS
    bedtools subtract -a $ARM_BED_PATH -b $RESULTS_PATH.headerless.bed \
    | perl -lane '@a=split;$a[1] = ++$a[1];$a[2] = --$a[2]; $a[3]="0.5"; $a[4]="1.0"; $a[5]="0.00"; $a[6]="2.0"; foreach my $column (7..8) {$a[$column]="1.0"}; $a[9]="99"; foreach my $column (10..$ENV{NUM_COLUMNS}) {$a[$column]="NA"}; print join "\t", @a;' > $RESULTS_PATH.temp

elif [[ "$MODE" == *"sequenza"* ]]; then
    # Make variable available to use within perl
    export NUM_COLUMNS
    bedtools subtract -a $ARM_BED_PATH -b $RESULTS_PATH.headerless.bed \
    | perl -lane '@a=split;$a[1] = ++$a[1];$a[2] = --$a[2]; foreach my $column (3..5) {$a[$column]="NA"}; $a[6]="1"; foreach my $column (7..8) {$a[$column]="NA"}; $a[9]="99"; foreach my $column (10..11) {$a[$column]="1"}; $a[$ENV{NUM_COLUMNS}]="NA"; print join "\t", @a;' > $RESULTS_PATH.temp

fi


# Merge the initial seg file with the missing segments and rearrange columns to match the style of seg files
cat $RESULTS_PATH.headerless.bed $RESULTS_PATH.temp | sort -k1,1 -k2,2n -V | perl -lane 'print if ($F[2]-$F[1])>1;' > $RESULTS_PATH.merged.seg

# Now, remove blacklisted regions (centromeres and p arm telomeres) from this file.

if [[ "$MODE" == *"SEG"* ]]; then
    # After substracting blacklisted regions, seg file needs the column with sample ID to be shifted back to first position
    bedtools subtract -a $RESULTS_PATH.merged.seg -b $BLACKLIST_PATH | perl -pale 'BEGIN { $"="\t"; } $_ = "@F[$#F,0..$#F-1]"' | perl -lane 'print if ($F[3]-$F[2])>1;' > $RESULTS_PATH.deblacklisted.seg

    # Return back the header with dummy_segment column added
    cat $RESULTS_PATH.header | perl -p -e 'chomp; $_ .= "\tdummy_segment\n"' > $RESULTS_PATH

    # Populate dummy_segment column indicating whether the segment is filled in (1) or is an
    # original output (0), then change the "99" to "0.0"
    cat $RESULTS_PATH.deblacklisted.seg | perl -F"\t" -lane 'print join("\t", @F, (($F[4] == 99 && $F[5] == 99) ? 1 : 0))' | perl -lane '@a=split;$a[4]="0.0"; $a[5]="0.0"; print join "\t", @a;' >> $RESULTS_PATH
elif [[ "$MODE" == *"sequenza"* ]] ; then
    # Sequenza txt file already has first 3 columns as bed-like and does not need column rearrangement
    bedtools subtract -a $RESULTS_PATH.merged.seg -b $BLACKLIST_PATH | perl -lane 'print if ($F[2]-$F[1])>1;' > $RESULTS_PATH.deblacklisted.seg

    # Return back the header with dummy_segment column added
    cat $RESULTS_PATH.header | perl -p -e 'chomp; $_ .= "\tdummy_segment\n"' > $RESULTS_PATH

    # Populate dummy_segment column indicating whether the segment is filled in (1) or is an
    # original output (0), then change the "99" to "2"
    cat $RESULTS_PATH.deblacklisted.seg | perl -F"\t" -lane 'print join("\t", @F, (($F[9] == 99) ? 1 : 0))' | perl -lane '@a=split;$a[9]="2"; print join "\t", @a;' >> $RESULTS_PATH
elif [[ "$MODE" == *"subclones"* ]]; then
    # Subclones file already has first 3 columns as bed-like and does not need column rearrangement
    bedtools subtract -a $RESULTS_PATH.merged.seg -b $BLACKLIST_PATH | perl -lane 'print if ($F[2]-$F[1])>1;' > $RESULTS_PATH.deblacklisted.seg

    # Return back the header with dummy_segment column added
    cat $RESULTS_PATH.header | perl -p -e 'chomp; $_ .= "\tdummy_segment\n"' > $RESULTS_PATH

    # Populate dummy_segment column indicating whether the segment is filled in (1) or is an
    # original output (0), then change the "99" to "1.0"
    cat $RESULTS_PATH.deblacklisted.seg | perl -F"\t" -lane 'print join("\t", @F, (($F[9] == 99) ? 1 : 0))' | perl -lane '@a=split;$a[9]="1.0"; print join "\t", @a;' >> $RESULTS_PATH
fi


# Cleanup
rm $RESULTS_PATH.temp
rm $RESULTS_PATH.header
rm $RESULTS_PATH.headerless.bed
rm $RESULTS_PATH.merged.seg
rm $RESULTS_PATH.deblacklisted.seg
