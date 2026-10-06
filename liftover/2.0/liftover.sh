#!/usr/bin/env bash

# This script will use UCSC liftOver to convert seg or bed files between hg19/hg38 genome builds.
# The scripts accepts input files with any number of columns.
# The argument 'mode' will determine whether it is `SEG` or 'BED' files to convert.
# The argument 'chain' is expecting the path to liftOver chain file.
# The optional argument 'header' accepts 'YES' or 'NO' to indicate whether header is present in the input file. It defaults to 'NO'
# The optional argument <minmatch> expects a float in the range 0..1 to pass to the liftover minMatch parameter. It defaults to 0.95
# The optional CREATE_LOG enables the tracking of the cumulative segment size before and after liftover (currently only for SEG mode). Defaults to 'NO'

# Usage:
# bash liftover.sh <mode> <input file> <output file> <chain> [header] [minmatch] [create_log]
# Example:
# bash liftover.sh SEG some_hg19_file.seg lifted_hg38_file.seg /some/path/hg19ToHg38.over.chain.gz YES 0.9775 YES
#

# KNOWN LIMITATIONS:
# Lifting large segments commonly fails. In SEG mode, this is dealt with (mostly)
# by breaking each segment into many chunks and merging the chunks afterward. There will still be some loss of segments in this process.
# This chunking process is currently *not* implemented for BED mode.
# Also, this script adds the chr prefix if not present, otherwise liftOver silently outputs empty file
set -eu
IFS=$'\n\t'

# --- helpful helpers ---
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
warn() { printf 'WARN: %s\n' "$*" >&2; }
# Show line on failure (nice for debugging)
trap 'die "Command failed (exit $?) at line $LINENO: ${BASH_COMMAND}"' ERR

export LC_ALL=C
export LANG=C

# Read variables to store the arguments from command line
MODE="$1"
INPUT_FILE="$2"
OUTPUT_FILE="$3"
CHAIN="$4"
HEADER={$5:-"NO"}
MINMATCH="${6:-0.95}"
CREATE_LOG=${7:-"NO"}

# If the input file has header, save it temporarily in a separate file
if [[ "$HEADER" == *"YES"* ]]; then
    head -1 $INPUT_FILE > $OUTPUT_FILE.header
    tail -n+2 $INPUT_FILE > $OUTPUT_FILE.nohead
elif [[ "$HEADER" == *"NO"* ]]; then
    cat $INPUT_FILE > $OUTPUT_FILE.nohead
else
    echo "You specified header $HEADER, which is not recognized. Please specify YES or NO."
    exit 1 # terminate and indicate error
fi

# If SEG mode, rearrange columns and collapse extra columns together
if [[ "$MODE" == *"SEG"* ]]; then
    CHR_COL=2 # this is used near the end to make output chr prefix match that of the input file
    cat "$OUTPUT_FILE.nohead" \
    | perl -ne '
        use strict; use warnings;
        our $S_I; our $E_I;
        BEGIN {
            $S_I=2; $E_I=3;
        }

        chomp; my @a = split /\t/;

        # normalize chr
        my $chr = $a[1];
        $chr = ($chr eq "23") ? "X" : $chr;
        $chr = ($chr eq "24") ? "Y" : $chr;
        $chr = "chr$chr" unless $chr =~ /^chr/;

        my ($s1,$e1) = ($a[$S_I], $a[$E_I]);
        # remove decimal point and trailing digits from coordinates
        # these should theoretically never happen but this protects against it
        $s1 =~ s/\..+//;
        $e1 =~ s/\..+//;

        my $extra = join "|", @a;
        print join("\t", $chr, $s1, $e1, "$extra\n");

    ' > "$OUTPUT_FILE.collapsed"

# Split the segments into chunks
    cat $OUTPUT_FILE.collapsed \
       | perl -ne '
       # split [start,end) into chunks of size <= $chunk_size (BED semantics)
       BEGIN { $chunk_size = 250000;our $segid = 1; }
       chomp; @a = split /\t/;
       my ($chr,$s,$e) = @a[0,1,2];


    # skip empty ranges or where end is less than start
    next if !defined $s || !defined $e || $e <= $s;

    while ($s + $chunk_size < $e) {
        my $ce = $s + $chunk_size;
        print "$chr\t$s\t$ce\t$a[3]|SEGMENT_$segid\n";

        $s = $ce;                # no +1 for BED half-open
        $chunk++;
    }
    print "$chr\t$s\t$e\t$a[3]|SEGMENT_$segid\n";   # remainder, where $ce+250000 > $e, but $ce < $e still
    $segid++;
  ' > "$OUTPUT_FILE.chunked" && rm $OUTPUT_FILE.nohead && rm $OUTPUT_FILE.collapsed
elif [[ "$MODE" == *"BED"* ]]; then
    CHR_COL=1
    cat "$OUTPUT_FILE.nohead" \
    | perl -ne '
        use strict; use warnings;
        our $S_I; our $E_I; our $chr_I;
        BEGIN {
            $S_I=1; $E_I=2; $chr_I = 0;
        }

        chomp; my @a = split /\t/;

        # normalize chr
        my $chr = $a[$chr_I];
        $chr = ($chr eq "23") ? "X" : $chr;
        $chr = ($chr eq "24") ? "Y" : $chr;
        $chr = "chr$chr" unless $chr =~ /^chr/;

        my ($s,$e) = ($a[$S_I], $a[$E_I]);
        # remove decimal point and trailing digits from coordinates
        # these should theoretically never happen but this protects against it
        $s =~ s/\..+//;
        $e =~ s/\..+//;

        my $ncol=3;
        print "$chr\t$s\t$e\t";
        print join "|", @a[$ncol..$#a];
        print "\n";
    ' > "$OUTPUT_FILE.collapsed"

  cat $OUTPUT_FILE.collapsed \
  | perl -ne '
    # split [start,end) into chunks of size <= $chunk_size (BED semantics)
    BEGIN {$| = 1; $chunk_size = 250000; our $segid = 1;}
    chomp; @a = split /\t/;
    my ($chr,$s,$e) = @a[0,1,2];
    my $chunk = 1;

    # skip empty ranges or where end is less than start
    next if !defined $s || !defined $e || $e <= $s;

    while ($s + $chunk_size < $e) {
        my $ce = $s + $chunk_size;
        print "$chr\t$s\t$ce\t$a[3]|SEGMENT_$segid\n";

        $s = $ce;                # no +1 for BED half-open
        $chunk++;
    }
    print "$chr\t$s\t$e\t$a[3]|SEGMENT_$segid\n";   # remainder, where $ce+250000 > $e, but $ce < $e still
    $segid++;
  ' > "$OUTPUT_FILE.chunked" && rm $OUTPUT_FILE.nohead && rm $OUTPUT_FILE.collapsed
else
    echo "You specified mode $MODE, which is not supported. Please provide SEG or BED file."
    exit 1 # terminate and indicate error
fi

UNMAPPED="${OUTPUT_FILE%.*}.unmapped.bed"

liftOver -minMatch=$MINMATCH $OUTPUT_FILE.chunked  $CHAIN $OUTPUT_FILE.chunked.lift.bed $UNMAPPED 2> /dev/null
rm $OUTPUT_FILE.chunked


# Merge adjacent segments with the same collapsed column value, which includes the SEGMENT_# made for each segment above.
sort -V -k1,1 -k2,2n $OUTPUT_FILE.chunked.lift.bed  \
    | perl -s -F'\t' -ane ' # we pass some variables into the script at the end of this code block
    BEGIN{use strict;}
    next if /^\s*$/;                   # skip blank lines

    my $n = @F;                        # number of columns

    my ($chr,$s,$e) = @F[0,1,2];

    my $lastcol = @F[$#F];  # all the concatenated fields from original file plus SEGMENT_#

    if (!defined $cur_id) {
        # seed current run
        ($cur_chr,$cur_s,$cur_e,$cur_id) = ($chr,$s,$e,$lastcol);
        $cur_tail = $lastcol;          # keep column 4
        next;
    }

    if (   $lastcol  eq $cur_id
        && $chr eq $cur_chr
        && $s <= $cur_e + 1            # contiguous (allowing off-by-one stitching)
        ) {
        # extend current run
        $cur_e = $e;

    } else {
        # flush previous run
        print join("\t", $cur_chr, $cur_s, $cur_e, $cur_id);
        # start a new run
        ($cur_id,$cur_chr,$cur_s,$cur_e) = ($lastcol,$chr,$s,$e);
    }

    END {
        if (defined $cur_id) {
        print join("\t", $cur_chr, $cur_s, $cur_e, $cur_id);
        }
    }
    '   | perl -ne 's/\|SEGMENT_\S+$//;s/\|/\t/g;print;' > $OUTPUT_FILE.merged_noheader # removes SEGMENT_# and converts collapsed column back to tab sep

    rm $OUTPUT_FILE.chunked.lift.bed


LOG="$INPUT_FILE.log"
if [[ "$CREATE_LOG" == *"YES"* ]]; then
    # check the total segment size before and after our chopping/lifting and write to a log file
    if [[ "$MODE" == *"SEG"* ]]; then
        cat $INPUT_FILE \
        | perl -ne '
        chomp;
        @a=split("\t");
        next if $a[1] =~ /chrom/;
        $seg_len = $a[3] - $a[2];
        $total+=$seg_len;
        $totals{$a[1]}+=$seg_len;
        END {for(sort keys %totals){print "$_\t$totals{$_}\n";}; print "Cumulative segment length before liftOver:\t$total\n"};
        ' > $LOG
        cat $OUTPUT_FILE \
        | perl -ne 'chomp;
        @a=split("\t");
        next if $a[1] =~ /chrom/;
        $seg_len = $a[3] - $a[2];
        $total+=$seg_len;
        $totals{$a[1]}+=$seg_len;
        END {for(sort keys %totals){print "$_\t$totals{$_}\n";}; print "Cumulative segment length after liftOver:\t$total\n"};
        ' >> $LOG
    fi
fi
# Make output chr prefix match that of the input file
DELIM=$'\t'
HAS_CHR=$(
  awk -v c="$CHR_COL" -v FS="$DELIM" '
    NR==1 { next }                     # skip header
    /^#/ { next }                      # (optional) skip comment lines
    { v = $c; sub(/\r$/, "", v) }      # trim CR if Windows line endings
    v ~ /^chr/ { found=1; exit }       # found any data row starting with "chr"
    END { print found ? 1 : 0 }        # print once
  ' "$INPUT_FILE"
)
echo "chr prefix on column $CHR_COL of input file? $HAS_CHR"

# Next, if the input file had header, merge it back to the lifted file
if [[ "$HEADER" == *"YES"* ]]; then
    if [[ "$MODE" == *"SEG"* ]]; then
      if [[ "$HAS_CHR" == "1" ]]; then
        sort -V -k1,1 -k2,2n $OUTPUT_FILE.merged_noheader \
        | perl -ane 'print "$F[3]\t$F[0]\t$F[1]\t$F[2]\t"  , join("\t", @F[7..$#F]), "\n";' \
        > $OUTPUT_FILE.sort.noheader
      else
        sort -V -k1,1 -k2,2n $OUTPUT_FILE.merged_noheader \
        | perl -ane '$F[0] =~ s/chr//; print "$F[3]\t$F[0]\t$F[1]\t$F[2]\t"  , join("\t", @F[7..$#F]), "\n";' \
        > $OUTPUT_FILE.sort.noheader
      fi
    else
        sort -V -k1,1 -k2,2n $OUTPUT_FILE.merged_noheader \
        > $OUTPUT_FILE.sort.noheader
    fi
    rm $OUTPUT_FILE.merged_noheader
    cat $OUTPUT_FILE.header $OUTPUT_FILE.sort.noheader > $OUTPUT_FILE \
    && rm $OUTPUT_FILE.sort.noheader \
    && rm $OUTPUT_FILE.header # this only exissts if input has header, so clean up this temp file here
else
    if [[ "$MODE" == *"SEG"* ]]; then
        cat $OUTPUT_FILE.merged_noheader \
        | sort -V -k1,1 -k2,2n -V \
        | perl -ane 'print "$F[3]\t$F[0]\t$F[1]\t$F[2]\t"  , join("\t", @F[7..$#F]), "\n";' \
        > $OUTPUT_FILE && rm $OUTPUT_FILE.merged_noheader

    else
        sort -V -k1,1 -k2,2n -V $OUTPUT_FILE.merged_noheader \
        > $OUTPUT_FILE && rm $OUTPUT_FILE.merged_noheader

    fi
fi



echo "Done lifting $INPUT_FILE to $OUTPUT_FILE"
