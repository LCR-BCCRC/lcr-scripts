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

set -Eeuo pipefail
IFS=$'\n\t'

# --- helpful helpers ---
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
warn() { printf 'WARN: %s\n' "$*" >&2; }

require_cmd() {
  for c in "$@"; do
    command -v "$c" >/dev/null 2>&1 || die "Required command not found in PATH: $c"
  done
}

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
#echo "===HEADER: $HEADER MINMATCH: $MINMATCH CREATE_LOG: $CREATE_LOG==="

# Second, if the input file has header, save it temporarily in a separate file
# Prepend with chr if chromosomes are not prepended, othervise liftOver silently outputs empty file

if [[ "$HEADER" == *"YES"* ]]; then
    #echo "Header is specified as $HEADER, handling it separately ..."
    head -1 $INPUT_FILE > $OUTPUT_FILE.header
    tail -n+2 $INPUT_FILE > $OUTPUT_FILE.nohead
elif [[ "$HEADER" == *"NO"* ]]; then
    #echo "Header is specified as $HEADER, the first entry of $MODE file will not be treated as header. Processing ..."
    cat $INPUT_FILE > $OUTPUT_FILE.nohead
else
    echo "You specified header $HEADER, which is not recognized. Please specify YES or NO."
    exit 1 # terminate and indicate error
fi


# First, check that proper mode is specified, rearrange columns for seg file, and collapse extra columns together
if [[ "$MODE" == *"SEG"* ]]; then
    #echo "Running in $MODE mode ..."
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
        $chr = "chr$chr" unless $chr =~ /^chr/;
        

        # convert to 0-based BED coords
        my ($s1,$e1) = ($a[$S_I]-1, $a[$E_I]);  # [s1,e1)
        # remove decimal point and trailing digits from coordinates
        # these should theoretically never happen but this protects against it
        $s1 =~ s/\..+//;
        $e1 =~ s/\..+//;
        
        my $extra = join "|", @a;
        print join("\t", $chr, $s1, $e1, "$extra\n");

    ' > "$OUTPUT_FILE.collapsed"

cat $OUTPUT_FILE.collapsed \
| perl -ne '
    # split [start,end) into chunks of size <= $chunk_size (BED semantics)
    BEGIN { $chunk_size = 150000;our $segid = 1; }
    chomp; @a = split /\t/;
    my ($chr,$s,$e) = @a[0,1,2];


    # guard bad/empty ranges
    next if !defined $s || !defined $e || $e <= $s;

    while ($s + $chunk_size < $e) {
        my $ce = $s + $chunk_size;
        print "$chr\t$s\t$ce\t$a[3]|SEGMENT_$segid\n";

        $s = $ce;                # no +1 for BED half-open
        $chunk++;
    }
    print "$chr\t$s\t$e\t$a[3]|SEGMENT_$segid\n";   # remainder (only if $s < $e)
    $segid++;
  ' > "$OUTPUT_FILE.chunked" && rm $OUTPUT_FILE.nohead && rm $OUTPUT_FILE.collapsed
elif [[ "$MODE" == *"BED"* ]]; then
    #echo "Running in the $MODE mode ..."
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
    echo "wrote $OUTPUT_FILE.collapsed"

cat $OUTPUT_FILE.collapsed \
| perl -ne '
    # split [start,end) into chunks of size <= $chunk_size (BED semantics)
    BEGIN {$| = 1; $chunk_size = 150000; our $seg_id = 1;}
    chomp; @a = split /\t/;
    my ($chr,$s,$e) = @a[0,1,2];
    my $chunk = 1;

    # guard bad/empty ranges
    next if !defined $s || !defined $e || $e <= $s;
    
    while ($s + $chunk_size < $e) {
        my $ce = $s + $chunk_size;
        print "$chr\t$s\t" . ($ce-1) . "\t$a[3]|SEGMENT_$seg_id\n";

        $s = $ce;                # no +1 for BED half-open
    }
    
    print "$chr\t$s\t$e\t$a[3]|SEGMENT_$seg_id\n";   # remainder (only if $s < $e)
    $seg_id++;
  ' > "$OUTPUT_FILE.chunked" && rm $OUTPUT_FILE.nohead && rm $OUTPUT_FILE.collapsed
else
    echo "You specified mode $MODE, which is not supported. Please provide SEG or BED file."
    exit 1 # terminate and indicate error
fi

UNMAPPED="${OUTPUT_FILE%.*}.unmapped.bed"

liftOver -minMatch=$MINMATCH $OUTPUT_FILE.chunked  $CHAIN $OUTPUT_FILE.chunked.lift.bed $UNMAPPED 2> /dev/null
rm $OUTPUT_FILE.chunked

# Now, split back all concatenated columns into the separate ones and rearrange back if it is SEG file
# Also merge all segments that are adjacent and share the same values 

sort -k1,1 -k2,2n $OUTPUT_FILE.chunked.lift.bed  \
    | perl -s -F'\t' -ane ' # we pass some variables into the script at the end of this code block
    BEGIN{use strict;}
    next if /^\s*$/;                     # skip blank lines

    # If there is a header (first field literally "ID"), print & skip
    if ($. == 1 && $F[0] eq "ID") { print join("\t", @F); next; }
    
    my $n = @F;                          # number of columns
    
    my ($chr,$s,$e) = @F[0,1,2];

    
    my $lastcol = @F[$#F];  #all the concatenated fields from original file

    if (!defined $cur_id) {
        # seed current run
        ($cur_chr,$cur_s,$cur_e,$cur_id) = ($chr,$s,$e,$lastcol);
        $cur_tail = $lastcol;         # keep column 4
        next;
    }

    if (   $lastcol  eq $cur_id
        && $chr eq $cur_chr
        && $s <= $cur_e + 1              # contiguous (allowing off-by-one stitching)
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
    '   | perl -ne 's/\|SEGMENT_\S+$//;s/\|/\t/g;print;' > $OUTPUT_FILE.merged_noheader

    rm $OUTPUT_FILE.chunked.lift.bed

        # merge bed and ensure all columns beyond 3 are identical after removing CHUNK_ID


LOG="$INPUT_FILE.log"
if [[ "$CREATE_LOG" == *"YES"* ]]; then
    #check the total segment size before and after our chopping/lifting and write to a log file
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

# Next, if the input file had header, merge it back to the lifted file
if [[ "$HEADER" == *"YES"* ]]; then
    if [[ "$MODE" == *"SEG"* ]]; then
        sort -k1,1 -k2,2n $OUTPUT_FILE.merged_noheader \
        | perl -ane 'print "$F[3]\t$F[0]\t$F[1]\t$F[2]\t"  , join("\t", @F[$#F-1..$#F]), "\n";' \
        > $OUTPUT_FILE.sort.noheader 
    else
        sort -k1,1 -k2,2n $OUTPUT_FILE.merged_noheader \
        > $OUTPUT_FILE.sort.noheader 
    fi
    rm $OUTPUT_FILE.merged_noheader
    cat $OUTPUT_FILE.header $OUTPUT_FILE.sort.noheader > $OUTPUT_FILE \
    && rm $OUTPUT_FILE.sort.noheader \
    && rm $OUTPUT_FILE.header # this is only specific if input has header, so clean up this temp file here
else
    if [[ "$MODE" == *"SEG"* ]]; then
        cat $OUTPUT_FILE.merged_noheader \
        | sort -k1,1 -k2,2n -V \
        | perl -ane 'print "$F[3]\t$F[0]\t$F[1]\t$F[2]\t"  , join("\t", @F[$#F-1..$#F]), "\n";' \
        > $OUTPUT_FILE && rm $OUTPUT_FILE.merged_noheader
        #|
        #| perl -ne 's/\|/\t/g;print;' \
    else
        sort -k1,1 -k2,2n -V $OUTPUT_FILE.merged_noheader \
        > $OUTPUT_FILE && rm $OUTPUT_FILE.merged_noheader

    fi
fi



echo "Done lifting $INPUT_FILE to $OUTPUT_FILE"
