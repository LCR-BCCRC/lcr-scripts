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

set -e

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
        ## Apply blacklisted_hg38.bed ##
        use strict; use warnings;
        our $S_I; our $E_I; 
        BEGIN {
            $S_I=2; $E_I=3; 
        }

        chomp; my @a = split /\t/;
        # pass header (should no longer be needed)
        #if (/seg\.mean|log\.ratio/) { print join("\t",$a[1],$a[2],$a[3], join("|",@a)),"\n"; next; }

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
    BEGIN { $chunk_size = 150000 }
    chomp; @a = split /\t/;
    my ($chr,$s,$e) = @a[0,1,2];
    my $chunk = 1;

    # guard bad/empty ranges
    next if !defined $s || !defined $e || $e <= $s;

    while ($s + $chunk_size < $e) {
        my $ce = $s + $chunk_size;
        print "$chr\t$s\t$ce\t$a[3]|CHUNK_$chunk\n";

        $s = $ce;                # no +1 for BED half-open
        $chunk++;
    }
    print "$chr\t$s\t$e\t$a[3]|CHUNK_$chunk\n";   # remainder (only if $s < $e)

  ' > "$OUTPUT_FILE.chunked" && rm $OUTPUT_FILE.nohead
elif [[ "$MODE" == *"BED"* ]]; then
    #echo "Running in the $MODE mode ..."
    cat $OUTPUT_FILE.nohead \
    | perl -ne '@a=split("\t"); $a[0] = "chr$a[0]" unless $a[0] =~ /^chr/; $ncol=3;while($ncol>0){$first = shift @a;print "$first\t";$ncol--;}print join "|", @a;' \
    > $OUTPUT_FILE.collapsed && rm $OUTPUT_FILE.nohead
else
    echo "You specified mode $MODE, which is not supported. Please provide SEG or BED file."
    exit 1 # terminate and indicate error
fi

UNMAPPED="${OUTPUT_FILE%.*}.unmapped"

# Running liftOver
if [[ "$MODE" == *"SEG"* ]]; then
    #run liftover on chunked seg data
    #echo "liftOver -minMatch=$MINMATCH $OUTPUT_FILE.chunked  $CHAIN $OUTPUT_FILE.chunked.lift.bed $UNMAPPED"
    liftOver -minMatch=$MINMATCH $OUTPUT_FILE.chunked  $CHAIN $OUTPUT_FILE.chunked.lift.bed $UNMAPPED 2> /dev/null

    cat $OUTPUT_FILE.chunked.lift.bed \
        | perl -ne '@a=split;@b = split /\|/, $a[3];
                   print "$b[0]\t$a[0]\t$a[1]\t$a[2]\t" . join("\t",@b[$#b-2..$#b-1]) . "\n";' \
        > $OUTPUT_FILE.lifted-temp.bed 
    rm $OUTPUT_FILE.chunked.lift.bed
    rm $OUTPUT_FILE.collapsed
    rm $OUTPUT_FILE.chunked
    #echo "DONE running liftover, created $OUTPUT_FILE.lifted-temp.bed"
else
    # run liftOver
    #echo "Running UCSC liftOver with minimum match of converted regions set to $MINMATCH ..."
    liftOver -minMatch=$MINMATCH $OUTPUT_FILE.collapsed $CHAIN $OUTPUT_FILE.lifted-temp.bed $UNMAPPED && rm $OUTPUT_FILE.collapsed 2> /dev/null

fi


# Now, split back all concatenated columns into the separate ones and rearrange back if it is SEG file
# Also merge all segments that are adjacent and share the same values 
if [[ "$MODE" == *"SEG"* ]]; then
    #echo "merging adjacent segments from chunking"
    # Ensure sorted by ID, chrom, then start
    cp $OUTPUT_FILE.lifted-temp.bed $OUTPUT_FILE.unmerged
    sort -k1,1 -k2,2V -k3,3n $OUTPUT_FILE.unmerged | \
    perl -F'\t' -ane '
    BEGIN { $EPS = 1e-6; }               # tolerance for float equality (seg.mean)
    next if /^\s*$/;                     # skip blank lines

    # If there is a header (first field literally "ID"), print & skip
    if ($. == 1 && $F[0] eq "ID") { print join("\t", @F); next; }

    my $n = @F;                          # number of columns
    my ($id,$chr,$s,$e) = @F[0,1,2,3];
    my ($last2,$last1) = @F[$n-2,$n-1];  # last two columns (e.g. seg.mean, C)

    # Float-equality helper (for the penultimate column, typically seg.mean)
    sub same_float { my ($a,$b,$eps)=@_; return abs($a-$b) <= $eps; }

    if (!defined $cur_id) {
        # seed current run
        ($cur_id,$cur_chr,$cur_s,$cur_e) = ($id,$chr,$s,$e);
        @cur_tail = @F[4..$#F];            # keep columns 5..end from the first chunk
        ($cur_last2,$cur_last1) = ($last2,$last1);
        next;
    }

    # Check if we should merge with current run
    if (   $id  eq $cur_id
        && $chr eq $cur_chr
        && $last1 eq $cur_last1
        && ( $last2 eq $cur_last2 || same_float($last2,$cur_last2,$EPS) )
        && $s <= $cur_e + 1              # contiguous (allowing off-by-one stitching)
        ) {
        # extend current run
        $cur_e = $e;
    } else {
        # flush previous run
        print join("\t", $cur_id,$cur_chr,$cur_s,$cur_e,@cur_tail);
        # start a new run
        ($cur_id,$cur_chr,$cur_s,$cur_e) = ($id,$chr,$s,$e);
        @cur_tail = @F[4..$#F];
        ($cur_last2,$cur_last1) = ($last2,$last1);
    }

    END {
        if (defined $cur_id) {
        print join("\t", $cur_id,$cur_chr,$cur_s,$cur_e,@cur_tail);
        }
    }
    ' > $OUTPUT_FILE.merged_noheader
    rm $OUTPUT_FILE.lifted-temp.bed
    rm $OUTPUT_FILE.unmerged
    #rm $UNMAPPED
fi

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
    
    sort -k1,1 -k2,2n $OUTPUT_FILE.merged_noheader > $OUTPUT_FILE.sort.noheader \
    && rm $OUTPUT_FILE.merged_noheader
    cat $OUTPUT_FILE.header $OUTPUT_FILE.sort.noheader > $OUTPUT_FILE \
    && rm $OUTPUT_FILE.sort.noheader \
    && rm $OUTPUT_FILE.header # this is only specific if input has header, so clean up this temp file here
else
    if [[ "$MODE" == *"SEG"* ]]; then
        cat $OUTPUT_FILE.merged_noheader \
        | sort -k1,1 -k2,2n -V \
        | perl -ne 's/\|/\t/g;print;' \
        > $OUTPUT_FILE && rm $OUTPUT_FILE.merged_noheader
    else
        sort -k1,1 -k2,2n -V $OUTPUT_FILE.lifted-temp.bed \
        > $OUTPUT_FILE && rm $OUTPUT_FILE.lifted-temp.bed

    fi
fi



echo "Done lifting $INPUT_FILE to $OUTPUT_FILE"
