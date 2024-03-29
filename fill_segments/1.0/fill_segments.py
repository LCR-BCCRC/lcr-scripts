#!/usr/bin/python3

"""
This script will fill empty segments in seg files by creating a dummy segment with log2(CN)=0 and no LOH.
This is an important step preparing seg files for GISTIC after lifting over from one genome build to another.
This script assumes that CN is in log2 ratio.
It requires seg file and chromosome arms file as mandatory inputs. The path to output file also must be specified.

Example:

python3 fill_segments.py --input </path/to/input_file>.seg --chromArm </path/to/chrom_arm_file>.tsv --output </path/to/output_file>.seg
"""

# import required modules
import pandas as pd
import argparse
import simplejson
import sys



def main():
    # initiate the parser and handle arguments from command line
    args = parse_args()
    input_file = args.input
    output_file = args.output
    chrom_file = args.chromArm

    # determine the format of input file
    input_format = input_file[-3:]
    
    # check arguments given in command line
    check_arguments(args, input_format)

    # create a dictionary containing coordinates of chromosome arms
    arm_chrom = load_chrom_arm(chrom_file)
    # get the order of chromosomes
    chrom_order = list(arm_chrom.keys()) + ["buffer"]

    # initialize empty variable for the new segments
    columns_new = []
    columns_edges = []

    # initialize list to store all segments, since it is faster than concatenating pd df with large number of segments
    seg_filled = []

    # assign new vaeriables for ease of referense when creating new segments
    # '0.00' assigned so that it will be easy to identify fill-in segments, and grep if necessary
    empty_loh = str("0.0")
    empty_logr = str("0.00")


    # fill segments
    seg = open(input_file, 'r+') 
    lines=seg.readlines()

    # first, get header of the file
    header=lines[0].rstrip("\n").rstrip("\r").split("\t")

    print("Filling missing segments and smoothing centromeres...")
    # next, go through each segment, skipping the header
    for i in range(1,len(lines)-1):

        # read 2 segments at a time to compare coordinates of end of previous sefment, and start of the next segments
        columns_first = (lines[i].rstrip("\n").rstrip("\r")).split("\t")
        columns_second = (lines[i+1].rstrip("\n").rstrip("\r")).split("\t")
        if columns_first[5]=='-inf': columns_first[5]=str('-10')
        if columns_second[5]=='-inf': columns_second[5]=str('-10')

        # insert empty segment from the beginning of chromosome of the first segment in file to complete the telomeric region of first chromosome
        if i==1: 
            columns_new = [columns_first[0], columns_first[1], str(arm_chrom[columns_first[1]]['p']['start']), str(int(columns_first[2])-1), empty_loh, empty_logr]
            seg_filled.append(columns_new)
            seg_filled.append(columns_first)
            # deal with fencepost problem
            if (int(columns_first[3]) > arm_chrom[columns_first[1]]['p']['end'] and int(columns_first[3]) < arm_chrom[columns_first[1]]['q']['start']):            
                columns_first[3] = str(arm_chrom[columns_first[1]]['p']['end'])
            seg_filled.append(columns_first)

            if (chrom_order[chrom_order.index(columns_second[1])] == chrom_order[chrom_order.index(columns_first[1])+1]):
              missing_arm = chrom_order[chrom_order.index(columns_first[1])]
              columns_edges = [columns_first[0], columns_first[1], str(arm_chrom[missing_arm]['q']['start']), str(arm_chrom[missing_arm]['q']['end']), empty_loh, empty_logr]
              seg_filled.append(columns_edges)
              seg_filled.append(columns_second)
              continue        

        # scenario 1: when it is same sample and same chromosome
        if (columns_first[0]==columns_second[0] and columns_first[1]==columns_second[1]):

            # handle very rare overlapping segments (occurs ~ 0.008%)
            if (int(columns_first[3]) > int(columns_second[2])):
                columns_first[3] = int(columns_second[2])-1
                seg_filled.append(columns_first)
                pass

            # for segments in p arm
            if (int(columns_second[2]) < arm_chrom[columns_second[1]]['p']['end']):
                # create empty segment to fill in
                columns_new = [columns_first[0], columns_first[1], str(int(columns_first[3])+1), str(int(columns_second[2])-1), empty_loh, empty_logr]
                seg_filled.append(columns_new)
                next_segment = (lines[i+1].rstrip("\n").rstrip("\r")).split("\t")
                if (int(columns_second[3]) < arm_chrom[columns_second[1]]['p']['end'] and int(columns_second[3]) < int(next_segment[2])):
                    seg_filled.append(columns_second)
                seg_filled.append(columns_second)    

            # deal with centromeres
            # I already know that this is same sample, and same chromosome
            elif (int(columns_first[2]) < arm_chrom[columns_first[1]]['p']['end'] and int(columns_second[2]) > arm_chrom[columns_second[1]]['p']['end']):

                # first lets deal with end of p arm: segment 1 might end before centromere, or within centromere
                if int(columns_first[3]) < arm_chrom[columns_first[1]]['p']['end']:
                    columns_new = [columns_first[0], columns_first[1], str(int(columns_first[3])+1), str(arm_chrom[columns_first[1]]['p']['end']), empty_loh, empty_logr]
                    seg_filled.append(columns_new)
                # if it extends into centromere, cut segment 1 at the end of p arm
                else:
                    columns_first[3] = str(arm_chrom[columns_first[1]]['p']['end'])
                    seg_filled.append(columns_first)

                # now lets deal with start of q arm: it might start within or after centromere
                if int(columns_second[2]) < arm_chrom[columns_second[1]]['q']['start']:
                    columns_second[2] = str(arm_chrom[columns_second[1]]['q']['start'])
                    seg_filled.append(columns_second)
                    next_segment = (lines[i+1].rstrip("\n").rstrip("\r")).split("\t")
                    if (int(next_segment[2]) < arm_chrom[next_segment[1]]['q']['start'] and int(next_segment[3]) > arm_chrom[next_segment[1]]['q']['start']):
                        next_segment[2] = str(arm_chrom[next_segment[1]]['q']['start'])
                        seg_filled.append(next_segment)

                # possible edge cases around centromere
                else:
                    columns_new = [columns_second[0], columns_second[1], str(arm_chrom[columns_second[1]]['q']['start']), str(int(columns_second[2])-1), empty_loh, empty_logr]
                    previous_segment = (lines[i].rstrip("\n").rstrip("\r")).split("\t")
                    if (previous_segment[0]==columns_second[0] and int(previous_segment[3])>arm_chrom[columns_second[1]]['q']['start']):
                        columns_edges = [columns_second[0], columns_second[1], str(arm_chrom[columns_second[1]]['q']['start']), str(int(previous_segment[3])), columns_first[4], columns_first[5]]
                        seg_filled.append(columns_edges)
                        columns_new = [columns_edges[0], columns_edges[1], str(int(columns_edges[3])+1), str(int(columns_second[2])-1), empty_loh, empty_logr]
                    seg_filled.append(columns_new)
                    seg_filled.append(columns_second)

            # for segments in q arm
            elif (int(columns_first[2]) > arm_chrom[columns_second[1]]['q']['start']):
                # create empty segment to fill in
                columns_new = [columns_first[0], columns_first[1], str(int(columns_first[3])+1), str(int(columns_second[2])-1), empty_loh, empty_logr]
                seg_filled.append(columns_new)
                seg_filled.append(columns_second)

            # some segments are completely within centromere. drop them
            elif (int(columns_first[2]) > arm_chrom[columns_first[1]]['p']['end'] and int(columns_first[2]) < arm_chrom[columns_first[1]]['q']['start']):
                if (int(columns_first[3]) > arm_chrom[columns_first[1]]['q']['start']):
                  columns_new = [columns_first[0], columns_first[1], str(arm_chrom[columns_first[1]]['q']['start']), str(int(columns_first[3])), columns_first[4], columns_first[5]]
                  columns_edges = [columns_first[0], columns_first[1], str(int(columns_first[3])+1), str(int(columns_second[2])-1), empty_loh, empty_logr]
                  seg_filled.append(columns_new)
                else:
                  columns_edges = [columns_second[0], columns_second[1], str(arm_chrom[columns_second[1]]['q']['start']), str(int(columns_second[2])-1), empty_loh, empty_logr]
                  if (int(columns_second[2]) > arm_chrom[columns_second[1]]['p']['end'] and int(columns_second[2]) < arm_chrom[columns_second[1]]['q']['start']):
                    columns_second[2] = str(arm_chrom[columns_second[1]]['q']['start'])
                    seg_filled.append(columns_second)
                seg_filled.append(columns_edges)
                if (int(columns_second[2]) > arm_chrom[columns_second[1]]['p']['end'] and int(columns_second[3]) < arm_chrom[columns_second[1]]['q']['start']):
                  pass # this just drops the segment from output if it is within centromere
                else:
                  if (int(columns_second[2])>arm_chrom[columns_second[1]]['q']['start']):
                    seg_filled.append(columns_second)
            elif (int(columns_second[2]) > arm_chrom[columns_second[1]]['p']['end'] and int(columns_second[2]) < arm_chrom[columns_second[1]]['q']['start']):
                pass # this is handled later

            # did I miss anything? it is possible some edge cases were not considered at time of script development
            else:
              print(columns_first[0], columns_second[0], columns_first[1], columns_second[1], columns_first[2], columns_second[2])
              raise ValueError ("Other sort of way. This is an edge case that needs debugging!")

        # scenario 2: same sample, but going over to the new chromosome
        elif (columns_first[0]==columns_second[0] and columns_first[1]!=columns_second[1]):

            # very rare cases when whole chromosome is missing, identify them here
            if (chrom_order[chrom_order.index(columns_second[1])] != chrom_order[chrom_order.index(columns_first[1])+1]):
              missing_chrom = chrom_order[chrom_order.index(columns_first[1])+1]
              missing_p = [columns_first[0], missing_chrom, str(arm_chrom[missing_chrom]['p']['start']), str(arm_chrom[missing_chrom]['p']['end']), empty_loh, empty_logr]
              missing_q = [columns_first[0], missing_chrom, str(arm_chrom[missing_chrom]['q']['start']), str(arm_chrom[missing_chrom]['q']['end']), empty_loh, empty_logr]
              seg_filled.append(missing_p)
              seg_filled.append(missing_q)
            
            # first, are there any segments in the p arm? that means second segments starts all the way in centromere or q arm
            if (int(columns_first[3]) > arm_chrom[columns_first[1]]['q']['start']):
              if (int(columns_first[2]) > arm_chrom[columns_first[1]]['p']['end'] and int(columns_first[2]) < arm_chrom[columns_first[1]]['q']['start']):
                  previous_segment = (lines[i-1].rstrip("\n").rstrip("\r")).split("\t")
                  if (chrom_order[chrom_order.index(previous_segment[1])] != chrom_order[chrom_order.index(columns_first[1])]):
                    columns_edges = [columns_first[0], columns_first[1], str(arm_chrom[columns_first[1]]['p']['start']), str(arm_chrom[columns_first[1]]['p']['end']), empty_loh, empty_logr]
                    seg_filled.append(columns_edges)                  
                  pass
              elif (int(columns_first[2]) < arm_chrom[columns_first[1]]['p']['end'] and int(columns_first[3]) > arm_chrom[columns_first[1]]['q']['start']):
                  previous_segment = (lines[i-1].rstrip("\n").rstrip("\r")).split("\t")
                  if (chrom_order[chrom_order.index(previous_segment[1])] != chrom_order[chrom_order.index(columns_first[1])-1]):
                      columns_edges = [columns_first[0], columns_first[1], columns_first[2], str(arm_chrom[columns_first[1]]['p']['end']), columns_first[4], columns_first[5]]
                      columns_first[2] = arm_chrom[columns_first[1]]['q']['start']
                      seg_filled.append(columns_edges)
                      seg_filled.append(columns_first)
              else:
                  columns_edges = [columns_first[0], columns_first[1], str(int(columns_first[3])+1), str(arm_chrom[columns_first[1]]['q']['end']), empty_loh, empty_logr]
                  seg_filled.append(columns_edges)
            if (int(columns_second[2]) > arm_chrom[columns_second[1]]['p']['end']):
              if (chrom_order[chrom_order.index(columns_second[1])] != chrom_order[chrom_order.index(columns_first[1])+1]):
                seg_filled.append(missing_p)
                seg_filled.append(missing_q)
              if (int(columns_first[2]) < arm_chrom[columns_first[1]]['q']['start'] and int(columns_first[2]) > arm_chrom[columns_first[1]]['p']['end']):
                columns_first[2] = arm_chrom[columns_first[1]]['q']['start']
                seg_filled.append(columns_first)                  
              columns_new = [columns_second[0], columns_second[1], str(arm_chrom[columns_second[1]]['p']['start']), str(arm_chrom[columns_second[1]]['p']['end']), empty_loh, empty_logr]
              seg_filled.append(columns_new)
              if (int(columns_second[2]) > arm_chrom[columns_second[1]]['q']['start']):
                columns_edges = [columns_second[0], columns_second[1], str(arm_chrom[columns_second[1]]['q']['start']), str(int(columns_second[2])-1), empty_loh, empty_logr]
                seg_filled.append(columns_edges)
              else:
                columns_second[2] = str(arm_chrom[columns_second[1]]['q']['start'])
              seg_filled.append(columns_second)

            # are there any segments in the q arm? that means first segment ends before start of q arm
            elif (int(columns_first[3]) < arm_chrom[columns_first[1]]['q']['start']):
              columns_new = [columns_first[0], columns_first[1], str(arm_chrom[columns_first[1]]['q']['start']), str(arm_chrom[columns_first[1]]['q']['end']), empty_loh, empty_logr]
              seg_filled.append(columns_first)
              if (int(columns_first[3]) < arm_chrom[columns_first[1]]['p']['end']):
                  columns_edges = [columns_first[0], columns_first[1], str(int(columns_first[3])+1), str(arm_chrom[columns_first[1]]['p']['end']), empty_loh, empty_logr]
                  seg_filled.append(columns_edges)
              seg_filled.append(columns_new)
              seg_filled.append(columns_second)              

            # are there any segments that starts in p arm and span centromere? if so, maintain loh flag and logr, but cut out centromere
            elif (int(columns_second[2]) < arm_chrom[columns_second[1]]['p']['end'] and int(columns_second[3]) > arm_chrom[columns_second[1]]['q']['start']):
              previous_segment = (lines[i].rstrip("\n").rstrip("\r")).split("\t")
              if "X" not in str(columns_second[1]):
                next_segment = (lines[i+2].rstrip("\n").rstrip("\r")).split("\t")
                columns_new = [columns_second[0], columns_second[1], str(int(columns_second[2])+1), str(arm_chrom[columns_second[1]]['p']['end']), columns_second[4], columns_second[5]]
                columns_edges = [columns_second[0], columns_second[1], str(arm_chrom[columns_second[1]]['q']['start']), str(int(columns_second[3])), columns_second[4], columns_second[5]]
                if (columns_new[0]==previous_segment[0] and columns_new[1]!=previous_segment[1]):
                    columns_new[2]=str(arm_chrom[columns_new[1]]['p']['start'])
                if (columns_second[0]==next_segment[0] and columns_second[1]!=next_segment[1]):
                    seg_filled.append(columns_new)
                    seg_filled.append(columns_edges)

            # in other cases, there are segments both in p and q arms
            else:
              columns_edges = [columns_second[0], columns_second[1], str(arm_chrom[columns_first[1]]['p']['start']), str(int(columns_second[2])-1), empty_loh, empty_logr]
              if (int(columns_first[2]) > arm_chrom[columns_second[1]]['p']['end']):
                if (int(columns_first[2]) > arm_chrom[columns_first[1]]['p']['end'] and int(columns_first[2]) < arm_chrom[columns_first[1]]['q']['start']):
                  columns_first[2] = arm_chrom[columns_first[1]]['q']['start']
                  seg_filled.append(columns_first)                
                columns_new = [columns_first[0], columns_first[1], str(int(columns_first[3])+1), str(arm_chrom[columns_first[1]]['q']['end']),  empty_loh, empty_logr]
                seg_filled.append(columns_new)
              if (chrom_order[chrom_order.index(columns_second[1])] != chrom_order[chrom_order.index(columns_first[1])+1]):
                seg_filled.append(missing_p)
                seg_filled.append(missing_q)
              seg_filled.append(columns_edges)
              if (int(columns_second[3]) < arm_chrom[columns_second[1]]['p']['end']):
                seg_filled.append(columns_second)


        # scenario 3: new sample, obviously new chromosome
        else:
            previous_segment = (lines[i].rstrip("\n").rstrip("\r")).split("\t")
            columns_edges = [columns_first[0], columns_first[1], str(int(previous_segment[3])+1), str(arm_chrom[columns_first[1]]['q']['end']), empty_loh, empty_logr]
            columns_new = [columns_second[0], columns_second[1], str(arm_chrom[columns_first[1]]['p']['start']), str(int(columns_second[2])-1), empty_loh, empty_logr]
            seg_filled.append(columns_edges)
            seg_filled.append(columns_new)
            seg_filled.append(columns_second)

    seg.close()

    # make df from list of lists and convert chromosome coordinates to integers
    seg_filled_df = pd.DataFrame(seg_filled, columns = header)
    seg_filled_df["start"] = seg_filled_df["start"].astype(int)
    seg_filled_df["end"] = seg_filled_df["end"].astype(int)

    # remove any inverted segments, if there are
    print("Checking and removing inverted segments...")
    seg_filled_df = seg_filled_df[(seg_filled_df["end"]>seg_filled_df["start"])]

    # remove any duplicated segments, if there are
    print("Checking and removing duplicated segments...")
    seg_filled_df = seg_filled_df.groupby((seg_filled_df["end"] != seg_filled_df["end"].shift(-2)).cumsum().values).first()    
    seg_filled_df = seg_filled_df.groupby((seg_filled_df["end"] != seg_filled_df["end"].shift(-1)).cumsum().values).first()
    seg_filled_df = seg_filled_df.groupby((seg_filled_df["start"] != seg_filled_df["start"].shift(-1)).cumsum().values).first()

    # save to the output file specified by user
    print("Saving to file...")
    seg_filled_df.to_csv(output_file, header=True, index=False, sep="\t")
    print("Done!")


# Create nested dictionary to store shromosome arms coordinates. It is adopted from Chris's implementation in other script that summarizes CNVs
def load_chrom_arm(chrom_file):
    arm_chrom = {}
    required_cols = ["chromosome", "start", "end", "arm"]
    header_cols = {}

    i = 0
    with open(chrom_file) as f:
        for line in f:
            i += 1
            line = line.rstrip("\n").rstrip("\r")  # Remove line endings
            cols = line.split("\t")

            # Skip empty lines
            if not line:
                continue

            # If we haven't parsed the header yet, assume this is the first line of the file (aka the header)
            if not header_cols:
                j = 0
                for col in cols:
                    if col in required_cols:
                        header_cols[col] = j
                    j += 1

                # Check to make sure all required columns are found
                for col in required_cols:
                    if col not in header_cols:
                        raise AttributeError("Unable to locate column %s in the chromosome arm positions file \'%s\'" % (col, chrom_file))
                # If we get this far, the header is valid
                continue
            
            if cols[0] not in arm_chrom:
                arm_chrom[cols[0]] = {}
            if cols[3]:
                if cols[3] not in arm_chrom[cols[0]]:
                    arm_chrom[cols[0]][cols[3]]={}
                arm_chrom[cols[0]][cols[3]]['start'] = int(cols[1])
                arm_chrom[cols[0]][cols[3]]['end'] = int(cols[2])
    return arm_chrom


# Check that required arguments are provided, and the input is in .seg format
def check_arguments(args, input_format):
    if input_format == 'seg' and not all([args.input, args.output, args.chromArm]):
        raise ValueError ('Must specify input .seg file, output file, and file listing coordinates of chromosome arms.')
    elif input_format != 'seg':
        raise ValueError ('Input file must be in .seg format')
    else:
      pass


# Parse arguments from command line
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input",
                        help="Imput file in .seg format to fill segments", required=True)
    parser.add_argument("--output",
                        help="Resulting file after filling missing segments", required=True)
    parser.add_argument("--chromArm",
                        help="File with coordinates of chromosme arms for a given genome build", required=True)

    # ignore everything else that is not required by this script
    args, unknown = parser.parse_known_args()
    # return arguments provided by user
    return args


if __name__ == '__main__':
    main()
