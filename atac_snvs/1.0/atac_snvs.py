#!/usr/bin/env python



##### ATTRIBUTION #####

# MBB 498 Project - Rachel LaFrance



##### ABOUT #####
# This is a Python script intended to identify SNVs from ATAC-Seq BAM files


# FUNCTIONS:

# extract_barcode : will extract CB_Values (barcodes) from the BAM file
# process_reads : requires BAM file and reference genome FASTA file, identifies SNVs
# write_tsv : tsv file
    # snv_file = all SNVs detected
    # snv_multi = multiple barcoded reads per SNV
# maf_comparison : filters MAF file by "SNP", then compares the generated tsv file against it to confirm matches
# find_matching_barcodes : looks for any barcodes that picked up on multiple SNVs for further validation



##### SETUP #####
import pysam
from pyfaidx import Fasta 
from collections import defaultdict
import pandas as pd



### input and output files ###
def main():
    # the destination of the patient BAM file and reference genome to run process_reads function
    bam_file = "" #input file name
    reference_fasta = "reference_genome.fa"
    output = process_reads(bam_file, reference_fasta)

    # name the output files
    snv_file = "mpileup.tsv"
    snv_multi = "multiple_barcodes.tsv"
    write_tsv(snv_file, snv_multi, output)

    # the destination for patient MAF, desired file chosen from write_tsv function above, and named output
    maf_comparison("data/99-13280.maf", snv_file, "MAF_comparison.tsv")

    # use the maf_comparison file above as input
    comparison_file = "MAF_comparison.tsv" 
    output_file = "matched_barcodes.tsv"
    find_matching_barcodes(comparison_file, output_file)



# Extract barcodes from BAM using the CB tag
def extract_barcodes(read):
    try:
        return read.get_tag("CB") 
    except KeyError:
        return None #to bypass KeyError
    


# Identify any SNVs from barcoded reads
def process_reads(bam_file, reference_fasta):
    # open the bam file and reference genome
    bam = pysam.AlignmentFile(bam_file, "rb")
    fasta = Fasta(reference_fasta) 

    # use defaultdict - nested dict to store output data:
    # corresponding to chromosome -> position -> reference_base -> query_base -> set(CB barcodes)
    output = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(set))))

    for read in bam:
        # checks if read maps to reference, and doesn't map to secondary sites and split alignments
        if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
            chrom = bam.get_reference_name(read.reference_id) # chromosome name
            reference_base = fasta[chrom][read.reference_start].seq.upper() # reference base at read start position
            query_base = read.query_sequence.upper()[0]  # query_base : nucleotide from BAM file

            # identify SNVs : remove same reads, and any "N" bases
            if query_base != reference_base and query_base != 'N':
                cb_tag = extract_barcodes(read) # extract the CB barcode 
                output[chrom][read.reference_start][reference_base][query_base].add(cb_tag) # outputs the chr, location, ref and query base, with CB

    bam.close()
    fasta.close()

    return output



# OPTIONS:
# 1 - snv_file : Output tsv file of all SNVs detected from the BAM
# 2 - snv_multi : Output tsv file with more than one barcoded read per SNV, narrowed results 
def write_tsv(snv_file, snv_multi, output):

    with open(snv_file, "w") as output_f:
        output_f.write("Chromosome\tLocation\tReference_Base\tQuery_Base\tNumber_of_Barcodes\tCB_Values\n") # header
    
        with open(snv_multi, "w") as s_output:
            s_output.write("Chromosome\tLocation\tReference_Base\tQuery_Base\tNumber_of_Barcodes\tCB_Values\n")

            # Iterate through the chr -> pos -> ref base _> query in data
            for chrom, positions in output.items():
                for position, data in positions.items():
                    for reference_base, query_bases in data.items():
                        for query_base, barcodes in query_bases.items():
                            num_barcodes = len(barcodes) # counts number of barcodes
                            valid_barcodes = [barcode for barcode in barcodes if barcode is not None]  # to filter out "None" values
                            if valid_barcodes:  # to only proceed with valid barcodes 
                                cb_values = ",".join(valid_barcodes) # makes string
                                output_line = f"{chrom}\t{position}\t{reference_base}\t{query_base}\t{num_barcodes}\t{cb_values}\n"

                                if num_barcodes > 1: # only output more than one barcode
                                    s_output.write(output_line) # snv_multi

                                output_f.write(output_line) # snv_file



# Compares all results from snv_file checked against the MAF file, and adds additional columns from the MAF to the file for more info
def maf_comparison(maf_file, snv_file, output_file):
    # using pandas to read MAF file
    maf_data = pd.read_csv(maf_file, sep='\t', comment='#', header=0, dtype=str)

    # filtering the data based on column 10 and "SNP" values
    filtered_data = maf_data[maf_data['Variant_Type'] == "SNP"].copy()  # avoid SettingWithCopyWarning

    # convert 'Chromosome' to string in order to avoid DtypeWarning message
    filtered_data['Chromosome'] = filtered_data['Chromosome'].astype(str)

    # open snv_output file with pandas
    mpileup_data = pd.read_csv(snv_file, sep='\t', dtype=str)
    columns_to_include = maf_data.columns.tolist() + mpileup_data.columns[2:].tolist()

    # initialize
    merged_data = pd.DataFrame(columns=columns_to_include)

    for _, mpileup_row in mpileup_data.iterrows():
        chromosome = mpileup_row['Chromosome']
        location = float(mpileup_row['Location'])
        location_max = location + 1 # to account for the script - location/position off by 1 compared to MAF

        # search MAF data based on conditions - matching with mpileup
        matching_rows = filtered_data[(filtered_data['Chromosome'] == chromosome) &
                                      (filtered_data['Start_Position'].astype(float) <= location_max) &
                                      (filtered_data['End_Position'].astype(float) >= location) &
                                        (filtered_data['Tumor_Seq_Allele2'] == mpileup_row['Query_Base'])]
    
        if not matching_rows.empty:
            # merge the matching rows with additional columns from mpileup data, ignore index to prevent unnecessary tabs and incorrect formatting
            merged_row = pd.concat([matching_rows.iloc[0], mpileup_row[2:]], axis=0, ignore_index=False)
            merged_data = pd.concat([merged_data, merged_row.to_frame().T], ignore_index=False)

    # merged data to a new tsv file --> "MAF_comparison"
    merged_data.to_csv(output_file, sep='\t', index=False)



# Iterates through CB_Values to identify if multiple SNVs are picked up from the same barcode
def find_matching_barcodes(comparison_file, output_file):
    matching_CB_and_Chr_values = {}
    
    with open(comparison_file, 'r') as comp_f:
        next(comp_f)  # to skip header
        for line in comp_f:
            fields = line.strip().split('\t')
            output_CB_value = fields[48]  # get CB_Values (49th column)
            chromosome = fields[4]  # get Chromsome (5th)

            # Add the row to the list corresponding to its CB_Value and Chromosome
            key = (output_CB_value, chromosome)
            if key not in matching_CB_and_Chr_values:
                matching_CB_and_Chr_values[key] = []
            matching_CB_and_Chr_values[key].append(fields[:13])  # add columns from the MAF for more information

    with open(output_file, 'w') as matching_f:
        matching_f.write("CB_Value\tChromosome\tMAF_Info\n")
        for (CB_value, chromosome), rows in matching_CB_and_Chr_values.items():
            if len(rows) > 1:
                for row in rows:
                    matching_f.write("{}\t{}\t{}\n".format(CB_value, chromosome, '\t'.join(row)))



if __name__ == "__main__":
    main()