#!/usr/bin/env Rscript

# Usage:
#   Rscript generate_smr_inputs.R --genome <path/to/genome/master/seg> --capture <path/to/genome/master/seg> --output_path <path/to/output/folder>
#            --all_sample_sets <path/to/sample_sets> --case_set <case_set>
#
#   Example: Rscript generate_smr_inputs.R --genome /projects/rmorin/projects/gambl-repos/gambl-sgillis/results/gambl/gistic2-1.0/00-inputs/genome--projection/all--grch37.seg --capture /projects/rmorin/projects/gambl-repos/gambl-sgillis/results/gambl/gistic2-1.0/00-inputs/capture--projection/all--grch37.seg 
#             --output_dir /projects/rmorin_scratch/sgillis_temp/lcr-scripts/generate_smr_inputs/1.0/ --all_sample_sets /projects/rmorin/projects/gambl-repos/gambl-sgillis/data/metadata/level3_samples_subsets.tsv --case_set FLs_with_LSARP_Trios
#
# Notes:
#   Adapted from generate_smg_inputs/1.0/generate_smg_inputs.R.
#   This script is intended for generating input to SMR modules in LCR-modules (gistic2).
#   It expects to be provided with the tab-deliminated file where sample subsets for a particular analysis are specified, where the first column (sample_id)
#   defines the unique sample ID, and each column indicates whether this ID is included (1) or not (0) in a particular case set. 
#   The column name for the case set will be used as the naming of the output seg file at the user-provided location.
#   As of right now it creates input for gistic2 using genome and capture seg files from cnv_master
#   It can be expanded to include other seq_types and to format inputs for other SMR tools.

# Load packages -----------------------------------------------------------
cat("Loading packages...")
suppressWarnings(
suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(tidyverse)
})
)

# Parse command-line arguments -----------------------------------------------------
parser <- ArgumentParser(description="Generates inputs for Significantly Mutated Regions (SMR) tools. 
Uses file defining sample sets and a case set name to generate a seg file from genome and capture seg files. 
Currently prepares input for gistic2. Genome data takes precedence over capture data.")

parser$add_argument("--genome", "-g", nargs=1, type= 'character', help="Path to the genome--projection/all--{projection}.seg file")
parser$add_argument("--capture", "-c", nargs=1, type= 'character', help="Path to the capture--projection/all--{projection}.seg file")
parser$add_argument("--projection", "-p", nargs=1, type= 'character', required=TRUE, help="Genome build projection. e.g. hg38 or grch37")
parser$add_argument("--output_dir", "-o", nargs=1, type= 'character', required=TRUE, help="Path to write the combined seg file for the case set")
parser$add_argument("--all_sample_sets", nargs=1, type= 'character', required=TRUE, help="Tab delimited file where the first column is sample ID 
                                                        and the rest of the columns are named after case sets. 
                                                        Samples will have a 1 in a column if they are part of that case set, 
                                                        and zero otherwise.")
parser$add_argument("--case_set", nargs=1, required=TRUE, help="Name of the case set to subset the region data to. Must match a column name in ALL_SAMPLE_SETS")

# Gets the args as a named list
args <- parser$parse_args()

# Check that at least one of -g or -c is given
if ( (length(args$genome)==0) && (length(args$capture)==0) ){
  stop(paste("Exiting because at least one genome or capture seg file was not provided"))
}

# Check existance of sample set file -----------------------------------------------------
if (file.exists(args$all_sample_sets)) {
  full_case_set = suppressMessages(read_tsv(args$all_sample_sets))
} else {
  stop(paste("Exiting because sample sets file", args$all_sample_sets, "is not found."))
}

full_case_set =
  full_case_set %>% rename_at(vars(matches(
    "sample_id", ignore.case = TRUE
  )),
  ~ "ID")

# Get sample IDs of the case_set
case_set_samples =
  full_case_set %>%
  dplyr::filter(!!sym(args$case_set) == 1) %>%
  pull(ID)

# Load genome seg file and get regions for the caseset -------------------
if (length(args$genome) != 0){ # if -g value providced
  cat("Loading genome seg...")
  if (!file.exists(args$genome)) {
    stop(paste("Exiting because genome data seg file", args$genome, "is not found."))
  } else {
    cat("Finding available data for samples in requested case set...")
    genome_seg <- suppressMessages(read_tsv(args$genome, col_types = cols())) %>%
      filter(ID %in% case_set_samples)
  }
}

# Load capture seg file and get regions for the caseset -------------------
if (length(args$capture) != 0){ # if -c value provided
  cat("Loading capture seg...")
  if (!file.exists(args$capture)) {
    stop(paste("Exiting because capture data seg file", args$capture, "is not found."))
  } else if (length(args$genome) != 0){ # if -g provided
    cat("Removing samples already in genome data and finding available data for samples in requested case set...")
    capture_seg <- suppressMessages(read_tsv(args$capture, col_types = cols())) %>%
      filter(!ID %in% unique(genome_seg$ID)) %>%
      filter(ID %in% case_set_samples)
  } else { # if -g not provided
    cat("Finding available data for samples in requested case set...")
    capture_seg <- suppressMessages(read_tsv(args$capture, col_types = cols())) %>%
      filter(ID %in% case_set_samples)
  }
}

# Merge genome and capture data where applicable -------------------
if ( (length(args$genome)!=0) && (length(args$capture)!=0) ){ # if both -g and -c have values provided
  cat("Merging genome and capture data ...")
  full_seg <- rbind(genome_seg, capture_seg)
} else if ( length(args$genome)!=0 && (length(args$capture)==0) ){ # only -g provided
  full_seg <- genome_seg
} else { # only -c provided, since it would have exited earlier if both weren't given
  full_seg <- capture_seg
}

# Report missing samples -------------------
missing_samples <- setdiff(case_set_samples,
                          unique(full_seg$ID))

if (length(missing_samples)==0) {
  cat(paste("Found regions for all samples. ", length(case_set_samples), "samples will be used in the resulting seg file."))
} else {
  cat(paste("WARNING: ", length(missing_samples), " samples will not be available for the analysis."))
  cat("Did not find regions for these samples in the combine seg data:")
  cat(missing_samples)
}

# Write out final seg file -------------------
cat("Writing combined seg data to file...")
write_tsv(full_seg, paste0(args$output_dir, "/", args$case_set, "--", args$projection, ".seg"))

cat("DONE!")