#!/usr/bin/env Rscript

### Generate inputs for Significantly Mutated Rregions (SMR) tools ###
# Uses file defining sample set to generate a seg file
# and prepare it to be used with SMR tools e.g. gistic2

# Usage:
#   Rscript generate_smr_inputs.R --genome <path/to/genome/master/seg> --capture <path/to/genome/master/seg> --output_path <path/to/output/folder>
#            --all_sample_sets <path/to/sample_sets> --case_set <case_set>
#
#   Example: Rscript generate_smr_inputs.R --genome /projects/rmorin/projects/gambl-repos/gambl-sgillis/results/gambl/gistic2-1.0/00-inputs/genome--projection/all--grch37.seg --capture /projects/rmorin/projects/gambl-repos/gambl-sgillis/results/gambl/gistic2-1.0/00-inputs/capture--projection/all--grch37.seg 
#             --output_dir /projects/rmorin_scratch/sgillis_temp/lcr-scripts/generate_smr_inputs/1.0/ --all_sample_sets /projects/rmorin/projects/gambl-repos/gambl-sgillis/data/metadata/level3_samples_subsets.tsv --case_set FLs_with_LSARP_Trios
#
# Notes:
#   Adapted from generate_smg_inputs/1.0/generate_smg_inputs.R.
#   This script is intended for use with the SMR modules in LCR-modules (gistic2).
#   It expects to be provided with the tab-deliminated file where sample subsets
#   for a particular analysis are specified, where first column (sample_id/Tumor_Sample_Barcode)
#   defines the unique sample ID, and each column indicates whether this ID is included (1) or not (0)
#   in a particular subset. The column name for the subset will be used as the naming of the
#   output seg file at the user-provided location.
#

# Load packages -----------------------------------------------------------
message("Loading packages...")
suppressWarnings(
suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(tidyverse)
})
)

# Parse command-line arguments -----------------------------------------------------
# TO DO : fill in script description in ArgumentParser()
parser <- ArgumentParser()
parser$add_argument("--genome", "-g", nargs=1, help="Path to the genome--projection/all--{projection}.seg file")
parser$add_argument("--capture", "-c", nargs=1, help="Path to the capture--projection/all--{projection}.seg file")
parser$add_argument("--output_dir", "-o", nargs=1, help="Path to write the combined seg file for the case set")
parser$add_argument("--all_sample_sets", nargs=1, help="Tab delimited file where the first column is sample ID 
                                                        and the rest of the columns are named after case sets. 
                                                        Samples will have a 1 in a column if they are part of that case set, 
                                                        and zero otherwise.")
parser$add_argument("--case_set", nargs=1, help="Name of the case set to subset the region data to. Must match a column name in all_sample_sets file")

# Gets the args as a named list
args <- parser$parse_args()

# Print args for debugging
print(paste("genome:",args$genome,
            "capture:",args$capture,
            "output_dir:",args$output_dir,
            "all_sample_sets:",args$all_sample_sets,
            "case_set:",args$case_set))


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
#if (args$case_set){
  case_set_samples =
    full_case_set %>%
    dplyr::filter(!!sym(args$case_set) == 1) %>%
    pull(ID)
# } else {
#   stop(paste("Case_set is not specified."))
# }

# Load genome seg file and get regions for the  caseset-------------------
message("Loading genome seg and finding available data for samples in requested case set...")
if (!file.exists(args$genome)) {
  stop(paste("Exiting because genome data seg file", args$genome, "is not found."))
} else {
  genome_seg <- suppressMessages(read_tsv(args$genome, col_types = cols())) %>%
    filter(ID %in% case_set_samples)
}

# Load capture seg file -------------------
message("Loading capture seg and removing samples already in genome data...")
if (!file.exists(args$capture)) {
  stop(paste("Exiting because genome data seg file", args$capture, "is not found."))
} else {
  capture_seg <- suppressMessages(read_tsv(args$capture, col_types = cols())) %>%
    filter(!ID %in% unique(genome_seg$ID)) %>%
    filter(ID %in% case_set_samples)
}

# Merge genome and capture data -------------------
message("Merging genome and capture data...")
full_seg <- rbind(genome_seg, capture_seg)

# Report missing samples -------------------
missing_samples <- setdiff(case_set_samples,
                          unique(full_seg$ID))

if (length(missing_samples)==0) {
  message(paste("Found regions for all samples. ", length(case_set_samples), "samples will be used in the resulting seg file."))
} else {
  message(paste("WARNING: ", length(missing_samples), " samples will not be available for the analysis."))
  message("Did not find regions for these samples in the combine seg data:")
  message(cat(missing_samples))
}

# Write out final seg file -------------------
message("Writing combined seg data to file...")
write_tsv(full_seg, paste0(args$output_dir, "/", args$case_set, ".seg"))

message("DONE!")