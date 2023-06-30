#!/usr/bin/env Rscript

### Generate inputs for Significantly Mutated Rregions (SMR) tools ###
# Uses file defining sample set to generate a seg file
# and prepare it to be used with SMR tools e.g. gistic2

# Usage:
#   Rscript generate_smr_inputs.R <path/to/master/seg> <path/to/sample_sets> <path/to/output/folder> <case_set> <mode>
#   Example: generate_smr_inputs.R results/gambl/gistic2-1.0/00-inputs/capture--projection/all-grch37.seg results/gambl/gistic2-1.0/00-inputs/genome--projection/all-grch37.seg data/metadata/level3_samples_subsets.tsv "FLs_with_LSARP_Trios" gistic2
#
# Notes:
#   Adapted from generate_smg_inputs/1.0/generate_smg_inputs.R
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
  library(data.table)
  library(tidyverse)
})
)

# Determine arguments -----------------------------------------------------

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE) %>% as.list()
arg_names <- c("master_seg", "all_sample_sets", "output_path", "case_set", "mode")
# if there are multiple seg files passed, collapse them into one list
args = c(
        list(unlist(args[1:(length(args)-4)])),
        args[(length(args)-3):length(args)]
        )
args <- setNames(args, arg_names[1:length(args)])

# Print args for debugging
print(paste("master_seg:",args$master_seg,
            "all_sample_sets:",args$all_sample_sets,
            "output_path:",args$output_path,
            "case_set",args$case_set,
            "mode",args$mode))

# Ensure consistent naming of the sample ID column-------------------
if (file.exists(args$all_sample_sets)) {
  full_case_set = suppressMessages(read_tsv(args$all_sample_sets))
} else {
  message(paste("Warning: case set is requested, but the case set file", full_case_set_path, "is not found."))
  stop("Exiting because did not found case set-defining file")
}

full_case_set =
  full_case_set %>% rename_at(vars(matches(
    "sample_id", ignore.case = TRUE
  )),
  ~ "Tumor_Sample_Barcode")

# get case set as defined in the file
this_subset_samples =
  full_case_set %>%
  dplyr::filter(!!sym(args$case_set) == 1) %>%
  pull(Tumor_Sample_Barcode)

# Load master seg files and get regions for the subset-------------------
message("Loading master seg and finding available data for samples in requested subset...")
if (length(args$master_seg)>1){
  message("More than one seg file is supplied. Concatenating them into single file.")
  master_seg =
    tibble(filename = args$master_seg) %>% # create a data frame
    # holding the file names
    mutate(file_contents = map(filename, # read files into
                             ~ read_tsv(args$master_seg, col_types = cols())) # a new data column
    ) %>%
    unnest(cols = c(file_contents)) %>%
    select(-filename)
} else {
  master_seg = suppressWarnings(read_tsv(args$master_seg, col_types = cols()))
}