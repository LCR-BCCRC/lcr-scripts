#!/usr/bin/env Rscript

### Generate inputs for Significantly Mutated Rregions (SMR) tools ###
# Uses file defining sample set to generate a seg file
# and prepare it to be used with SMR tools e.g. gistic2

# Usage:
#   Rscript generate_smg_inputs.R <path/to/master/maf> <path/to/sample_sets> <path/to/output/folder> <case_set> <mode> <include non-coding>
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
