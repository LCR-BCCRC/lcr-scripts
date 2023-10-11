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
cat("Loading packages... \n")
suppressWarnings(
suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(tidyverse)
  library(stringr)
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
  stop(paste("Exiting because at least one genome or capture seg file was not provided. \n"))
}

# Check existance of sample set file -----------------------------------------------------
if (file.exists(args$all_sample_sets)) {
  full_case_set = suppressMessages(read_tsv(args$all_sample_sets))
} else {
  stop(paste("Exiting because sample sets file", args$all_sample_sets, "is not found. \n"))
}

full_case_set =
  full_case_set %>% rename_at(vars(matches(
    "sample_id", ignore.case = TRUE
  )),
  ~ "ID")

# Get sample IDs of the case_set
case_set_samples <-full_case_set %>%
  dplyr::filter(!!sym(args$case_set) == 1) %>%
  pull(ID)

# Load genome seg file and get regions for the caseset -------------------
if (length(args$genome) != 0){ # if -g value providced
  cat("Loading genome seg... \n")
  if (!file.exists(args$genome)) {
    stop(paste("Exiting because genome data seg file", args$genome, "is not found. \n"))
  } else {
    cat("Finding available data for samples in requested case set... \n")
    genome_seg <- suppressMessages(read_tsv(args$genome, col_types = cols())) %>%
      filter(ID %in% case_set_samples)
  }
}

# Load capture seg file and get regions for the caseset -------------------
if (length(args$capture) != 0){ # if -c value provided
  cat("Loading capture seg... \n")
  if (!file.exists(args$capture)) {
    stop(paste("Exiting because capture data seg file", args$capture, "is not found. \n"))
  } else if (length(args$genome) != 0){ # if -g provided
    cat("Removing samples already in genome data and finding available data for samples in requested case set... \n")
    capture_seg <- suppressMessages(read_tsv(args$capture, col_types = cols())) %>%
      filter(!ID %in% unique(genome_seg$ID)) %>%
      filter(ID %in% case_set_samples)
  } else { # if -g not provided
    cat("Finding available data for samples in requested case set... \n")
    capture_seg <- suppressMessages(read_tsv(args$capture, col_types = cols())) %>%
      filter(ID %in% case_set_samples)
  }
}

# Merge genome and capture data where applicable -------------------
if ( (length(args$genome)!=0) && (length(args$capture)!=0) ){ # if both -g and -c have values provided
  cat("Merging genome and capture data ... \n")
  full_seg <- rbind(genome_seg, capture_seg)
} else if ( length(args$genome)!=0 && (length(args$capture)==0) ){ # only -g provided
  full_seg <- genome_seg
} else { # only -c provided, since it would have exited earlier if both weren't given
  full_seg <- capture_seg
}

# Sort by chrom, start, end
full_seg <- full_seg %>%
  arrange(ID, chrom, start, end)

# Filter to only canonical chromosomes -------------------
  cat("Filtering to only canonical chromosomes... \n")
if (args$projection %in% "hg38"){
  full_seg <- full_seg %>%
    filter(str_detect(chrom, regex("chr[XY\\d]+$", ignore_case = TRUE)))
} else if (args$projection %in% "grch37"){
  full_seg <- full_seg %>%
    filter(str_detect(chrom, regex("^[XY\\d]+$", ignore_case = TRUE)))
}

# Remove possible overlaps -------------------
cat("Resolving overlapping regions... \n")
check_overlap = function(seg) {
    highest_end = 0
    overlap <- c()
    region_sizes <- c()
    for (i in 1:nrow(seg)) {
        if (i>1 && seg$ID[i] == seg$ID[i-1] && seg$chrom[i] == seg$chrom[i-1]) {
            if (seg$start[i] >= highest_end) {
                overlap[i] = "NOToverlap"
                region_sizes[i] = "FALSE"
            }else{
                overlap[i] = "overlap"
                region_sizes[i] = (seg$end[i] - seg$start[i])
            }
            if (seg$end[i] > highest_end) {
                highest_end = seg$end[i]
            }
        } else {
            highest_end = seg$end[i]
            overlap[i] = "NA"
            region_sizes[i] = (seg$end[i] - seg$start[i])
        }
    }
    seg <- seg %>% mutate(overlap_status = overlap, region_size = region_sizes)
    return(seg)
}

solve_overlap = function(seg) {
    num_overlap = which(seg$overlap_status == "overlap")
    num_pre_overlap_sorted = (unique(sort(c(num_overlap-1,num_overlap))))
    non_overlap = seg[-num_pre_overlap_sorted,]
    seg <- seg[num_pre_overlap_sorted,]
    for (i in 1:nrow(seg)) {
        if (seg$overlap_status[i] == "overlap") {
            if (seg$end[i] < seg$end[i-1]){
                new_row1 <- data.frame(ID = seg$ID[i],
                                       chrom = seg$chrom[i],
                                       start = seg$start[i-1],
                                       end = seg$start[i],
                                       LOH_flag = seg$LOH_flag[i-1],
                                       log.ratio = seg$log.ratio[i-1],
                                       overlap_status = "NOToverlap",
                                       region_size = "FALSE")
                new_row2 <- data.frame(ID = seg$ID[i],
                                       chrom = seg$chrom[i],
                                       start = seg$end[i],
                                       end = seg$end[i-1],
                                       LOH_flag = seg$LOH_flag[i-1],
                                       log.ratio = seg$log.ratio[i-1],
                                       overlap_status = "NOToverlap",
                                       region_size = "FALSE")
                seg <- seg[-(i-1), ]
                seg <- rbind(seg, new_row1, new_row2)
            }else if (seg$start[i-1] == seg$start[i]) {
                new_row  <- data.frame(ID = seg$ID[i],
                                       chrom = seg$chrom[i],
                                       start = seg$end[i-1],
                                       end = seg$end[i],
                                       LOH_flag = seg$LOH_flag[i],
                                       log.ratio = seg$log.ratio[i],
                                       overlap_status = "NOToverlap",
                                       region_size = "FALSE")
                seg <- seg[-(i), ]
                seg <- rbind(seg, new_row)
            }else if (seg$region_size[i] < seg$region_size[i-1]) {
                new_row1 <- data.frame(ID = seg$ID[i],
                                       chrom = seg$chrom[i],
                                       start = seg$start[i-1],
                                       end = seg$start[i],
                                       LOH_flag = seg$LOH_flag[i-1],
                                       log.ratio = seg$log.ratio[i-1],
                                       overlap_status = "NOToverlap",
                                       region_size = "FALSE")
                new_row2 <- data.frame(ID = seg$ID[i],
                                       chrom = seg$chrom[i],
                                       start = seg$start[i],
                                       end = seg$end[i-1],
                                       LOH_flag = seg$LOH_flag[i],
                                       log.ratio = seg$log.ratio[i],
                                       overlap_status = "NOToverlap",
                                       region_size = "FALSE")
                new_row3 <- data.frame(ID = seg$ID[i],
                                       chrom = seg$chrom[i],
                                       start = seg$end[i-1],
                                       end = seg$end[i],
                                       LOH_flag = seg$LOH_flag[i],
                                       log.ratio = seg$log.ratio[i],
                                       overlap_status = "NOToverlap",
                                       region_size = "FALSE")
                seg <- seg[-c(i, i-1), ]
                seg <- rbind(seg, new_row1, new_row2, new_row3)
            }else if (seg$region_size[i] > seg$region_size[i-1]) {
                new_row1 <- data.frame(ID = seg$ID[i],
                                       chrom = seg$chrom[i],
                                       start = seg$start[i-1],
                                       end = seg$start[i],
                                       LOH_flag = seg$LOH_flag[i-1],
                                       log.ratio = seg$log.ratio[i-1],
                                       overlap_status = "NOToverlap",
                                       region_size = "FALSE")
                new_row2 <- data.frame(ID = seg$ID[i],
                                       chrom = seg$chrom[i],
                                       start = seg$start[i],
                                       end = seg$end[i-1],
                                       LOH_flag = seg$LOH_flag[i-1],
                                       log.ratio = seg$log.ratio[i-1],
                                       overlap_status = "NOToverlap",
                                       region_size = "FALSE")
                new_row3 <- data.frame(ID = seg$ID[i],
                                       chrom = seg$chrom[i],
                                       start = seg$end[i-1],
                                       end = seg$end[i],
                                       LOH_flag = seg$LOH_flag[i],
                                       log.ratio = seg$log.ratio[i],
                                       overlap_status = "NOToverlap",
                                       region_size = "FALSE")
                seg <- seg[-c(i, i-1), ]
                seg <- rbind(seg, new_row1, new_row2, new_row3)
            }
        }
        seg <- seg %>%
        arrange(ID, chrom, start, end)
    }
    seg = seg %>% arrange(ID, chrom, start, end) 
    seg = check_overlap(seg)
    while("overlap" %in% seg$overlap_status){
        seg = check_overlap(solve_overlap(seg))
    }
    seg = rbind(non_overlap, seg) %>% 
      arrange(ID, chrom, start, end) %>%
      select(-overlap_status, -region_size) %>%
      filter(!start == end)
    return(seg)
}

full_seg_checked <- check_overlap(full_seg)
# if no overlaps, do not run solve function
if ("overlap" %in% full_seg_checked$overlap_status) {
  full_seg <- solve_overlap(full_seg_checked)
} else {
  full_seg <- full_seg_checked %>%
  select(-overlap_status, -region_size)
}

# Report missing samples -------------------
missing_samples <- setdiff(case_set_samples,
                          unique(full_seg$ID))

if (length(missing_samples)==0) {
  cat(paste("Found regions for all samples. ", length(case_set_samples), "samples will be used in the resulting seg file. \n"))
} else {
  cat(paste("WARNING: ", length(missing_samples), " samples will not be available for the analysis. \n"))
  cat("Did not find regions for these samples in the combine seg data: \n")
  cat(missing_samples)
}

# Write out final seg file -------------------
cat("Writing combined seg data to file... \n")
write_tsv(full_seg, paste0(args$output_dir, "/", args$case_set, "--", args$projection, ".seg"))

cat("DONE!")