# /usr/bin/env Rscript

# Usage:
#   Rscript generate_smr_inputs.R --genome <path/to/genome/master/seg> --capture <path/to/genome/master/seg> --output_path <path/to/output/folder>
#            --all_sample_sets <path/to/sample_sets> --case_set <case_set>
#
#   Example: Rscript generate_smr_inputs.R --genome /projects/rmorin/projects/gambl-repos/gambl-sgillis/results/gambl/gistic2-1.0/00-inputs/genome--projection/all--grch37.seg --capture /projects/rmorin/projects/gambl-repos/gambl-sgillis/results/gambl/gistic2-1.0/00-inputs/capture--projection/all--grch37.seg
#             --output_dir /projects/rmorin_scratch/sgillis_temp/lcr-scripts/generate_smr_inputs/1.0/ --all_sample_sets /projects/rmorin/projects/gambl-repos/gambl-sgillis/data/metadata/level3_subsetting_categories.tsv --case_set FLs_with_LSARP_Trios
#
# TO DO: update this description
# Notes:
#   Adapted from generate_smg_inputs/1.0/generate_smg_inputs.R.
#   This script is intended for generating input to SMR modules in LCR-modules (gistic2).
#   It expects to be provided with the tab-deliminated file where sample subsets for a particular analysis are specified, where the first column (sample_id)
#   defines the unique sample ID, and each column indicates whether this ID is included (1) or not (0) in a particular case set.
#   The column name for the case set will be used as the naming of the output seg file at the user-provided location.
#   As of right now it creates input for gistic2 using genome and capture seg files from cnv_master
#   It can be expanded to include other seq_types and to format inputs for other SMR tools.

# Log both the stdout and stderr
log <- file(snakemake@log[[1]], open="wt")
sink(log ,type = "output")
sink(log, type = "message")

# Load packages -----------------------------------------------------------
cat("Loading packages... \n")
suppressWarnings(
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(stringr)
  library(digest)
})
)

# Input snakemake variables
subsetting_categories_file <- snakemake@input[["subsetting_categories"]]

case_set <- snakemake@wildcards[["case_set"]]
projection <- snakemake@wildcards[["projection"]]
launch_date <- snakemake@wildcards[["launch_date"]]

output_dir <- snakemake@config[["lcr-modules"]][["gistic2"]][["dirs"]][["prepare_seg"]]

seq_type <- unlist(snakemake@params[["seq_type"]])
metadata_str <- snakemake@params[["metadata"]]

metadata <- data.frame(sample_id=metadata_str[c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)],
          seq_type=metadata_str[c(FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE)],
          genome_build=metadata_str[c(FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE)],
          cohort=metadata_str[c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE)],
          pathology=metadata_str[c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE)],
          unix_group=metadata_str[c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE)],
          time_point=metadata_str[c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)])
# Format NA
metadata <- metadata %>%  
  mutate_all(~na_if(., ''))

full_subsetting_categories <- suppressMessages(read_tsv(subsetting_categories_file))

# Get subsetting values for the case_set
subsetting_values <- full_subsetting_categories %>%
  filter(sample_set == case_set)

cat("Subsetting values:\n")
print(subsetting_values)

# Function for getting the sample ids
subset_samples <- function(categories, metadata) {
    samples <- metadata %>%
            select(sample_id,  seq_type, genome_build,
            cohort, pathology, time_point, unix_group) %>%
            filter(seq_type %in% unlist(strsplit(categories$seq_type, ",")),
            genome_build %in% unlist(strsplit(categories$genome_build, ",")),
            cohort %in% unlist(strsplit(categories$cohort, ",")),
            pathology %in% unlist(strsplit(categories$pathology, ",")),
            unix_group %in% unlist(strsplit(categories$unix_group, ",")),
            if (is.na(time_point) || categories$time_points == "primary-only") time_point %in% c(NA,"A") 
            else if (is.na(time_point) || categories$time_points == "all") time_point %in% c(NA,"A","B","C","D","E","G","F","J","H") 
            else time_point %in% c("B","C","D","E","G","F","J","H")
            ) %>% 
            pull(sample_id)

    return(samples)
}

# Get sample ids of the case_set
case_set_samples <- subset_samples(subsetting_values, metadata)

# Get seg file paths depending on seq type
seg_files <- snakemake@input[["seg"]]
if ("genome" %in% seq_type && !("capture" %in% seq_type)) { # genome only
  cat("Loading genome seg...\n")
  full_seg <- suppressMessages(read_tsv(seg_files[str_detect(seg_files, "genome")])) %>%
    filter(ID %in% case_set_samples)

} else if (!("genome" %in% seq_type) && "capture" %in% seq_type) { # capture only
  cat("Loading capture seg...\n")
  full_seg<- suppressMessages(read_tsv(seg_files[str_detect(seg_files, "capture")])) %>%
    filter(ID %in% case_set_samples)

} else if ("genome" %in% seq_type && "capture" %in% seq_type) { # both
  cat("Loading genome seg...\n")
  genome_seg <- suppressMessages(read_tsv(seg_files[str_detect(seg_files, "genome")])) %>%
    filter(ID %in% case_set_samples)

  cat("Loading capture seg...\n")
  capture_seg <- suppressMessages(read_tsv(seg_files[str_detect(seg_files, "capture")])) %>%
    filter(!ID %in% unique(genome_seg$ID)) %>%
    filter(ID %in% case_set_samples)

  full_seg <- rbind(genome_seg, capture_seg)
}

# Sort by chrom, start, end
full_seg <- full_seg %>%
  arrange(ID, chrom, start, end)

# Filter to only canonical chromosomes -------------------
cat("Filtering to only canonical chromosomes... \n")
if (projection %in% "hg38"){
  full_seg <- full_seg %>%
    filter(str_detect(chrom, regex("chr[XY\\d]+$", ignore_case = TRUE)))
} else if (projection %in% "grch37"){
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

# Report missing samples and calculate the md5sum-------------------
missing_samples <- setdiff(case_set_samples,
                          unique(full_seg$ID))

if (length(missing_samples)==0) {
  cat(paste("Found regions for all samples. ", length(case_set_samples), "samples will be used in the resulting seg file. \n"))
  final_sample_set <- case_set_samples
  md5sum <- digest(final_sample_set)
} else {
  cat(paste("WARNING: ", length(missing_samples), " samples will not be available for the analysis. \n"))
  cat("Did not find regions for these samples in the combine seg data: \n")
  print(missing_samples)
  final_sample_set <- full_seg %>% unique(full_seg$ID)
  md5sum <- digest(final_sample_set)
}

full_output_dir <- paste0(output_dir, case_set, "--", projection, "--", launch_date)

# Check if output dir extists, create if not
if (!dir.exists(file.path(full_output_dir))){
  cat("Output directory for case_set and launch date combo does not exist. Creating it...\n")
  cat(full_output_dir,"\n")
  dir.create(file.path(full_output_dir), recursive = TRUE)
} else {
  cat("Output directory for case_set and launch date combo exists.\n")
  cat(full_output_dir,"\n")
}

# Write out final seg file -------------------
cat("Writing combined seg data to file... \n")
write_tsv(full_seg, paste0(full_output_dir, "/", md5sum, ".seg"))

# Write out md5sum file -------------------
cat("Writing sample ids to file... \n")
write_tsv(data.frame(final_sample_set), paste0(full_output_dir, "/", md5sum, "_sample_ids.txt"))

# Writing empty file for snakemake checkpoint rule output
file.create(paste0(full_output_dir, "/done"))

cat("DONE!")
sink()