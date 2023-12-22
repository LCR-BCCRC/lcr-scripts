# /usr/bin/env Rscript

# Description:
#   Adapted from generate_smg_inputs/1.0/generate_smg_inputs.R.
#   This script is intended for generating input to SMR modules in LCR-modules (gistic2).
#   It expects to be run as part of a snakemake workflow which provides a file with categories
#   to subset the metadata in order to get samples IDs of interest. The snakemake workflow will also provide
#    launch_date, projection, and output directory values.
#   As of right now this script creates input for gistic2 using genome and capture seg files from cnv_master,
#   including resolving overlapping regions.
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

# Determine arguments from snakemake -----------------------------------------------------
subsetting_categories_file <- snakemake@input[["subsetting_categories"]]
full_subsetting_categories <- suppressMessages(read_tsv(subsetting_categories_file, comment="#"))
output_dir <- dirname(snakemake@output[[1]])

sample_set <- snakemake@wildcards[["sample_set"]]
projection <- snakemake@wildcards[["projection"]]
launch_date <- snakemake@wildcards[["launch_date"]]

meta <- snakemake@params[['metadata']]
meta_cols <- snakemake@params[['metadata_cols']]

cat("Arguments from snakemake...\n")
cat(paste("Sample sets file:", subsetting_categories_file, "\n"))
cat(paste("Output directory:", output_dir, "\n"))
cat(paste("Sample set:", sample_set, "\n"))
cat(paste("Launch date:", launch_date, "\n"))

# pandas df from snakemake is passed as a character vector
# This converts it into a dataframe
num_rows <- length(meta)/length(meta_cols)
meta_matrix <- t(matrix(meta, nrow = 25, ncol = num_rows))
# Convert to dataframe and name columns
metadata <- as.data.frame(meta_matrix)
colnames(metadata) <- meta_cols
# Format NA
metadata <- metadata %>%
  mutate_all(~na_if(., ''))

# Get subsetting values for the sample_set
# Renaming the variable required to subset the df correctly
case_set <- sample_set
subsetting_values <- full_subsetting_categories %>%
  filter(sample_set == case_set)

cat("Subsetting values:\n")
print(subsetting_values)

# Function for getting the sample ids
subset_samples <- function(categories, meta) {

  if ("time_point" %in% names(categories)){
    if(categories$time_point == "primary_only"){
      # Make a vector of acceptable values to store in subsetting_values list
      categories$time_point <- "NA,A,1"
    } else if (categories$time_point == "non-primary-only"){
      # Make a vector mutually exclusive with the one above
      categories$time_point <- paste(unique(eta$time_point[!meta$time_point %in% c(NA, "A", "1")]), collapse=",")
    } else if(categories$time_point == "all"){
      categories$time_point <- paste(unique(meta$time_point), collapse=",")
    }
  }

  for (col in colnames(categories)[-1]){
      subset_values <- unlist(strsplit(categories[[col]], ","))
      subset_values[subset_values == "NA"] <- NA
      meta <- meta %>%
          filter(.data[[col]] %in% subset_values)
  }

  samples <- meta %>%
    pull(sample_id)

  return(samples)
}

# Get sample ids of the sample_set
this_sample_set <- subset_samples(subsetting_values, metadata)

# Get seg file paths depending on seq type
seg_files <- snakemake@input[["seg"]]
if ("genome" %in% subsetting_values$seq_type && !("capture" %in% subsetting_values$seq_type)) { # genome only
  cat("Loading genome seg...\n")
  full_seg <- suppressMessages(read_tsv(seg_files[str_detect(seg_files, "genome")])) %>%
    filter(ID %in% this_sample_set)

} else if (!("genome" %in% subsetting_values$seq_type) && "capture" %in% subsetting_values$seq_type) { # capture only
  cat("Loading capture seg...\n")
  full_seg<- suppressMessages(read_tsv(seg_files[str_detect(seg_files, "capture")])) %>%
    filter(ID %in% this_sample_set)

} else if ("genome" %in% subsetting_values$seq_type && "capture" %in% subsetting_values$seq_type) { # both
  cat("Loading genome seg...\n")
  genome_seg <- suppressMessages(read_tsv(seg_files[str_detect(seg_files, "genome")])) %>%
    filter(ID %in% this_sample_set)

  cat("Loading capture seg...\n")
  capture_seg <- suppressMessages(read_tsv(seg_files[str_detect(seg_files, "capture")])) %>%
    filter(!ID %in% unique(genome_seg$ID)) %>%
    filter(ID %in% this_sample_set)

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

# Check if output dir extists, create if not
if (!dir.exists(file.path(output_dir))){
  cat("Output directory for sample_set and launch date combo does not exist. Creating it...\n")
  cat(output_dir,"\n")
  dir.create(file.path(output_dir), recursive = TRUE)
} else {
  cat("Output directory for sample_set and launch date combo exists.\n")
  cat(output_dir,"\n")
}

# Report missing samples and calculate the md5sum-------------------
missing_samples <- setdiff(this_sample_set,
                          unique(full_seg$ID))

if (length(missing_samples)==0) {
  cat(paste("Found regions for all samples. ", length(this_sample_set), "samples will be used in the resulting seg file. \n"))
  final_sample_set <- this_sample_set
  md5sum <- digest(final_sample_set)
} else {
  cat(paste("WARNING: ", length(missing_samples), " samples will not be available for the analysis. \n"))
  cat("Writing missing sample ids to file... \n")
  write_tsv(data.frame(missing_samples), paste0(output_dir, "/", md5sum, "_missing_sample_ids.txt"))
  final_sample_set <- unique(full_seg$ID)
  md5sum <- digest(final_sample_set)
}

# Write out final seg file -------------------
cat("Writing combined seg data to file... \n")
write_tsv(full_seg, paste0(output_dir, "/", md5sum, ".seg"))

# Write out md5sum file -------------------
cat("Writing sample ids to file... \n")
write_tsv(data.frame(final_sample_set), paste0(output_dir, "/", md5sum, "_sample_ids.txt"))

# Writing empty file for snakemake checkpoint rule output
file.create(paste0(output_dir, "/done"))

cat("DONE!")
sink()