#!/usr/bin/env Rscript

# Description:
#   Adapted from generate_smg_inputs/1.0/generate_smg_inputs.R.
#   This script is intended for use with the SMG (MutSig2CV, dNdS, HotMAPS,
#   OncodriveFML) and SMR (gistic2) modules in in LCR-modules.
#   It expects to be run as part of a snakemake workflow which provides a file with categories
#   to subset the metadata in order to get samples IDs of interest. The snakemake workflow will also provide
#   seq_type, launch_date, and output directory values.
#   The expected input file type for SMR modules is seg and for SMG modules it is maf.

# Log both the stdout and stderr
log <- file(snakemake@log[[1]], open="wt")
sink(log ,type = "output")
sink(log, type = "message")

# Load packages -----------------------------------------------------------
message("Loading packages...")
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
launch_date <- snakemake@wildcards[["launch_date"]]

mode <- snakemake@params[["mode"]]
if ("maf" %in% names(snakemake@input)){
  include_non_coding <- snakemake@params[["include_non_coding"]]
} else if ("seg" %in% names(snakemake@input)){
  projection <- snakemake@wildcards[["projection"]]
}

meta <- snakemake@params[['metadata']]
meta_cols <- snakemake@params[['metadata_cols']]

cat("Arguments from snakemake...\n")
cat(paste("Sample sets file:", subsetting_categories_file, "\n"))
cat(paste("Output directory:", output_dir, "\n"))
cat(paste("Sample set:", sample_set, "\n"))
cat(paste("Launch date:", launch_date, "\n"))
cat(paste("Mode:", mode, "\n"))
if ("maf" %in% names(snakemake@input)){
  cat(paste("Include non-coding:", include_non_coding, "\n"))
}
if ("seg" %in% names(snakemake@input)){
  cat(paste("Projection:", projection, "\n"))
}

# Determine sample ids in sample set -----------------------------------------------------
# pandas df from snakemake is passed as a character vector
# This converts it into a dataframe
num_rows <- length(meta)/length(meta_cols)
meta_matrix <- t(matrix(meta, nrow = 25, ncol = num_rows))
# Format NA
meta_matrix[meta_matrix==""] <- NA
# Convert to dataframe and name columns
metadata <- as.data.frame(meta_matrix)
colnames(metadata) <- meta_cols

# Get subsetting values for this sample_set
# Renaming the variable required to subset the df correctly
case_set <- sample_set
subsetting_values <- full_subsetting_categories %>%
  filter(sample_set == case_set) %>%
  as.list(.)

# Split comma sep values
subsetting_values <- lapply(subsetting_values, function(x) unlist(str_split(x, ",")))
# Change "NA" character to NA object
subsetting_values <- lapply(subsetting_values,function(x) ifelse(x=="NA",NA,x))

cat("Subsetting values (list):\n")
print(subsetting_values)

# Function for getting the sample ids
subset_samples <- function(categories, meta) {

  if ("time_point" %in% names(categories) && length(categories$time_point)==1){
    if(categories$time_point == "primary-only"){
      # Make a vector of acceptable values to store in subsetting_values list
      categories$time_point <- c(NA, "A", "1")
    } else if (categories$time_point == "non-primary-only"){
      # Make a vector mutually exclusive with the one above
      categories$time_point <- unique(meta$time_point[!meta$time_point %in% c(NA, "A", "1")])
    } else if(categories$time_point == "all"){
      categories$time_point <- unique(meta$time_point)
    }
  }

  for (col in names(categories)[-1]){
    if(length(categories[[col]]) == 1 & is.na(categories[[col]])) { # excludes the NA case in time_point
      next
    } else { 
      meta <- meta %>%
          filter(.data[[col]] %in% categories[[col]])
    }
  }

  samples <- meta %>%
    pull(sample_id)

  return(samples)
}

# Get sample ids of the sample_set
this_subset_samples <- subset_samples(subsetting_values, metadata)

# Load master input files and get variants for the sample set-------------------
if ("maf" %in% names(snakemake@input)){
  input_files <- snakemake@input[["maf"]]
} else if ("seg" %in% names(snakemake@input)) {
  input_files <- snakemake@input[["seg"]]
}

if ("genome" %in% subsetting_values$seq_type && !("capture" %in% subsetting_values$seq_type)) { # genome only
  cat("Loading genome input file ...\n")
  genome_input <- suppressMessages(read_tsv(input_files[str_detect(input_files, "genome")]))
  cat("Subsetting to sample set...\n")

  if ("maf" %in% names(snakemake@input)){
    subset_input <- genome_input %>%
      filter(Tumor_Sample_Barcode %in% this_subset_samples)

  } else if ("seg" %in% names(snakemake@input)) {
    subset_input <- genome_input %>%
      filter(ID %in% this_sample_set)
  }
} else if (!("genome" %in% subsetting_values$seq_type) && "capture" %in% subsetting_values$seq_type) { # capture only
  cat("Loading capture input file...\n")
  capture_input <- suppressMessages(read_tsv(input_files[str_detect(input_files, "capture")]))

  cat("Subsetting to sample set...\n")
  if ("maf" %in% names(snakemake@input)){
    subset_input <- capture_input %>%
      filter(Tumor_Sample_Barcode %in% this_subset_samples)

  } else if ("seg" %in% names(snakemake@input)) {
    subset_input <- capture_input %>%
      filter(ID %in% this_sample_set)
  }
} else if ("genome" %in% subsetting_values$seq_type && "capture" %in% subsetting_values$seq_type) { # both
  cat("Loading genome input file ...\n")
  genome_input <- suppressMessages(read_tsv(input_files[str_detect(input_files, "genome")]))

  cat("Loading capture input file...\n")
  capture_input <- suppressMessages(read_tsv(input_files[str_detect(input_files, "capture")]))

  cat("Subsetting to sample set...\n")
  if ("maf" %in% names(snakemake@input)){
    genome_subset <- genome_input %>%
      filter(Tumor_Sample_Barcode %in% this_subset_samples)

    capture_subset <- capture_input %>%
      filter(!Tumor_Sample_Barcode %in% unique(genome_subset$Tumor_Sample_Barcode)) %>%
      filter(Tumor_Sample_Barcode %in% this_subset_samples)

    subset_input <- rbind(genome_subset, capture_subset)

  } else if ("seg" %in% names(snakemake@input)) {
    genome_subset <- genome_input %>%
      filter(ID %in% this_sample_set)

    capture_subset <- capture_input %>%
      filter(!ID %in% unique(genome_subset$ID)) %>%
      filter(ID %in% this_subset_samples)

    subset_input <- rbind(genome_subset, capture_subset)
  }
}

# Check if output dir extists, create if not-------------------
if (!dir.exists(file.path(output_dir))){
  cat("Output directory for sample_set and launch date combo does not exist. Creating it...\n")
  cat(output_dir,"\n")
  dir.create(file.path(output_dir), recursive = TRUE)
} else {
  cat("Output directory for sample_set and launch date combo exists.\n")
  cat(output_dir,"\n")
}

# Report missing samples and calculate the md5sum-------------------
if ("maf" %in% names(snakemake@input)){
  missing_samples = setdiff(this_subset_samples,
                          unique(subset_input$Tumor_Sample_Barcode))
} else if ("seg" %in% names(snakemake@input)){
  missing_samples = setdiff(this_subset_samples,
                          unique(subset_input$ID))
}
if (length(missing_samples)==0) {
  cat(paste("Success! Found variants for all samples.", length(this_subset_samples), "samples will be used in the analysis\n"))
  md5sum <- digest(this_subset_samples)
  final_sample_set <- this_subset_samples
} else {
  cat(paste("WARNING:", length(missing_samples), "will not be available for the analysis.\n"))
  cat("Writing missing sample ids to file... \n")
  if ("maf" %in% names(snakemake@input)){
    final_sample_set <- unique(subset_input$Tumor_Sample_Barcode)
  } else if ("seg" %in% names(snakemake@input)){
    final_sample_set <- unique(subset_input$ID)
  }
  md5sum <- digest(final_sample_set)
  write_tsv(data.frame(missing_samples), paste0(output_dir, "/", md5sum, "_missing_sample_ids.txt"))
}

# Subset for coding only if user requested -------------------
if ("maf" %in% names(snakemake@input)) {
  if (include_non_coding) {
    cat("Proceeding with non-coding mutations...\n")
  } else{
    cat("Excluding non-coding mutations...\n")
    coding_class = c("Frame_Shift_Del",
                  "Frame_Shift_Ins",
                  "In_Frame_Del",
                  "In_Frame_Ins",
                  "Missense_Mutation",
                  "Nonsense_Mutation",
                  "Nonstop_Mutation",
                  "Silent",
                  "Splice_Region",
                  "Splice_Site",
                  "Targeted_Region",
                  "Translation_Start_Site")
      subset_input <-subset_input %>%
        dplyr::filter(Variant_Classification %in% coding_class)
  }
}

# Format input according to the requirements of each individual tool --------------------------------------
cat("Preparing input data to be used with", mode, "\n")

if (mode == "MutSig2CV") {
  # check that the maf file is in grch37-based coordinates
  if (grepl("38", subset_input$NCBI_Build[1])) {
    cat("Requested mode is MutSig2CV, but the supplied file is in the hg38-based coordinates.\n")
    cat("Unfortunatelly, MutSig only works for grch37-based maf files.\n")
    stop("Please supply the mutation data in grch37-based version.")
  }

  subset_input <- subset_input %>%
    rename("chr"="Chromosome",
           "pos"="Start_Position",
           "gene"="Hugo_Symbol",
           "patient"="Tumor_Sample_Barcode",
           "ref_allele"="Reference_Allele",
           "newbase"="Tumor_Seq_Allele2",
           "type"="Variant_Classification",
           "classification"="Variant_Type")

  subset_input <- subset_input %>%
    select(chr, pos, gene, patient, ref_allele, newbase, type, classification)

  grouping_column <- "patient"
}

if (mode == "dNdS") {
  # check that the maf file is in grch37-based coordinates
  if (grepl("38", subset_input$NCBI_Build[1])) {
    cat("Requested mode is dNdS, but the supplied file is in the hg38-based coordinates.\n")
    cat("Unfortunatelly, dNdS is configured to only work for grch37-based maf files.\n")
    stop("Please supply the mutation data in grch37-based version.")
  }

  subset_input <- subset_input %>%
      select(Tumor_Sample_Barcode,
              Chromosome,
              Start_Position,
              Reference_Allele,
              Tumor_Seq_Allele2) %>%
      `names<-`(c("sampleID",
                  "chr",
                  "pos",
                  "ref",
                  "mut"))

  grouping_column <- "sampleID"
}

if (mode == "FishHook") {
  subset_input <- subset_input %>%
    select(Hugo_Symbol,
           Tumor_Sample_Barcode,
           Chromosome,
           Start_Position,
           End_Position,
           Variant_Classification,
           Strand)

  grouping_column <- "Tumor_Sample_Barcode"
}

if (mode == "HotMAPS") {
  if (grepl("38", subset_input$NCBI_Build[1])) {
    cat("Requested mode is HotMAPS, but the supplied file is in the hg38-based coordinates.\n")
    cat("Unfortunately, HotMAPS is configured to only work for grch37-based maf files.\n")
    stop("Please supply the mutation data in grch37-based version.")
  }

  grouping_column <-"Tumor_Sample_Barcode"

  subset_input <- subset_input %>% unique()
}

if (mode == "gistic2") {
  # Sort by chrom, start, end
  subset_input <- subset_input %>%
    arrange(ID, chrom, start, end)

  # Filter to only canonical chromosomes -------------------
  cat("Filtering to only canonical chromosomes... \n")
  if (projection %in% "hg38"){
    subset_input <- subset_input %>%
      filter(str_detect(chrom, regex("chr[XY\\d]+$", ignore_case = TRUE)))
  } else if (projection %in% "grch37"){
    subset_input <- subset_input %>%
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
            } else if (seg$start[i-1] == seg$start[i]) {
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
            } else if (seg$region_size[i] < seg$region_size[i-1]) {
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
            } else if (seg$region_size[i] > seg$region_size[i-1]) {
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

  subset_input_checked <- check_overlap(subset_input)
  # if no overlaps, do not run solve function
  if ("overlap" %in% subset_input_checked$overlap_status) {
    subset_input <- solve_overlap(subset_input_checked)
  } else {
    subset_input <- subset_input_checked %>%
    select(-overlap_status, -region_size)
  }

  subset_input <- subset_input_checked
}

# Write out appropriate files based on inputs -------------------
if ("maf" %in% names(snakemake@input)){
  # Prepare maf file contents for documentation purposes
  contents <- subset_input %>%
    group_by(across(all_of(grouping_column))) %>%
    summarise(N_mutations=n()) %>%
    ungroup %>%
    mutate(non_coding_included=include_non_coding)

  cat("Writing resulting maf to file...\n")
  write_tsv(subset_input, paste0(output_dir, "/", md5sum, ".maf"))

  cat("Writing maf contents grouped by sample to file...\n")
  write_tsv(contents, paste0(output_dir, "/", md5sum, ".maf.content"))

} else if ("seg" %in% names(snakemake@input)){
  cat("Writing combined seg data to file... \n")
  write_tsv(subset_input, paste0(output_dir, "/", md5sum, ".seg"))

  cat("Writing sample ids to file... \n")
  write_tsv(data.frame(final_sample_set), paste0(output_dir, "/", md5sum, "_sample_ids.txt"))
}

# Writing empty file for snakemake checkpoint rule output
file.create(paste0(output_dir, "/done"))

cat("DONE!")
sink()
