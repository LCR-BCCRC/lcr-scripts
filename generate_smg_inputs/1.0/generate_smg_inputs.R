#!/usr/bin/env Rscript

# Description:
#   Adapted from generate_smg_inputs/1.0/generate_smg_inputs.R.
#   This script is intended for use with the SMG modules in LCR-modules (MutSig2CV, dNdS, HotMAPS,
#   OncodriveFML).
#   It expects to be run as part of a snakemake workflow which provides a file with categories
#   to subset the metadata in order to get samples IDs of interest. The snakemake workflow will also provide
#   seq_type, launch_date, and output directory values.

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

include_non_coding <- snakemake@params[["include_non_coding"]]
mode <- snakemake@params[["mode"]]
meta <- snakemake@params[['metadata']]
meta_cols <- snakemake@params[['metadata_cols']]

cat("Arguments from snakemake...\n")
cat(paste("Sample sets file:", subsetting_categories_file, "\n"))
cat(paste("Output directory:", output_dir, "\n"))
cat(paste("Sample set:", sample_set, "\n"))
cat(paste("Launch date:", launch_date, "\n"))
cat(paste("Include non-coding:", include_non_coding, "\n"))
cat(paste("Mode:", mode, "\n"))

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
      meta <- meta %>%
          filter(.data[[col]] %in% categories[[col]])
  }

  samples <- meta %>%
    pull(sample_id)

  return(samples)
}

# Get sample ids of the sample_set
this_subset_samples <- subset_samples(subsetting_values, metadata)

# Load master mafs and get mutations for the case set-------------------
maf_files <- snakemake@input[["maf"]]
if ("genome" %in% subsetting_values$seq_type && !("capture" %in% subsetting_values$seq_type)) { # genome only
  cat("Loading genome maf...\n")
  subset_maf <- suppressMessages(read_tsv(maf_files[str_detect(maf_files, "genome")])) %>%
    filter(Tumor_Sample_Barcode %in% this_subset_samples)

} else if (!("genome" %in% subsetting_values$seq_type) && "capture" %in% subsetting_values$seq_type) { # capture only
  cat("Loading capture maf...\n")
  subset_maf <- suppressMessages(read_tsv(maf_files[str_detect(maf_files, "capture")])) %>%
    filter(Tumor_Sample_Barcode %in% this_subset_samples)

} else if ("genome" %in% subsetting_values$seq_type && "capture" %in% subsetting_values$seq_type) { # both
  cat("Loading genome maf...\n")
  genome_maf <- suppressMessages(read_tsv(maf_files[str_detect(maf_files, "genome")])) %>%
    filter(Tumor_Sample_Barcode %in% this_subset_samples)

  cat("Loading capture maf...\n")
  capture_maf <- suppressMessages(read_tsv(maf_files[str_detect(maf_files, "capture")])) %>%
    filter(!Tumor_Sample_Barcode %in% unique(genome_maf$Tumor_Sample_Barcode)) %>%
    filter(Tumor_Sample_Barcode %in% this_subset_samples)

  subset_maf <- rbind(genome_maf, capture_maf)
}

# subset for coding only if user requested
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
    subset_maf =
      subset_maf %>%
      dplyr::filter(Variant_Classification %in% coding_class)
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
missing_samples = setdiff(this_subset_samples,
                          unique(subset_maf$Tumor_Sample_Barcode))

if (length(missing_samples)==0) {
  cat(paste("Success! Found mutations for all samples.", length(this_subset_samples), "patients will be used in the analysis\n"))
  md5sum <- digest(this_subset_samples)
  final_sample_set <- this_subset_samples
} else {
  cat(paste("WARNING:", length(missing_samples), "will not be available for the analysis.\n"))
  cat("Writing missing sample ids to file... \n")
  final_sample_set <- unique(subset_maf$Tumor_Sample_Barcode)
  md5sum <- digest(final_sample_set)
  write_tsv(data.frame(missing_samples), paste0(output_dir, "/", md5sum, "_missing_sample_ids.txt"))
}

# Format maf according to the requirements of each individual tool --------------------------------------
cat("preparing maf file to be used with", mode, "\n")
# MutSig2CV
if (mode == "MutSig2CV") {
  # check that the maf file is in grch37-based coordinates
  if (grepl("38", subset_maf$NCBI_Build[1])) {
    cat("Requested mode is MutSig2CV, but the supplied file is in the hg38-based coordinates.\n")
    cat("Unfortunatelly, MutSig only works for grch37-based maf files.\n")
    stop("Please supply the mutation data in grch37-based version.")
  }

  subset_maf =
    subset_maf %>%
    rename("chr"="Chromosome",
           "pos"="Start_Position",
           "gene"="Hugo_Symbol",
           "patient"="Tumor_Sample_Barcode",
           "ref_allele"="Reference_Allele",
           "newbase"="Tumor_Seq_Allele2",
           "type"="Variant_Classification",
           "classification"="Variant_Type")

subset_maf =
  subset_maf %>%
    select(chr, pos, gene, patient, ref_allele, newbase, type, classification)

  grouping_column = "patient"

}

# dNdS
if (mode == "dNdS") {
  # check that the maf file is in grch37-based coordinates
  if (grepl("38", subset_maf$NCBI_Build[1])) {
    cat("Requested mode is dNdS, but the supplied file is in the hg38-based coordinates.\n")
    cat("Unfortunatelly, dNdS is configured to only work for grch37-based maf files.\n")
    stop("Please supply the mutation data in grch37-based version.")
  }

  subset_maf =
    subset_maf %>%
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

  grouping_column = "sampleID"

}

if (mode == "FishHook") {

  subset_maf =
    subset_maf %>%
    select(Hugo_Symbol,
           Tumor_Sample_Barcode,
           Chromosome,
           Start_Position,
           End_Position,
           Variant_Classification,
           Strand)

  grouping_column = "Tumor_Sample_Barcode"

}

if (mode == "HotMAPS") {
  if (grepl("38", subset_maf$NCBI_Build[1])) {
    cat("Requested mode is HotMAPS, but the supplied file is in the hg38-based coordinates.\n")
    cat("Unfortunately, HotMAPS is configured to only work for grch37-based maf files.\n")
    stop("Please supply the mutation data in grch37-based version.")
  }

  grouping_column = "Tumor_Sample_Barcode"

  subset_maf = subset_maf %>% unique()
}

# Prepare maf file contents for documentation purposes
contents = subset_maf %>%
  group_by(across(all_of(grouping_column))) %>%
  summarise(N_mutations=n()) %>%
  ungroup %>%
  mutate(non_coding_included=include_non_coding)

# Write out final maf file -------------------
cat("Writing resulting maf to file...\n")
write_tsv(subset_maf, paste0(output_dir, "/", md5sum, ".maf"))

# Write out contents file -------------------
cat("Writing maf contents grouped by sample to file...\n")
write_tsv(contents, paste0(output_dir, "/", md5sum, ".maf.content"))

# Writing empty file for snakemake checkpoint rule output
file.create(paste0(output_dir, "/done"))

cat("DONE!")
sink()
