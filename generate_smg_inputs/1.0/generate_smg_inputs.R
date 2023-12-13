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
subsetting_categories_file <- snakemake@input[["sample_sets"]]
output_dir <- dirname(snakemake@output[[1]])

case_set <- snakemake@wildcards[["sample_set"]]
launch_date <- snakemake@wildcards[["launch_date"]]

seq_type <- unlist(snakemake@params[["seq_type"]])
include_non_coding <- snakemake@params[["include_non_coding"]]
mode <- snakemake@params[["mode"]]

cat("Arguments from snakemake...\n")
cat(paste("Sample sets file:", subsetting_categories_file, "\n"))
cat(paste("Output directory:", output_dir, "\n"))
cat(paste("Sample set:", case_set, "\n"))
cat(paste("Launch date:", launch_date, "\n"))
cat(paste("Seq type(s):", seq_type, "\n"))
cat(paste("Include non-coding:", include_non_coding, "\n"))
cat(paste("Mode:", mode, "\n"))



# pandas df from snakemake is passed as a list of lists object
# This converts the lists to columns of a dataframe
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

# Get samples in subset and ensure consistent naming of the sample ID column-------------------
if (file.exists(subsetting_categories_file)) {
  full_subsetting_categories <- suppressMessages(read_tsv(subsetting_categories_file, comment="#"))
} else {
  cat(paste("Warning: case set is requested, but the subsetting categories file", subsetting_categories_file, "is not found.\n"))
  stop("Exiting because did not find subsetting categories file")
}

# Get subsetting values for this sample_set
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
this_subset_samples <- subset_samples(subsetting_values, metadata)

# Load master mafs and get mutations for the case set-------------------
maf_files <- snakemake@input[["maf"]]
if ("genome" %in% seq_type && !("capture" %in% seq_type)) { # genome only
  cat("Loading genome maf...\n")
  subset_maf <- suppressMessages(read_tsv(maf_files[str_detect(maf_files, "genome")])) %>%
    filter(Tumor_Sample_Barcode %in% this_subset_samples)

} else if (!("genome" %in% seq_type) && "capture" %in% seq_type) { # capture only
  cat("Loading capture maf...\n")
  subset_maf <- suppressMessages(read_tsv(maf_files[str_detect(maf_files, "capture")])) %>%
    filter(Tumor_Sample_Barcode %in% this_subset_samples)

} else if ("genome" %in% seq_type && "capture" %in% seq_type) { # both
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

# report any missing samples and calculate md5sum for samples
missing_samples = setdiff(this_subset_samples,
                          unique(subset_maf$Tumor_Sample_Barcode))

if (length(missing_samples)==0) {
  cat(paste("Success! Found mutations for all samples.", length(this_subset_samples), "patients will be used in the analysis\n"))
  md5sum <- digest(this_subset_samples)
  final_sample_set <- this_subset_samples
} else {
  cat(paste("WARNING:", length(missing_samples), "will not be available for the analysis.\n"))
  cat("Did not find mutations for these samples in the master input maf:\n")
  cat(missing_samples)
  cat("\n")
  final_sample_set <- subset_maf %>% pull(Tumor_Sample_Barcode)
  md5sum <- digest(final_sample_set)
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


# Prepare maf file contents for documentation purposes
contents = subset_maf %>%
  group_by(across(all_of(grouping_column))) %>%
  summarise(N_mutations=n()) %>%
  ungroup %>%
  mutate(non_coding_included=include_non_coding)


# Output --------------------------------------------------------------------------------------
# Check if output dir extists, create if not
if (!dir.exists(file.path(output_dir))){
  cat("Output directory for sample_set and launch date combo does not exist. Creating it...\n")
  cat(output_dir,"\n")
  dir.create(file.path(output_dir), recursive = TRUE)
} else {
  cat("Output directory for sample_set and launch date combo exists.\n")
  cat(output_dir,"\n")
}

cat("Writing resulting maf to file...\n")
write_tsv(subset_maf, paste0(output_dir, "/", md5sum, ".maf"))
write_tsv(contents, paste0(output_dir, "/", md5sum, ".maf.content"))

cat("Writing samples corresponding to md5sum to file...\n")
write_tsv(data.frame(final_sample_set), paste0(output_dir, "/", md5sum, "_sample_ids.txt"))

# Writing empty file for snakemake checkpoint rule output
file.create(paste0(output_dir, "/done"))

cat("DONE!")
sink()
