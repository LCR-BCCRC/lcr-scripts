# /usr/bin/env Rscript

# Set up log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

# Load packages -----------------------------------------------------------
cat("Loading packages... \n")
suppressWarnings(
    suppressPackageStartupMessages({
        library(tidyverse)
    })
)

# Interpret snakemake variables
grch37_filled_file <- snakemake@input[["grch37_filled"]]
hg38_filled_file <- snakemake@input[["hg38_filled"]]

output_grch37_file <- snakemake@output[["grch37_reviewed"]]
output_hg38_file <- snakemake@output[["hg38_reviewed"]]

review_dummy_segments <- function(
    input_file,
    output_file,
    projection
){
    incoming_seg <- suppressMessages(
        read_tsv(
            input_file
        )
    )

    message(paste0('The sample id is ', incoming_seg$ID[1]))
    message(paste0('The projection is ', projection))

    original_colnames <- colnames(incoming_seg)

    incoming_seg <- incoming_seg %>%
        mutate(
            CN = 2*2^log.ratio
        )

    message(
        'Inferring average log ratio for segments where dummy_segment is not 1 ...'
    )
    average_log <- incoming_seg %>%
        filter(!dummy_segment == 1) %>%
        mutate(
            length = end - start + 1,
            CN_seg = CN * length,
            logr_seg = log.ratio * length
        ) %>%
        group_by(ID) %>%
        summarise(
            mean = mean(CN),
            real_mean = sum(CN_seg) / sum(length),
            real_mean_logr = sum(logr_seg) / sum(length)
        ) # actual average per base

    message(paste0('The mean CN is ', average_log$mean))
    message(paste0('The mean log.ratio is ', average_log$real_mean_logr))

    incoming_seg <- left_join(
        incoming_seg,
        average_log %>%
            select(ID, real_mean, real_mean_logr),
            by = "ID"
    )

    reviewed_seg <- incoming_seg %>%
        mutate(
            log.ratio = ifelse(
                dummy_segment == 1,
                real_mean_logr,
                log.ratio
            )
        ) %>%
        select(
            all_of(original_colnames)
        )

    reviewed_seg %>%
        write_tsv(
            output_file
        )
}

review_dummy_segments(
    grch37_filled_file,
    output_grch37_file,
    "grch37"
)

review_dummy_segments(
    hg38_filled_file,
    output_hg38_file,
    "hg38"
)

message("DONE!")
sink()
