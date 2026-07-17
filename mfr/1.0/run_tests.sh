#!/usr/bin/env bash
set -euo pipefail

# Golden-file test for the mfr container: exercises the same operations the
# mfR lcr-module's src/python/extract_chrom_maf.py and src/R/cluster_foci.R
# scripts perform (tabix-based per-chromosome MAF extraction + coding-variant
# filtering, then hierarchical clustering of mutation positions into foci),
# using the same tools this Dockerfile installs. The logic below is a
# self-contained mirror of those scripts, not an import of them -- the actual
# script sources live in lcr-modules, a separate repo.
#
# Usage:
#   ./run_tests.sh                      # writes to tests/output (local golden-file generation)
#   ./run_tests.sh tests/docker_output  # used by CI to write Docker outputs for comparison
#
# Generating / updating golden files:
#   Run ./run_tests.sh from within the conda environment and commit tests/output/.

OUTDIR="${1:-tests/output}"
mkdir -p "$OUTDIR"

CHROM="chr1"
CODING_CLASSES="Missense_Mutation|Silent"

# ---- Stage 1: tabix-extract one chromosome's rows from each sample MAF and
# drop coding Variant_Classification rows (mirrors extract_chrom_maf.py) -----
echo "Extracting $CHROM from tests/input/*.maf.gz, dropping coding rows..."

HEADER=$(zcat tests/input/sample1.maf.gz | head -1)
VC_COL=$(echo "$HEADER" | tr '\t' '\n' | grep -nx "Variant_Classification" | cut -d: -f1)

{
    echo "$HEADER"
    for maf in tests/input/*.maf.gz; do
        tabix "$maf" "$CHROM" | awk -F'\t' -v col="$VC_COL" -v classes="$CODING_CLASSES" \
            'BEGIN { split(classes, arr, "|"); for (i in arr) drop[arr[i]] = 1 }
             !($col in drop)'
    done | sort -k3,3n
} > "$OUTDIR/$CHROM.noncoding.maf"

echo "Wrote $OUTDIR/$CHROM.noncoding.maf"

# ---- Stage 2: hierarchical clustering of unique positions into foci --------
# (mirrors cluster_one_chromosome() in cluster_foci.R; same defaults as
# lcr-modules' mfR/1.0/config/default.yaml: euclidean / centroid / h 1..100)
echo "Clustering positions into foci..."

PLOT_SCRATCH=$(mktemp -d)

Rscript --vanilla - "$OUTDIR/$CHROM.noncoding.maf" "$OUTDIR/$CHROM.foci.tsv" "$PLOT_SCRATCH/$CHROM.silhouette.pdf" <<'EOF'
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(cluster)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
in_maf   <- args[[1]]
out_tsv  <- args[[2]]
out_plot <- args[[3]]

pos_col       <- "Start_Position"
dist_method   <- "euclidean"
hclust_method <- "centroid"
h_min <- 1L
h_max <- 100L

maf <- read_tsv(in_maf, show_col_types = FALSE) %>% arrange(.data[[pos_col]])
positions <- sort(unique(maf[[pos_col]]))
n <- length(positions)

if (n < 2) {
  maf$group <- if (n == 0) integer(0) else 1L
  plot_df <- data.frame(h = integer(0), avg = numeric(0))
  best_h <- NA_integer_
} else {
  d  <- dist(positions, method = dist_method)
  hc <- hclust(d, method = hclust_method)

  h_values  <- seq.int(h_min, h_max)
  sil_means <- rep(NA_real_, length(h_values))

  for (i in seq_along(h_values)) {
    lab <- tryCatch(cutree(hc, h = h_values[i]), error = function(e) NULL)
    if (is.null(lab)) next
    # silhouette() needs >1 cluster AND at least one non-singleton cluster --
    # with every point in its own cluster (common at h_min before any merge
    # has happened) it returns a bare NA instead of a matrix, and indexing
    # that with [, "sil_width"] throws "incorrect number of dimensions".
    if (length(unique(lab)) > 1 && length(unique(lab)) < length(lab)) {
      # Belt-and-suspenders: the condition above rules out the two known
      # degenerate cases, but wrap anyway so any other silhouette() edge
      # case can't abort the job either (matches the cutree() tryCatch above).
      sil_means[i] <- tryCatch(
        mean(silhouette(lab, d)[, "sil_width"]),
        error = function(e) NA_real_
      )
    }
  }

  best_h <- if (all(is.na(sil_means))) h_min else h_values[which.max(sil_means)]
  labels <- tryCatch(cutree(hc, h = best_h), error = function(e) seq_along(positions))
  lab_df <- setNames(data.frame(positions, labels), c(pos_col, "group"))
  maf <- left_join(maf, lab_df, by = pos_col)
  plot_df <- data.frame(h = h_values, avg = sil_means)
}

write_tsv(maf, out_tsv, na = "")

p <- ggplot(plot_df, aes(x = h, y = avg)) +
  geom_line() +
  geom_point(shape = 16) +
  labs(x = "Height (h)", y = "Mean silhouette width", title = "mfr container test") +
  theme_minimal()
ggsave(out_plot, p, width = 7, height = 4)

message(sprintf("best_h = %s -> %s", as.character(best_h), out_tsv))
EOF

if [ ! -s "$PLOT_SCRATCH/$CHROM.silhouette.pdf" ]; then
    echo "FAIL: silhouette plot was not created" >&2
    exit 1
fi
echo "PASS: ggplot2/cairo rendering produced a non-empty plot (not compared -- PDFs aren't byte-stable across builds)"
rm -rf "$PLOT_SCRATCH"

echo "Wrote $OUTDIR/$CHROM.foci.tsv"
