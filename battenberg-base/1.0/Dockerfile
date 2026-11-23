# Battenberg dependencies without Battenberg itself. Not used by any module directly: it is the base
# for Battenberg test images, which add only the Battenberg commit under test (see ../README.md).
#
# Everything is pinned to what battenberg:1.2 and battenberg:1.3-test were built with (Nov 2026), so a
# test image differs from those only in the Battenberg commit:
#   base image digest, conda package versions, ASCAT commit, and CRAN packages from a dated snapshot.
# Change any pin only in a new version directory (battenberg-base/1.1), never in place.
FROM condaforge/mambaforge:latest@sha256:050753d138b6708128c9bb45055a9d943f06363f32f63a96f7dd5f304ada52a6
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        "r-base=4.5.3" \
        "r-remotes=2.5.0" \
        "r-optparse=1.8.2" \
        "impute2=2.3.2" \
        "cancerit-allelecount=4.3.0" \
        "perl=5.44.0" \
        "parallel=20260922" \
    && rm -rf /opt/conda/pkgs/*

# CRAN packages come from a dated Posit Package Manager snapshot (the day before battenberg:1.2 was
# built). The list is Battenberg's Imports from its DESCRIPTION, minus ASCAT and base-R packages.
RUN --mount=type=secret,id=GITHUB_PAT \
    if [ -s /run/secrets/GITHUB_PAT ]; then export GITHUB_PAT="$(cat /run/secrets/GITHUB_PAT)"; fi \
    && R --vanilla -q -e 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2026-11-21")); \
        remotes::install_github("morinlab/ascat/ASCAT@e3d6577909e059aa0aa873c08dd56c417abdf93a", lib = "/opt/conda/lib/R/library", dependencies = TRUE, upgrade = "never"); \
        install.packages(c("RColorBrewer", "ggplot2", "readr", "gtools", "gridExtra", "doParallel", "foreach"), lib = "/opt/conda/lib/R/library"); \
        missing <- setdiff(c("ASCAT", "RColorBrewer", "ggplot2", "readr", "gtools", "gridExtra", "doParallel", "foreach"), rownames(installed.packages())); \
        if (length(missing)) stop("not installed: ", paste(missing, collapse = " "))'

CMD ["/bin/bash"]
