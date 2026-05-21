FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        "r-sequenza=3.0.0" \
        r-remotes \
    && rm -rf /opt/conda/pkgs/*

# Pre-install the hg38-compatible GitHub fork of copynumber so the
# runtime check in run_sequenza.R is satisfied without network access.
RUN R --vanilla -q -e \
    'remotes::install_github("ShixiangWang/copynumber", lib="/opt/conda/lib/R/library", dependencies=TRUE, upgrade="never")'

CMD ["/bin/bash"]
