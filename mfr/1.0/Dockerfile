FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"
RUN mamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        "r-base>=4.1" \
        r-readr \
        r-dplyr \
        r-cluster \
        r-ggplot2 \
        "python>=3.9" \
        htslib \
    && rm -rf /opt/conda/pkgs/*
CMD ["/bin/bash"]
