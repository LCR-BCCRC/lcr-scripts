FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"
RUN mamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        'r-base>=4.1' \
        r-readr \
        r-dplyr \
        r-purrr \
        r-cluster \
        r-ggplot2 \
    && mamba clean --all --yes
CMD ["/bin/bash"]
