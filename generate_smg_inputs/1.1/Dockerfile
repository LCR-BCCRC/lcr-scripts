FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"
RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        "r-base=4.1.3" \
        "r-data.table=1.14.2" \
        "r-tidyverse=1.3.1" \
        "r-digest=0.6.29" \
    && mamba clean --all --yes
CMD ["/bin/bash"]
