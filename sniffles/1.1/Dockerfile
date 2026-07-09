FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"
RUN mamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        "sniffles=2.7.5" \
    && rm -rf /opt/conda/pkgs/*
