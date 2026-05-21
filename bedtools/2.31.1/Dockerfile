FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        "bedtools=2.31.1" \
    && rm -rf /opt/conda/pkgs/*

CMD ["/bin/bash"]
