FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        "mosdepth=0.3.13" \
    && rm -rf /opt/conda/pkgs/*

CMD ["/bin/bash"]
