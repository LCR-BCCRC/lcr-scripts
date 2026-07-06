FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        "hmftools-sage=5.0" \
        htslib \
    && rm -rf /opt/conda/pkgs/*

CMD ["/bin/bash"]
