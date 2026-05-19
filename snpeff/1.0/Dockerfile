FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        "snpeff=4.3.1t" \
        "bcftools=1.10.2" \
    && rm -rf /opt/conda/pkgs/*

CMD ["/bin/bash"]
