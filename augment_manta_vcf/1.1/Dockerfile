FROM condaforge/mambaforge:24.9.2-0
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        "cyvcf2=0.34.0" \
        "pyranges=0.1.7" \
    && mamba clean --all --yes

CMD ["/bin/bash"]
