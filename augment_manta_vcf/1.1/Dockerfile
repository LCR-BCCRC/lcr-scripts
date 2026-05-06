FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        "cyvcf2=0.20.0" \
        "pyranges=0.0.77" \
        "numpy=1.18.4" \
    && mamba clean --all --yes

CMD ["/bin/bash"]
