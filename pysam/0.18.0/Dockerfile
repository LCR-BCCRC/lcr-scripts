FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        "pysam=0.18.0" \
        "pandas" \
    && pip install "oncopipe==1.0.11" \
    && mamba clean --all --yes

CMD ["/bin/bash"]
