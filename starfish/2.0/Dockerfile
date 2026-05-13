FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        libxml2 \
        pysam \
        matplotlib \
        bcftools \
        rtg-tools \
    && mamba clean --all --yes

CMD ["/bin/bash"]
