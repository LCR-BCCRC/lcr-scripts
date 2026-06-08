FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        "vcf2maf=1.6.18" \
        "samtools" \
    && mamba clean -afy

CMD ["/bin/bash"]
