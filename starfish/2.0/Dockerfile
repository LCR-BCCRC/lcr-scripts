FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        pysam \
        matplotlib \
        bcftools \
        rtg-tools \
    && rm -rf /opt/conda/pkgs/*

RUN pip install pyvenn

CMD ["/bin/bash"]
