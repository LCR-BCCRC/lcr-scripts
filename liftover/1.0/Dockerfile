FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        ucsc-liftover \
        perl \
        curl \
        bc \
    && rm -rf /opt/conda/pkgs/*

CMD ["/bin/bash"]
