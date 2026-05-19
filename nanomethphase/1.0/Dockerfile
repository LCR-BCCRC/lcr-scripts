FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        "python=3.8" \
        "numpy=1.18.5" \
        "samtools=1.21" \
        "bioconductor-dss=2.48.0" \
        "r-base=4.3.3" \
        "r-sys" \
    && pip install nanomethphase==1.0 \
    && rm -rf /opt/conda/pkgs/*

CMD ["/bin/bash"]
