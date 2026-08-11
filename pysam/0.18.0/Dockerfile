FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

# pysam=0.18.0 requires Python <3.10 which conflicts with the Python 3.12
# pin in the mambaforge base.  Install without a version pin so mamba picks
# a Python-3.12-compatible pysam (>=0.18.0).
RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        "pysam>=0.18.0" \
        "pandas" \
    && pip install "oncopipe==1.0.11" \
    && mamba clean --all --yes

CMD ["/bin/bash"]
