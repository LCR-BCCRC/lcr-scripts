FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

# pysam=0.18.0 requires Python <3.10 which conflicts with the Python 3.12
# pin in the mambaforge base; install without a version pin so mamba picks
# a Python-3.12-compatible pysam (>=0.18.0).
# oncopipe==1.0.11 declares snakemake<5.19 which fails on Python 3.12;
# install with --no-deps since the container only needs oncopipe's Python
# utilities (op.absolute_symlink), not a full snakemake installation.
RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        "pysam>=0.18.0" \
        "pandas" \
    && pip install "oncopipe==1.0.11" --no-deps \
    && mamba clean --all --yes

CMD ["/bin/bash"]
