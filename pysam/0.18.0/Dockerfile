FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

# pysam=0.18.0 requires Python <3.10; install without a version pin so
# mamba picks a Python-3.12-compatible version (>=0.18.0).
# oncopipe==1.0.11 imports both yaml and snakemake at module load time and
# declares snakemake<5.19 which fails on Python 3.12.  Install
# snakemake-minimal via conda (any modern version satisfies the import),
# then oncopipe with --no-deps so pip does not enforce the <5.19 constraint.
RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        "pysam>=0.18.0" \
        "pandas" \
        "pyyaml" \
        "snakemake-minimal" \
    && pip install "oncopipe==1.0.11" --no-deps \
    && mamba clean --all --yes

CMD ["/bin/bash"]
