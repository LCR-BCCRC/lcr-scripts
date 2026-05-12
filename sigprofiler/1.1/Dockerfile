FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"
RUN mamba install --yes --name base \
        --channel conda-forge \
        "python=3.9" \
    && mamba clean --all --yes \
    && pip install --no-cache-dir \
        "sigprofilerextractor==1.1.3" \
        "sigprofilermatrixgenerator==1.1.31" \
        "sigprofilerplotting==1.1.17" \
        "torch==1.9.0"
CMD ["/bin/bash"]
