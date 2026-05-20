FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN apt-get update -q && apt-get install -y --no-install-recommends gcc && rm -rf /var/lib/apt/lists/* \
    && mamba create --yes --name env \
        --channel bioconda \
        --channel conda-forge \
        "python=3.8" \
        "numpy=1.18.5" \
        "samtools=1.21" \
        "bioconductor-dss=2.48.0" \
        "r-base=4.3.3" \
        "r-sys" \
        "pip" \
    && mamba clean --all --yes \
    && /opt/conda/envs/env/bin/pip install --no-cache-dir nanomethphase==1.0
ENV PATH=/opt/conda/envs/env/bin:$PATH
CMD ["/bin/bash"]
