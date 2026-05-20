FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"
RUN apt-get update -q && apt-get install -y --no-install-recommends gcc libc6-dev && rm -rf /var/lib/apt/lists/* \
    && mamba create --yes --name env \
        --channel conda-forge \
        --channel bioconda \
        "python=3.7" \
        "pip" \
        "htslib" \
        "curl" \
    && mamba clean --all --yes \
    && /opt/conda/envs/env/bin/pip install --no-cache-dir \
        "oncodrivefml==2.4.0"
ENV PATH=/opt/conda/envs/env/bin:$PATH
CMD ["/bin/bash"]
