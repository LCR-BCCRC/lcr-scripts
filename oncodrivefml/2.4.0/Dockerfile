FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"
RUN mamba create --yes --name env \
        --channel conda-forge \
        "python=3.7" \
        "pip" \
        "htslib" \
        "curl" \
    && mamba clean --all --yes \
    && /opt/conda/envs/env/bin/pip install --no-cache-dir \
        "oncodrivefml==2.4.0"
ENV PATH=/opt/conda/envs/env/bin:$PATH
CMD ["/bin/bash"]
