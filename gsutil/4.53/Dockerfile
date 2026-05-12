FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"
RUN mamba create --yes --name env \
        --channel conda-forge \
        "python=3.8" \
        "gsutil=4.53" \
    && mamba clean --all --yes
ENV PATH=/opt/conda/envs/env/bin:$PATH
CMD ["/bin/bash"]
