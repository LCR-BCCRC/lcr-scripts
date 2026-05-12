FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"
RUN mamba create --yes --name env \
        --channel bbglab \
        --channel defaults \
        "python=3.6" \
        "oncodriveclustl=1.1.3" \
    && mamba clean --all --yes
ENV PATH=/opt/conda/envs/env/bin:$PATH
ENV HOME=/tmp
CMD ["/bin/bash"]
