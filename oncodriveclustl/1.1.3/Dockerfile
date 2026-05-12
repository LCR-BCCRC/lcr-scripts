FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"
RUN mamba install --yes --name base \
        --channel bbglab \
        --channel defaults \
        "oncodriveclustl=1.1.3" \
    && mamba clean --all --yes
CMD ["/bin/bash"]
