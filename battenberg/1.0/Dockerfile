FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        r-base \
        r-devtools \
        "impute2=2.3.2" \
        cancerit-allelecount \
        parallel \
        perl \
    && rm -rf /opt/conda/pkgs/*

RUN R --vanilla -q -e 'devtools::install_github("morinlab/ascat/ASCAT", dependencies = TRUE, upgrade = "never")' \
    && R --vanilla -q -e 'devtools::install_github("morinlab/battenberg", dependencies = FALSE, upgrade = "never")'

CMD ["/bin/bash"]
