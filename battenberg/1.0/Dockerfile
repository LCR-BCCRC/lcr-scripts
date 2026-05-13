FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        --channel r \
        "r-base=3.6.3" \
        r-devtools \
        r-optparse \
        r-biocmanager \
        r-ggplot2 \
        r-doparallel \
        r-foreach \
        r-gridextra \
        r-reshape2 \
        r-gtools \
        r-rcolorbrewer \
        "impute2=2.3.2" \
        "cancerit-allelecount=4.0.2" \
        parallel \
        perl \
    && rm -rf /opt/conda/pkgs/*

RUN R --vanilla -q -e 'devtools::install_github("morinlab/ascat/ASCAT", dependencies = TRUE, upgrade = "never")' \
    && R --vanilla -q -e 'devtools::install_github("morinlab/battenberg", dependencies = FALSE, upgrade = "never")'

CMD ["/bin/bash"]
