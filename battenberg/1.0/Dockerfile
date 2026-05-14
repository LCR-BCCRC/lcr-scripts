FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        r-base \
        r-remotes \
        r-optparse \
        "impute2=2.3.2" \
        cancerit-allelecount \
        parallel \
        perl \
    && rm -rf /opt/conda/pkgs/*

RUN R --vanilla -q -e 'remotes::install_version("readr", version="1.4.0", repos="https://cran.r-project.org", lib="/opt/conda/lib/R/library", upgrade="never")' \
    && R --vanilla -q -e 'remotes::install_github("morinlab/ascat/ASCAT", lib="/opt/conda/lib/R/library", dependencies = TRUE, upgrade = "never")' \
    && R --vanilla -q -e 'remotes::install_github("morinlab/battenberg", lib="/opt/conda/lib/R/library", dependencies = TRUE, upgrade = "never")'

CMD ["/bin/bash"]
