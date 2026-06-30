FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        "r-base=4.3" \
        r-remotes \
        r-optparse \
        r-dplyr \
        r-ggplot2 \
        r-data.table \
        r-ggpubr \
        r-biocmanager \
        r-sass \
        r-tidyr \
        r-stringr \
        r-purrr \
        r-cowplot \
        r-patchwork \
        r-scales \
        r-magrittr \
        r-htmltools \
        r-knitr \
        r-rmarkdown \
        r-curl \
        r-httr \
        r-xml2 \
        r-jsonlite \
        r-crayon \
        r-cli \
        r-glue \
        r-fs \
        r-withr \
        zlib \
        libuv \
        libcurl \
        openssl \
        libxml2 \
        libgit2 \
    && rm -rf /opt/conda/pkgs/*

RUN R --vanilla -q -e 'remotes::install_github("caravagnalab/CNAqc", lib="/opt/conda/lib/R/library", dependencies = TRUE, upgrade = "never")'

CMD ["/bin/bash"]
