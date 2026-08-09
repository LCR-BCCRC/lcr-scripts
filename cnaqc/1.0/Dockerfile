FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        "r-base=4.3" \
        r-remotes \
        r-optparse \
        r-crayon \
        r-ggpubr \
        r-ggrepel \
        r-vcfr \
        r-clisymbols \
        r-rcolorbrewer \
        r-gtools \
        r-magrittr \
        r-dplyr \
        r-ggplot2 \
        r-cli \
        r-tibble \
        r-tidyr \
        r-akima \
        r-cowplot \
        r-ggsci \
        r-progress \
        r-readr \
        r-r.utils \
        r-scales \
        r-data.table \
        r-biocmanager \
        r-sass \
        bioconductor-rhtslib \
        bioconductor-rsamtools \
        bioconductor-genomicranges \
        bioconductor-variantannotation \
        bioconductor-annotationdbi \
        bioconductor-complexheatmap \
        zlib \
        libuv \
        libcurl \
        openssl \
        libxml2 \
        libgit2 \
        jq \
    && rm -rf /opt/conda/pkgs/*

RUN R --vanilla -q -e 'remotes::install_github("caravagnalab/CNAqc", lib="/opt/conda/lib/R/library", dependencies = TRUE, upgrade = "never")'

CMD ["/bin/bash"]
