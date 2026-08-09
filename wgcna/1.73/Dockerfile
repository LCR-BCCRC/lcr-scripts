FROM mambaorg/micromamba:2.8.1

LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN micromamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        --strict-channel-priority \
        r-base=4.4 \
        r-tidyverse \
        r-matrixstats \
        r-remotes \
        r-config \
        r-ggthemes \
        r-magrittr \
        r-stringr \
        bioconductor-deseq2=1.46.0 \
        bioconductor-limma \
        bioconductor-biocparallel=1.40.0 \
        r-wgcna=1.73 \
        git \
        jq \
        make \
        c-compiler \
        cxx-compiler \
        fortran-compiler \
        pkg-config \
    && micromamba clean --all --yes

RUN micromamba run --name base R --vanilla -q -e 'remotes::install_github("morinlab/GAMBLR.helpers@fdde14c", lib="/opt/conda/lib/R/library", dependencies = TRUE, upgrade = "never")'
CMD ["/bin/bash"]