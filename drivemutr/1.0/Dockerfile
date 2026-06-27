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
        r-broom \
        r-ggbeeswarm \
        r-iml \
        r-pre \
        r-remotes \
        r-rlang \
        bioconductor-biocparallel \
        bioconductor-biostrings \
        bioconductor-bsgenome.hsapiens.ucsc.hg19 \
        bioconductor-genomicranges \
        bioconductor-iranges \
        bioconductor-motifbreakr \
        bioconductor-motifdb \
        bioconductor-rtracklayer \
        git \
        make \
        c-compiler \
        cxx-compiler \
        fortran-compiler \
        pkg-config \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ARG GAMBLR_DATA_REF=master
RUN Rscript -e "remotes::install_github( \
        'morinlab/GAMBLR.data', \
        ref = '${GAMBLR_DATA_REF}', \
        dependencies = NA, \
        upgrade = 'never' \
    )"
CMD ["/bin/bash"]