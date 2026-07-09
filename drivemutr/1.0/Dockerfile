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
        r-biocmanager \
        bioconductor-biocparallel=1.40.0 \
        bioconductor-biostrings=2.74.0 \
        bioconductor-bsgenome.hsapiens.ucsc.hg19=1.4.3 \
        bioconductor-genomicranges \
        bioconductor-iranges \
        bioconductor-motifdb \
        bioconductor-rtracklayer \
        git \
        jq \
        make \
        c-compiler \
        cxx-compiler \
        fortran-compiler \
        pkg-config \
    && micromamba clean --all --yes
    
RUN micromamba run --name base R --vanilla -q -e 'BiocManager::install("motifbreakR", version = "3.20", ask = FALSE, update = FALSE)'

RUN micromamba run --name base R --vanilla -q -e 'remotes::install_github("morinlab/GAMBLR.data", lib="/opt/conda/lib/R/library", dependencies = TRUE, upgrade = "never")'

CMD ["/bin/bash"]