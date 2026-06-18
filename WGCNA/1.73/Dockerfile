FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"
RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        r-base \
        r-readr \
        r-dplyr \
        r-tibble \
        r-ggplot2 \
        r-tidyr \
        r-matrixstats \
        r-remotes \
        bioconductor-deseq2 \
        bioconductor-limma \
        bioconductor-biocparallel \
        bioconductor-wgcna \
    && mamba clean --all --yes

ARG GITHUB_PAT=""
RUN Rscript -e "remotes::install_github('morinlab/GAMBLR.helpers', ref = 'master', upgrade = 'never')"
CMD ["/bin/bash"]