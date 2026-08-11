FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

# Install R and all conda-available GAMBLR dependencies.
# Channels match gamblr.yaml (conda-forge + bioconda, no defaults).
RUN mamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        "r-base=4.1.3" \
        "r-nloptr=1.2.2.3" \
        "r-tidyverse" \
        "r-data.table" \
        "r-devtools" \
        "r-glue" \
        "r-config" \
        "r-dbi" \
        "r-rmariadb" \
        "r-knitr" \
        "r-rmarkdown" \
        "r-workflowr" \
        "r-cowplot" \
        "r-ggbeeswarm" \
        "r-ggpubr" \
        "r-ggrepel" \
        "r-ggthemes" \
        "r-plotly" \
        "r-rcircos" \
        "r-wordcloud" \
        "r-biomartr" \
        "r-testthat" \
        "bioconductor-complexheatmap" \
        "bioconductor-maftools" \
        "bioconductor-rtracklayer" \
        "bioconductor-sradb" \
    && mamba clean --all --yes

# Install GAMBLR from the commit SHA pinned in the vcf2maf 1.3 module config.
# dependencies=TRUE installs CRAN-only deps (g3viz, metaviz) that are not
# available via conda.  The library() call makes the layer fail loudly if
# the install silently completed with a non-zero R CMD INSTALL exit status.
RUN R --vanilla -q -e \
    'remotes::install_github("morinlab/GAMBLR@d8d643face5b67bdfe850aabcbfe26f82254a825", \
        lib="/opt/conda/lib/R/library", \
        dependencies=TRUE, \
        upgrade="never")' \
    && R --vanilla -q -e 'library(GAMBLR)'

CMD ["/bin/bash"]
