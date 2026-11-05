FROM rocker/tidyverse:4.2
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

# Install system dependencies
RUN apt-get update -q && apt-get install -y --no-install-recommends \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN install2.r --error \
    data.table \
    XML \
    ggpubr

# Install BioConductor packages
RUN R -e "BiocManager::install('Rhtslib')" \
    && R -e "BiocManager::install('BSgenome')" \
    && R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')" \
    && R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')" \
    && R -e "BiocManager::install('VariantAnnotation')"

# Install remote packages
RUN R -e 'options(timeout=9999999); devtools::install_github("mg14/mg14", dependencies = TRUE, upgrade = "never")' \
    && R -e 'options(timeout=9999999); devtools::install_github("gerstung-lab/MutationTimeR@e4e266a", dependencies = TRUE, upgrade = "never")' \
    && R -e 'options(timeout=9999999); devtools::install_github("morinlab/GAMBLR.helpers@v1.3", dependencies = TRUE, upgrade = "never")'

CMD ["/bin/bash"]