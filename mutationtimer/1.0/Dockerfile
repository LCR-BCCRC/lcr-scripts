FROM rocker/tidyverse:4.2
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

# Install system dependencies?
# Not sure what these are or if it's required since we are starting from the tidyverse image

# Install CRAN packages
RUN install2.r --error \
    data.table \
    xml \
    ggpubr \

# Install BioConductor packages
RUN R -e "BiocManager::install('VariantAnnotation')" \
    && R -e "BiocManager::install('BSgenome')" \
    && R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')" \
    && R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')" \

# Install remote packages
RUN R -e 'options(timeout=9999999); devtools::install_github("mg14/mg14", dependencies = TRUE, upgrade = "never")' \
    && R -e 'options(timeout=9999999); devtools::install_github("gerstung-lab/MutationTimeR@e4e266a", dependencies = TRUE, upgrade = "never")' \
    && R -e 'options(timeout=9999999); devtools::install_github("morinlab/GAMBLR.helpers@v1.3", dependencies = TRUE, upgrade = "never")' \
    && R -e 'options(timeout=9999999); install.packages("", dependencies = TRUE, upgrade = "never")'