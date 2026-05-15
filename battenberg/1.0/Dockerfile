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

# morinlab/battenberg updated 2026-05-15: fixed length>1 if() condition in
# segmentation.R that was silently permissive under R 3.x but a hard error in R 4.4+
RUN R --vanilla -q -e 'remotes::install_github("morinlab/ascat/ASCAT", lib="/opt/conda/lib/R/library", dependencies = TRUE, upgrade = "never")' \
    && R --vanilla -q -e 'remotes::install_github("morinlab/battenberg", lib="/opt/conda/lib/R/library", dependencies = TRUE, upgrade = "never")'

CMD ["/bin/bash"]
