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

# Use the build's GitHub token, if provided, for authenticated API calls (avoids the anonymous rate limit).
RUN --mount=type=secret,id=GITHUB_PAT \
    if [ -s /run/secrets/GITHUB_PAT ]; then export GITHUB_PAT="$(cat /run/secrets/GITHUB_PAT)"; fi \
    && R --vanilla -q -e 'remotes::install_github("morinlab/ascat/ASCAT", lib="/opt/conda/lib/R/library", dependencies = TRUE, upgrade = "never")' \
    && R --vanilla -q -e 'remotes::install_github("morinlab/battenberg@bbd1b1c", lib="/opt/conda/lib/R/library", dependencies = TRUE, upgrade = "never")'

CMD ["/bin/bash"]
