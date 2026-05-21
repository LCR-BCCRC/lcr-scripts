FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        "r-sequenza=3.0.0" \
        r-remotes \
    && rm -rf /opt/conda/pkgs/*

# Pre-install the hg38-compatible GitHub fork of copynumber so the
# runtime check in run_sequenza.R is satisfied without network access.
# Use install_git (plain git clone) to avoid the GitHub API rate limit
# that install_github hits when resolving HEAD without a PAT.
RUN R --vanilla -q -e \
    'remotes::install_git("https://github.com/ShixiangWang/copynumber.git", lib="/opt/conda/lib/R/library", dependencies=TRUE, upgrade="never")'

CMD ["/bin/bash"]
