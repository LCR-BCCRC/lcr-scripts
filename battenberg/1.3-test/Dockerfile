FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

# TEST IMAGE: battenberg:1.2 plus the experimental recentre_logr option (morinlab/battenberg branch
# recentre-logr). Opt-in only: point a project config's container_envs: battenberg at this tag.
# battenberg:1.2 is unchanged. If the patch changes, add a new tag (e.g. 1.3-test2) rather than
# rebuilding this one, so cached .simg copies never go stale.

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
    && R --vanilla -q -e 'remotes::install_github("morinlab/battenberg@3a8d3b9e8f14a0a83790e715e8c3db3a5fee69e2", lib="/opt/conda/lib/R/library", dependencies = TRUE, upgrade = "never")'

CMD ["/bin/bash"]
