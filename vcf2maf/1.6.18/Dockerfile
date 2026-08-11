FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

# Remove the python 3.12.* pin so mamba can resolve a Python version
# compatible with ensembl-vep=100 (which chains through t_coffee <=3.7).
# vcf2maf, ensembl-vep, and samtools are all Perl/C so Python version
# inside the container does not affect their runtime behaviour.
RUN rm -f /opt/conda/conda-meta/pinned \
    && mamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        "vcf2maf=1.6.18" \
        "ensembl-vep=100" \
        "samtools" \
    && mamba clean -afy

CMD ["/bin/bash"]
