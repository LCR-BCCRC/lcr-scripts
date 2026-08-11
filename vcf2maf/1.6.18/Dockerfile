FROM quay.io/biocontainers/ensembl-vep:100.4--pl526hecda079_0
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

# ensembl-vep=100 cannot be installed via modern conda: its bioconda recipe
# has a broken Perl dep tree (perl=5.26 required by ensembl-vep conflicts
# with perl>=5.32 required by perl-bioperl), independent of Python version.
# Using this pre-built biocontainer (VEP 100.4, Perl 5.26) avoids the
# conflict entirely.  vcf2maf is installed without its conda deps since VEP
# and the required Perl modules are already present; samtools is added
# separately.
RUN conda install --yes --no-deps \
        --channel bioconda \
        "vcf2maf=1.6.18" \
    && conda install --yes \
        --channel bioconda \
        --channel conda-forge \
        "samtools" \
    && conda clean -afy

CMD ["/bin/bash"]
