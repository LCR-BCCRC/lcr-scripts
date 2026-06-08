FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

RUN mamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        "vcf2maf=1.6.18" \
        "samtools" \
        "perl-bioperl" \
        "perl-set-intervaltree" \
        "perl-perlio-gzip" \
        "perl-archive-zip" \
    && mamba clean -afy

# VEP 100.2 from source — bypasses bioconda's Perl 5.26 packaging conflict
RUN git clone --branch release/100 --depth 1 \
        https://github.com/Ensembl/ensembl-vep.git /opt/vep \
    && ln -s /opt/vep/vep /opt/conda/bin/vep

CMD ["/bin/bash"]
