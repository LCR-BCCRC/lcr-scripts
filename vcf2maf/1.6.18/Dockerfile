FROM quay.io/biocontainers/ensembl-vep:100.4--pl526hecda079_0
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

# VEP 100.4 + Perl 5.26 are pre-installed in this biocontainer.
# The biocontainer is a minimal Debian image with no conda in PATH, so
# vcf2maf is installed directly from source and samtools via apt-get.
RUN curl -fsSL \
        https://raw.githubusercontent.com/mskcc/vcf2maf/v1.6.18/vcf2maf.pl \
        -o /usr/local/bin/vcf2maf.pl \
    && curl -fsSL \
        https://raw.githubusercontent.com/mskcc/vcf2maf/v1.6.18/maf2maf.pl \
        -o /usr/local/bin/maf2maf.pl \
    && chmod +x /usr/local/bin/vcf2maf.pl /usr/local/bin/maf2maf.pl \
    && apt-get update -qq \
    && apt-get install -y --no-install-recommends samtools \
    && rm -rf /var/lib/apt/lists/*

CMD ["/bin/bash"]
