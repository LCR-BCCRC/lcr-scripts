FROM quay.io/biocontainers/ensembl-vep:100.4--pl526hecda079_0
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

# VEP 100.4 + Perl 5.26 + samtools are pre-installed in this biocontainer.
# ensembl-vep=100 cannot be installed via modern conda (broken Perl dep tree),
# so we use this pre-built biocontainer as base and download vcf2maf scripts
# directly from the v1.6.18 GitHub tag.
RUN curl -fsSL \
        https://raw.githubusercontent.com/mskcc/vcf2maf/v1.6.18/vcf2maf.pl \
        -o /usr/local/bin/vcf2maf.pl \
    && curl -fsSL \
        https://raw.githubusercontent.com/mskcc/vcf2maf/v1.6.18/maf2maf.pl \
        -o /usr/local/bin/maf2maf.pl \
    && chmod +x /usr/local/bin/vcf2maf.pl /usr/local/bin/maf2maf.pl

CMD ["/bin/bash"]
