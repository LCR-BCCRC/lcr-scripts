FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

# Build tools needed to compile lofreq from source
RUN apt-get update -q && apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        build-essential \
        git \
        libtool \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# htslib is required at build time; bcftools is required at runtime by
# lofreq2_call_pparallel.py (the VCF-merging fix introduced in PR #109
# switched from lofreq vcfset to bcftools concat)
RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        "htslib" \
        "bcftools" \
        "samtools" \
    && mamba clean -afy

# Build lofreq from source at a specific post-v2.1.5 master commit.
# This picks up patches that were never included in any tagged release:
#   - PR #109: bcftools-based VCF merging in parallel mode (fixes failures
#              with multi-allelic variants), and masking of non-ACTG IUPAC
#              reference bases (e.g. 'W' at chr17:83129591 in hg38) as 'N'
#   - PR #126: fix AF computation
#   - PR #130: fix GSL include error
#   - PR #155: handle large BAM headers
ARG LOFREQ_COMMIT=46bd58d
RUN git clone https://github.com/CSB5/lofreq.git /opt/lofreq_src \
    && cd /opt/lofreq_src \
    && git checkout ${LOFREQ_COMMIT} \
    && ./bootstrap \
    && CPPFLAGS="-I/opt/conda/include" \
       LDFLAGS="-L/opt/conda/lib -Wl,-rpath,/opt/conda/lib" \
       PKG_CONFIG_PATH="/opt/conda/lib/pkgconfig" \
       ./configure --prefix=/usr/local --with-htslib=/opt/conda \
    && make \
    && make install \
    && cd / \
    && rm -rf /opt/lofreq_src

CMD ["/bin/bash"]
