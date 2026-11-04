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

# Build lofreq from source using the morinlab fork of CSB5/lofreq.
# This fork adds several fixes on top of v2.1.5 that were never released:
#   - fixed --continue mode in lofreq somatic: pre-existing normal VCFs in
#     the output directory are detected correctly without parsing log files
#   - support for directory path with trailing slash as -o argument
#   - bcftools-based VCF merging in lofreq call-parallel (fixes failures with
#     multi-allelic variants) and masking of non-ACTG IUPAC reference bases
#   - fix AF computation
ARG LOFREQ_COMMIT=2df3fca
RUN git clone https://github.com/morinlab/lofreq.git /opt/lofreq_src \
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
