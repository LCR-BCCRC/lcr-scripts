FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"
RUN mamba create --yes --name env \
        --channel conda-forge \
        --channel bioconda \
        --channel defaults \
        "blast=2.13.0" \
        "hisat2=2.2.1" \
        "kallisto=0.48.0" \
        "samtools=1.16.1" \
        "trinity=2.13.2" \
        "git" \
        "jq" \
    && mamba clean --all --yes \
    && /opt/conda/envs/env/bin/git clone --branch v1.0.1 --depth 1 \
        https://github.com/ForconiLab/IgSeqR.git /tmp/IgSeqR \
    && mkdir -p /opt/conda/envs/env/bin/data/igseqr \
    && cp /tmp/IgSeqR/bin/igseqr.sh /opt/conda/envs/env/bin/ \
    && chmod +x /opt/conda/envs/env/bin/igseqr.sh \
    && ln -sf /opt/conda/envs/env/bin/igseqr.sh /opt/conda/envs/env/bin/igseqr \
    && cp -r /tmp/IgSeqR/data/* /opt/conda/envs/env/bin/data/igseqr/ \
    && rm -rf /tmp/IgSeqR
ENV PATH=/opt/conda/envs/env/bin:$PATH
CMD ["/bin/bash"]
