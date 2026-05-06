FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

COPY augment_manta_vcf-1.yaml /tmp/env.yaml
RUN mamba env update --name base --file /tmp/env.yaml \
    && mamba clean --all --yes

CMD ["/bin/bash"]
