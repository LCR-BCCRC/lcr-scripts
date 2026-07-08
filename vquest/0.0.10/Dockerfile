FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"
RUN mamba create --yes --name env \
        --channel conda-forge \
        --channel defaults \
        "python=3.12" \
        "numpy" \
        "requests" \
        "pyyaml" \
        "pip" \
    && mamba clean --all --yes \
    && /opt/conda/envs/env/bin/pip install \
        "git+https://github.com/shawhahnlab/vquest.git@0c01f498ba01f4f289e58155c198f16dc7a277ac" \
        "lxml_html_clean==0.4.5"
ENV PATH=/opt/conda/envs/env/bin:$PATH
CMD ["/bin/bash"]
