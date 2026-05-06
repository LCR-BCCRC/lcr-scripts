# lcr-scripts

Collection of curated scripts from the Morin Lab, used primarily within [lcr-modules](https://github.com/LCR-BCCRC/lcr-modules) Snakemake pipelines.

## Repository structure

Each script lives under `<script_name>/<version>/` and typically contains:

- The script itself (`.py`, `.sh`, etc.)
- A conda environment YAML (`<script_name>-<major>.yaml`)
- Optionally, a `Dockerfile` for containerized use

## Docker containers

A GitHub Actions workflow (`.github/workflows/build-containers.yml`) automatically builds and pushes a container to the GitHub Container Registry (GHCR) whenever a `Dockerfile` is merged into `master`. Containers are tagged using the script name and version derived from the directory path:

```
ghcr.io/lcr-bccrc/lcr-scripts/<script_name>:<version>
```

For example, `augment_manta_vcf/1.1/Dockerfile` produces:

```
ghcr.io/lcr-bccrc/lcr-scripts/augment_manta_vcf:1.1
```

These images are referenced in lcr-modules via `container_envs` in each module's `config/default.yaml` and used with Snakemake's `--use-apptainer` flag.

### Adding a container for a new script

1. Create a `Dockerfile` in the appropriate `<script_name>/<version>/` directory.
2. The Dockerfile should copy the conda environment YAML and install it into the base environment, e.g.:

```dockerfile
FROM condaforge/mambaforge:latest
LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

COPY <script_name>-<major>.yaml /tmp/env.yaml
RUN mamba env update --name base --file /tmp/env.yaml \
    && mamba clean --all --yes

CMD ["/bin/bash"]
```

3. Open a PR and merge to `master` — the CI workflow will build and push the image automatically.
