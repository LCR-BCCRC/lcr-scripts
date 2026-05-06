# lcr-scripts

Collection of curated scripts from the Morin Lab, used primarily within [lcr-modules](https://github.com/LCR-BCCRC/lcr-modules) Snakemake pipelines.

## Repository structure

Each script lives under `<script_name>/<version>/` and typically contains:

- The script itself (`.py`, `.sh`, etc.)
- A conda environment YAML (`<script_name>-<major>.yaml`)
- Optionally, a `Dockerfile` for containerized use
- Optionally, a `run_tests.sh` and `tests/` directory for regression tests

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

RUN mamba install --yes --name base \
        --channel bioconda \
        --channel conda-forge \
        <package1> \
        <package2> \
    && mamba clean --all --yes

CMD ["/bin/bash"]
```

Install only the packages the script actually imports — avoid reusing the conda YAML directly, as
it is a fully-pinned export that may be incompatible with the Python version in the base image.

3. Open a PR and merge to `master` — the CI workflow will build and push the image automatically.

> **Note:** The conda YAML files in each script directory are symlinks into `envs/`. Docker cannot
> follow symlinks that escape the build context, so Dockerfiles must use the real path relative to
> the repo root (e.g. `COPY envs/<script_name>/<script_name>-<major>.yaml /tmp/env.yaml`) and the
> CI workflow sets the build context to `.` (repo root). Do not use the local symlink path.

## Testing

Scripts can include a `run_tests.sh` that runs the script on small curated inputs and writes outputs
to `tests/output/`. This follows a **golden-file** pattern:

1. **Establish a baseline (once):** Run `./run_tests.sh` from the script's version directory, then
   commit the files written to `tests/output/`. These become the expected outputs.

2. **Detect regressions (every subsequent run):** Run `./run_tests.sh` again — it overwrites
   `tests/output/`. Then check for unexpected changes with:

   ```bash
   git diff tests/output/
   ```

   No diff means the outputs are identical to the baseline. Any diff indicates a regression.

3. **Environment-variable header lines:** Some output files (e.g. VCFs) include header lines with
   absolute paths (`##cmdline`, `##regions_bed`) that differ between environments and will always
   appear in `git diff`. To compare only the meaningful content, filter them out:

   ```bash
   grep -v "^##cmdline\|^##regions_bed" tests/output/file.vcf | diff - tests/output/file.vcf
   ```

### CI integration

When a `Dockerfile` is present alongside a `run_tests.sh`, the CI workflow automatically:

1. Builds the Docker image from the `Dockerfile`.
2. Runs `run_tests.sh` **inside the container** — this validates that the container environment
   actually works, not just that the script runs locally.
3. Uploads the generated `tests/output/` files as workflow artifacts.

To lock in the baseline after the first successful CI run, download the artifacts, place the files
in `tests/output/`, and commit them.

### Adding tests for a new script

1. Create `tests/input/` with small, self-contained input files.
2. Create an empty `tests/output/` directory (add a `.gitkeep` so git tracks it).
3. Write `run_tests.sh` — start with `set -euo pipefail` so any failure exits non-zero.
   Put positional arguments **before** any option that uses `nargs="+"` (e.g. `--bed_regions`)
   to avoid argparse consuming them greedily.
4. Run `./run_tests.sh` once, inspect the outputs, and commit `tests/output/`.
