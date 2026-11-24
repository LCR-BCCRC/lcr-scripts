# battenberg-base

Battenberg's dependencies (R, ASCAT, the CRAN packages Battenberg imports, impute2, alleleCounter)
without Battenberg itself. No module uses this image. It exists so that Battenberg test images only
install the Battenberg commit under test: they build in about a minute, and differ from each other
and from `battenberg:1.2` only in that commit.

`battenberg-base:1.0` pins everything to what `battenberg:1.2` and `battenberg:1.3-test` were built
with: the mambaforge image digest, R 4.5.3 and the other conda packages, ASCAT commit `e3d6577`, and CRAN
packages from the 2026-11-21 Posit Package Manager snapshot. To change any of them, add a new version
directory (`battenberg-base/1.1/`); never rebuild an existing tag.

## Making a Battenberg test image

1. Push the Battenberg commit to `morinlab/battenberg` (any branch).
2. Add `battenberg/<new tag>/Dockerfile`, e.g. `battenberg/1.3-test2/Dockerfile`, with:

   ```dockerfile
   FROM ghcr.io/lcr-bccrc/lcr-scripts/battenberg-base:1.0
   LABEL org.opencontainers.image.source="https://github.com/LCR-BCCRC/lcr-scripts"

   # dependencies = FALSE: everything Battenberg imports is in the base, and installing fails loudly
   # if a new commit adds an import the base lacks (then make battenberg-base/1.1).
   RUN --mount=type=secret,id=GITHUB_PAT \
       if [ -s /run/secrets/GITHUB_PAT ]; then export GITHUB_PAT="$(cat /run/secrets/GITHUB_PAT)"; fi \
       && R --vanilla -q -e 'remotes::install_github("morinlab/battenberg@<full commit sha>", lib = "/opt/conda/lib/R/library", dependencies = FALSE, upgrade = "never")'

   CMD ["/bin/bash"]
   ```

   Copy `battenberg/1.3-test/run_tests.sh` next to it, adding a check for any new argument.
3. Open a PR (CI builds and tests the image), then merge it to publish
   `ghcr.io/lcr-bccrc/lcr-scripts/battenberg:<new tag>`.
4. Point a project config at it with `container_envs: battenberg: "docker://ghcr.io/lcr-bccrc/lcr-scripts/battenberg:<new tag>"`.

Use a new tag for every Battenberg commit rather than rebuilding one: Snakemake caches each image by
its URL, so a rebuilt tag leaves stale cached copies behind.

## Build order

A test image can only build once `battenberg-base:1.0` has been published, i.e. after the PR that adds
this directory is merged. Don't add a test image in the same PR as a new base version: CI builds
images in parallel, and PR builds never publish, so the test image would not find its base.
