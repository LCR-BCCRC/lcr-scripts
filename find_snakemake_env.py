#!/usr/bin/env python3
"""
find_snakemake_env.py
=====================
Locate the conda environment directory that Snakemake built from a given
environment YAML file.

Snakemake names conda environment directories using an MD5 hash of:
    1. The absolute real path of the conda environments directory
    2. The raw byte content of the environment YAML file

Usage
-----
    python find_snakemake_env.py <env_yaml> <conda_prefix>

Arguments
---------
env_yaml        Path to the conda environment YAML file.
conda_prefix    Directory where Snakemake stores conda environments
                (passed to Snakemake via --conda-prefix).

Example
-------
    python find_snakemake_env.py \\
        ~/lcr-scripts/envs/augment_manta_vcf/augment_manta_vcf-1.yaml \\
        /projects/rmorin_scratch/conda_environments

Notes
-----
Based on Snakemake's Env.hash property in snakemake/deployment/conda.py.
Verified against Snakemake 5.32.0; the hashing logic has been stable
across versions but may change in future releases.
"""

import sys
import os
import hashlib


def compute_hash(yaml_path, conda_prefix):
    env_dir = os.path.realpath(conda_prefix)
    with open(yaml_path, "rb") as f:
        content = f.read()
    md5 = hashlib.md5()
    md5.update(env_dir.encode())
    md5.update(content)
    return md5.hexdigest()


def find_env(yaml_path, conda_prefix):
    h = compute_hash(yaml_path, conda_prefix)

    # Snakemake appends an underscore to the hash
    full_path = os.path.join(conda_prefix, h + "_")
    short_path = os.path.join(conda_prefix, h[:8] + "_")

    if os.path.isdir(full_path):
        return full_path
    if os.path.isdir(short_path):
        return short_path
    return None


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <env_yaml> <conda_prefix>")
        sys.exit(1)

    yaml_path = sys.argv[1]
    conda_prefix = sys.argv[2]

    env_path = find_env(yaml_path, conda_prefix)
    if env_path:
        print(env_path)
    else:
        h = compute_hash(yaml_path, conda_prefix)
        print(
            f"Environment not found in {conda_prefix}\n"
            f"Expected directory: {h}_\n"
            f"The environment may not have been built yet, or was built with a "
            f"different --conda-prefix.",
            file=sys.stderr,
        )
        sys.exit(1)
