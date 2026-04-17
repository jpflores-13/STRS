#!/bin/bash

# setup_cooltools_env.sh
# One-time setup: creates a conda environment with hic2cool and cooltools.
# Run interactively on a login node — no SLURM needed.
# Usage: bash setup_cooltools_env.sh

set -euo pipefail

ENV_NAME="cooltools_env"

echo "Loading anaconda module..."
module load anaconda/2024.02

echo "Creating conda environment: ${ENV_NAME}"
conda create -y -n "${ENV_NAME}" python=3.10

echo "Activating environment..."
source activate "${ENV_NAME}"

echo "Installing hic2cool and cooltools via pip..."
pip install hic2cool cooltools

echo "Verifying installs..."
hic2cool --version
cooltools --version

echo ""
echo "Setup complete. To activate in future sessions:"
echo "  module load anaconda/2024.02 && source activate ${ENV_NAME}"