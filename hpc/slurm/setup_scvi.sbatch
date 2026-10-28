#!/bin/bash
#SBATCH --job-name=scvi_env
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=%x_%j.log
# ============================================================
# hpc/slurm/setup_scvi.sbatch
#
# Builds the Python environment "scvi" (Python 3.11, scvi-tools 1.4.2 and the
# torch it depends on) that the scVI steps call through reticulate. Uses the
# micromamba installed by setup_env.sbatch, so run that first. Prints SCVI_OK
# when scvi and torch import.
#
# Run:      sbatch hpc/slurm/setup_scvi.sbatch   (once, after setup_env.sbatch; log in the submit folder)
# ============================================================
set -euo pipefail
MM=$HOME/bin/micromamba; export MAMBA_ROOT_PREFIX=$HOME/micromamba
$MM create -y -n scvi -c conda-forge python=3.11 pip
$MM run -n scvi pip install --no-cache-dir "scvi-tools==1.4.2"
$MM run -n scvi python -c "import scvi, torch; print('scvi', scvi.__version__, 'torch', torch.__version__); print('SCVI_OK')"
