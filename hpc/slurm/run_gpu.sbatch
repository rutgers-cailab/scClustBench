#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --requeue
#SBATCH --output=logs/%x_%j.log
# ============================================================
# hpc/slurm/run_gpu.sbatch
#
# Runs one R script from pipeline/ that trains scVI, on a GPU node. R comes
# from the sccb environment and scvi-tools from the scvi Python environment
# through reticulate. The job stops at once (exit 3) if torch cannot see a
# GPU, so a misconfigured node does not silently train on CPU. The scvi
# environment was used with torch 2.11.0+cu128 on nodes with NVIDIA driver 570.
#
# Run:      SCVI_GPU=1 bash hpc/slurm/submit_all.sh, which calls
#           sbatch --job-name=13s_rods --cpus-per-task=8 --mem=128G --time=24:00:00 \
#             hpc/slurm/run_gpu.sbatch 13s_scvi_family.R POP=Rods
#           (submit from the repository root)
# Options:  SCVI_PYTHON  python with scvi-tools (default ~/micromamba/envs/scvi/bin/python)
#           SCCB_ENV     R environment prefix (default ~/micromamba/envs/sccb)
# ============================================================
set -euo pipefail
PY=${SCVI_PYTHON:-$HOME/micromamba/envs/scvi/bin/python}
# Fail fast if no CUDA device is visible to torch
$PY -c "import torch,sys; ok=torch.cuda.is_available(); print('torch', torch.__version__, 'cuda ok:', ok); sys.exit(0 if ok else 3)"
cd "$SLURM_SUBMIT_DIR"   # submit_all.sh submits from the repository root
# First argument is the script; every further VAR=value argument becomes an environment variable
SCRIPT=$1; shift; for kv in "$@"; do export "$kv"; done
export RETICULATE_PYTHON=$PY OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK PATH=${SCCB_ENV:-$HOME/micromamba/envs/sccb}/bin:$PATH
echo "[$(date)] $SCRIPT on $(hostname), env: $*"
Rscript "pipeline/$SCRIPT"
echo "[$(date)] done $SCRIPT"
