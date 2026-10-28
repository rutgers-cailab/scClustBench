#!/bin/bash
#SBATCH --requeue
#SBATCH --output=logs/%x_%j.log
# ============================================================
# hpc/slurm/run_scvi_cpu.sbatch
#
# Runs one R script from pipeline/ that trains scVI, on CPU nodes (the default
# in submit_all.sh; useful when the GPU queue is long). It uses the same scvi
# Python environment as run_gpu.sbatch; torch runs on the allocated cores.
# R's OpenBLAS is kept single-threaded: a multi-threaded OpenBLAS in the same
# process competes with torch for the cores and made rod training about 13
# times slower (178 s versus 14 s per epoch).
#
# Run:      bash hpc/slurm/submit_all.sh, which calls
#           sbatch --job-name=13s_rods --cpus-per-task=8 --mem=128G --time=24:00:00 \
#             hpc/slurm/run_scvi_cpu.sbatch 13s_scvi_family.R POP=Rods
#           (submit from the repository root)
# Options:  SCVI_PYTHON  python with scvi-tools (default ~/micromamba/envs/scvi/bin/python)
#           SCCB_ENV     R environment prefix (default ~/micromamba/envs/sccb)
# ============================================================
set -euo pipefail
PY=${SCVI_PYTHON:-$HOME/micromamba/envs/scvi/bin/python}
cd "$SLURM_SUBMIT_DIR"   # submit_all.sh submits from the repository root
# First argument is the script; every further VAR=value argument becomes an environment variable
SCRIPT=$1; shift; for kv in "$@"; do export "$kv"; done
# torch (OpenMP/MKL) gets every allocated core; R's OpenBLAS gets one thread (see above)
export RETICULATE_PYTHON=$PY OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=1 PATH=${SCCB_ENV:-$HOME/micromamba/envs/sccb}/bin:$PATH
echo "[$(date)] $SCRIPT on $(hostname) (CPU scVI), env: $*"
Rscript "pipeline/$SCRIPT"
echo "[$(date)] done $SCRIPT"
