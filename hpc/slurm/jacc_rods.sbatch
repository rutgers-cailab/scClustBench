#!/bin/bash
#SBATCH --job-name=12c_rods
#SBATCH --requeue
#SBATCH --cpus-per-task=10
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.log
# ============================================================
# hpc/slurm/jacc_rods.sbatch
#
# Array job for the Stage A bootstrap Jaccard of the rod partitions. Task N
# reads line N of hpc/slurm/12c_rods_tasks.txt (one partition column, e.g.
# PCA_KMeans_k2) and runs its 20 bootstrap re-clusterings on 10 cores through
# 12c_recompute_jaccard.R. A separate merge job (12c_recompute_jaccard.R with
# MERGE=1) collects the results once every task has ended.
#
# Run:      bash hpc/slurm/submit_all.sh, which calls
#           sbatch --array=1-56 hpc/slurm/jacc_rods.sbatch   (after step 04; from the repository root)
# Options:  SCCB_ENV  R environment prefix (default ~/micromamba/envs/sccb)
# ============================================================
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
COL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" hpc/slurm/12c_rods_tasks.txt)
# Single-threaded BLAS/OpenMP: the R code forks one worker per allocated core, so threads would oversubscribe
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PATH=${SCCB_ENV:-$HOME/micromamba/envs/sccb}/bin:$PATH
echo "[$(date)] rods $COL on $(hostname)"
DATASET=rods ONLY=$COL N_BOOT=20 Rscript pipeline/12c_recompute_jaccard.R
