#!/bin/bash
#SBATCH --requeue
#SBATCH --output=logs/%x_%j.log
# ============================================================
# hpc/slurm/run.sbatch
#
# Generic CPU job: runs one script from pipeline/ in the R environment, with
# any VAR=value arguments exported to it first. The special name FIGURES runs
# every paper figure step of run_pipeline.R instead. OpenMP and BLAS threads
# are set to the number of allocated CPUs. Job name, CPUs, memory and time
# limit are given on the sbatch command line by submit_all.sh.
#
# Run:      sbatch --job-name=04_rods --cpus-per-task=32 --mem=180G --time=36:00:00 \
#             hpc/slurm/run.sbatch 04_phase2_rod_benchmark.R BENCH_JACC=0
#           (submit_all.sh does this for every CPU step; submit from the repository root)
# Options:  SCCB_ENV  R environment prefix (default ~/micromamba/envs/sccb)
# ============================================================
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"   # submit_all.sh submits from the repository root
# First argument is the script; every further VAR=value argument becomes an environment variable
SCRIPT=$1; shift; for kv in "$@"; do export "$kv"; done
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK PATH=${SCCB_ENV:-$HOME/micromamba/envs/sccb}/bin:$PATH
echo "[$(date)] $SCRIPT on $(hostname), $SLURM_CPUS_PER_TASK cpus, env: $*"
if [ "$SCRIPT" = FIGURES ]; then Rscript run_pipeline.R --from fig1   # every paper figure
else Rscript "pipeline/$SCRIPT"; fi
echo "[$(date)] done $SCRIPT"
