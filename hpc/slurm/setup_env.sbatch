#!/bin/bash
#SBATCH --job-name=sccb_env
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=%x_%j.log
# ============================================================
# hpc/slurm/setup_env.sbatch
#
# Builds the R environment "sccb" (R 4.5 and every scClustBench dependency) in
# a user micromamba installation, so no administrator rights are needed.
# Installs micromamba into ~/bin if it is missing, creates the environment from
# conda-forge and bioconda, then installs torch, scDHA and ROGUE (pinned
# GitHub commits) from inside R. It finishes by checking that every package
# loads, printing ENV_OK, and saving sessionInfo() to ~/sccb_env_sessionInfo.txt.
#
# Run:      sbatch hpc/slurm/setup_env.sbatch   (once, before submit_all.sh; log in the submit folder)
# ============================================================
set -euo pipefail
export MAMBA_ROOT_PREFIX=$HOME/micromamba
mkdir -p $HOME/bin
# ---- micromamba (single static binary) ----
if [ ! -x $HOME/bin/micromamba ]; then
  curl -Ls -o $HOME/bin/micromamba https://github.com/mamba-org/micromamba-releases/releases/download/2.9.0-0/micromamba-linux-64 && chmod +x $HOME/bin/micromamba
fi
MM=$HOME/bin/micromamba
# ---- Conda packages: R, compilers and system libraries, CRAN and Bioconductor packages ----
$MM create -y -n sccb -c conda-forge -c bioconda --channel-priority strict \
  r-base=4.5 matplotlib-base compilers make cmake pkg-config curl openssl libxml2 hdf5 git \
  r-seurat r-seuratobject r-harmony r-mclust r-aricode r-igraph r-irlba r-fnn r-uwot \
  r-dplyr r-tidyr r-ggplot2 r-patchwork r-matrix r-cluster r-remotes r-biocmanager \
  r-rhpcblasctl r-reticulate r-tibble r-rcpp r-rcpparmadillo \
  bioconductor-singlecellexperiment bioconductor-scran bioconductor-scuttle \
  bioconductor-splatter bioconductor-scdblfinder bioconductor-scrnaseq
# ---- Packages installed from R (torch backend for scDHA, GitHub packages), then a load check ----
$MM run -n sccb Rscript -e '
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8)
if (!requireNamespace("torch", quietly = TRUE)) install.packages("torch")
Sys.setenv(TORCH_INSTALL = "1"); torch::install_torch()
remotes::install_github("duct317/scDHA@ac9060b70768", upgrade = "never")
remotes::install_github("PaulingLiu/ROGUE@6e1c8f9ee419", upgrade = "never")
pk <- c("Seurat","SeuratObject","harmony","mclust","aricode","igraph","torch","scDHA","ROGUE","splatter","scran","scRNAseq","scDblFinder","reticulate")
v <- sapply(pk, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING"))
print(v); writeLines(capture.output(sessionInfo()), file.path(Sys.getenv("HOME"), "sccb_env_sessionInfo.txt"))
stopifnot(all(v != "MISSING")); cat("ENV_OK\n")'
