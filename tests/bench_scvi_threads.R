# ============================================================
# tests/bench_scvi_threads.R
#
# Measures seconds per scVI training epoch on CPU for the rod population at
# 1, 4, 8 and 16 torch threads (2 epochs each), to choose the thread count for
# the scVI embedding family (step 13s).
#
# Inputs:   output/rds/07_rods_positive_control.rds (from Rscript run_pipeline.R 07_rods)
# Run:      RETICULATE_PYTHON=<python with scvi-tools> Rscript tests/bench_scvi_threads.R
# Options:  RETICULATE_PYTHON  python with scvi-tools and torch (required)
# ============================================================
source("config.R"); source("R/scClustBench.R")
# The 07 output stores the population object in $bc whatever the cell type (here the rods)
pc  <- readRDS(file.path(RDS_DIR, "07_rods_positive_control.rds")); obj <- sccb_prep(pc$bc, "LogNormalize", 3000)
torch <- reticulate::import("torch")   # same python session that scVI uses
for (th in c(1, 4, 8, 16)) {
  torch$set_num_threads(as.integer(th))
  t0 <- Sys.time(); lat <- sccb_embed(obj, "scvi", 1, batch = "sample", max_epochs = 2)
  if (is.null(lat)) stop("scVI training failed; check RETICULATE_PYTHON and the scvi-tools install")
  cat(sprintf("threads %2d: %.1f s per epoch\n", th, as.numeric(difftime(Sys.time(), t0, units = "secs")) / 2))
}
