# ============================================================
# tests/test_speedups.R
#
# Checks that the faster implementations in R/scClustBench.R (the cached Louvain
# resolution path, parallel bootstraps, parallel per-embedding clustering and
# parallel benchmark metrics) return exactly the same results as the serial
# reference implementation: R/scClustBench.R as of git commit a0163fe, stored
# in tests/reference_scClustBench_a0163fe.R.
#
# Run:      Rscript tests/test_speedups.R
# ============================================================
suppressPackageStartupMessages({ library(Seurat) })
source("config.R"); source("R/utils.R")
# Load the reference implementation into the environment `old`
old <- new.env(); sys.source("tests/reference_scClustBench_a0163fe.R", envir = old)
source("R/scClustBench.R")
# ---- Test data: three Gaussian blobs in 10 dimensions ----
set.seed(3); n <- 1500
lat <- rbind(matrix(rnorm(n * 10 / 3), ncol = 10), matrix(rnorm(n * 10 / 3, 2), ncol = 10), matrix(rnorm(n * 10 / 3, -2), ncol = 10))
rownames(lat) <- paste0("c", seq_len(nrow(lat)))
# ---- Louvain path cache gives the same labels as a separate search for each k ----
p <- sccb_louvain_path(lat, SEED)
for (k in c(2, 3, 5, 8, 4)) stopifnot(identical(p(k), old$sccb_louvain_to_k(lat, k, SEED)), identical(sccb_louvain_to_k(lat, k, SEED), p(k)))
ref <- old$sccb_louvain_to_k(lat, 3, SEED)
# ---- Parallel bootstraps and parallel R(k) give the same values as the serial reference ----
Sys.setenv(SLURM_CPUS_PER_TASK = "4")
for (m in c("louvain", "kmeans")) stopifnot(all.equal(.jacc_min(SCCB_METHODS[[m]], lat, ref, 3, 6, SEED), old$.jacc_min(old$SCCB_METHODS[[m]], lat, ref, 3, 6, SEED)))
o <- CreateSeuratObject(.dummy_counts(rownames(lat)))
lats <- list(a = lat, b = lat + rnorm(length(lat), sd = .3), c = lat + rnorm(length(lat), sd = .6))
r_new <- sccb_reproducibility(o, latents = lats, k_range = 2:6); r_old <- old$sccb_reproducibility(o, latents = lats, k_range = 2:6)
stopifnot(isTRUE(all.equal(r_new, r_old)))
cat("speed-up equivalence: OK\n")
# ---- Benchmark: parallel metrics and shared paths give the same table as the serial reference ----
Sys.setenv(SCCB_BENCH_WORKERS = "3")
o2 <- CreateSeuratObject(matrix(rpois(50 * nrow(lat), 1), 50, nrow(lat), dimnames = list(paste0("g", 1:50), rownames(lat))))
o2 <- NormalizeData(o2, verbose = FALSE)
b_new <- sccb_benchmark(o2, "pca", c("louvain", "kmeans"), 2:4, n_boot = 2, rogue = FALSE, latents = list(pca = lat))
b_old <- old$sccb_benchmark(o2, "pca", c("louvain", "kmeans"), 2:4, n_boot = 2, rogue = FALSE, latents = list(pca = lat))
stopifnot(isTRUE(all.equal(b_new, b_old, check.attributes = FALSE)), identical(attr(b_new, "labels"), attr(b_old, "labels")))
cat("benchmark equivalence: OK\n")
