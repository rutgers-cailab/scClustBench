# ============================================================
# examples/run_on_own_data.R
#
# Runs scClustBench on your own Seurat object and writes the benchmark,
# Stage A selection, R(k) curve, verdict and reported partition as CSVs.
# The input should hold raw counts in the RNA assay and be restricted to the
# population you want to test: a whole dataset, or one cell type when asking
# whether that type has reproducible subtypes.
#
# Inputs:   <cells.rds>, a Seurat object, and the name of a metadata column holding
#           the batch (sample or donor), or "none"
# Outputs:  <out_dir>/sccb_{scores,selection,curve,verdict,recommended,partition}.csv
# Run:      Rscript examples/run_on_own_data.R my_cells.rds sample sccb_out [k_max] [n_boot]
# Options:  k_max   largest k in the sweep k = 2..k_max (default 10)
#           n_boot  bootstrap resamples behind the Jaccard gate (default 20, as in the paper;
#                   5 gives a quick look)
#
# Run from the repository root (paths below are relative to it).
# ============================================================
source("R/scClustBench.R")
args   <- commandArgs(trailingOnly = TRUE)
path   <- if (length(args) >= 1) args[1] else "my_cells.rds"
batch  <- if (length(args) >= 2 && args[2] != "none") args[2] else NULL
out    <- if (length(args) >= 3) args[3] else "sccb_out"
k_max  <- if (length(args) >= 4) as.integer(args[4]) else 10L
n_boot <- if (length(args) >= 5) as.integer(args[5]) else 20L
if (!file.exists(path)) stop("no input '", path, "'. Usage: Rscript examples/run_on_own_data.R ",
                             "<cells.rds> <batch column|none> <out_dir> [k_max] [n_boot]. ",
                             "For a demo on simulated data run: Rscript tests/smoke_test.R", call. = FALSE)

# ---- Load the population ----
obj <- readRDS(path)
if (!is.null(batch) && !batch %in% colnames(obj@meta.data)) stop("batch column '", batch, "' not in metadata", call. = FALSE)

# ---- Optional: register extra methods and embeddings ----
# Anything with these signatures is benchmarked and scored like the built-ins.
SCCB_METHODS$hclust_ward <- function(latent, k, seed) cutree(hclust(dist(latent), "ward.D2"), k)  # O(n^2): small data only
SCCB_EMBEDDINGS$pca15 <- function(obj, seed, batch = NULL, ...) SCCB_EMBEDDINGS$pca(obj, seed, batch, npcs = 15)

# ---- Run Stage A and Stage B ----
res <- sccb_run(obj, batch = batch, k_range = 2:k_max,
                embeddings = c("pca", "pca15"),                 # what to benchmark (Stage A input)
                methods = c("kmeans", "gmm", "louvain", "leiden", if (ncol(obj) <= 5000) "hclust_ward"),
                # Stage B: independent embeddings for R(k). sccb_run() defaults to scDHA if installed,
                # else pca_perturb; listing several checks whether the verdict agrees across families.
                repro_embeddings = c(if (requireNamespace("scDHA", quietly = TRUE)) "scdha", "pca_perturb"),
                seeds = c(1, 2, 3, 7, 42), tau = 0.8, n_boot = n_boot, out_dir = out)

# ---- Results ----
# res$scores       every (embedding_method, k): balance, jacc_min, marker_spec, rogue, silhouette
# res$selection    Stage A: admissible partitions ranked by mean metric percentile
# res$curve        Stage B: R(k) = mean pairwise ARI across independent embeddings (+ smallest cluster)
# res$verdict      per repro embedding: ACCEPT (sub-structure reproducible; k* = reproducible granularity)
#                  or REJECT
# res$recommended  the answer: verdict, k*, and the configuration chosen at k* (top admissible at k*)
# res$partition    cell -> cluster for that configuration (also out_dir/sccb_partition.csv)
print(res$verdict[, c("embedding", "verdict", "k_star", "k_star_smallest", "k_star_at_edge", "max_R")])
if (!is.null(res$partition)) print(table(res$partition$cluster))
