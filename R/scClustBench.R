# ============================================================
# R/scClustBench.R
#
# The scClustBench API: benchmark candidate partitions of a Seurat object,
# screen and rank them (Stage A), measure the reproducibility curve R(k) across
# independently built embeddings (Stage B), and return a verdict and, on ACCEPT,
# the reported partition. The pipeline scripts call these functions on the retina
# data; examples/run_on_own_data.R shows how to run them on another dataset.
#
# Usage (from the repository root):
#   source("R/scClustBench.R")
#   res <- sccb_run(seurat_obj, batch = "sample", k_range = 2:10, out_dir = "sccb_out")
#
# Main functions:
#   sccb_prep()            normalise, select HVGs and scale once; embeddings reuse it
#   sccb_embed()           build one embedding (cells x dims) from SCCB_EMBEDDINGS
#   sccb_benchmark()       embeddings x methods x k -> five label-free metrics per partition
#   sccb_select()          Stage A (R/decision_rule.R): admissible + ranked partitions
#   sccb_reproducibility() R(k) = mean pairwise ARI across independently built embeddings
#   sccb_verdict()         Stage B (R/decision_rule.R): ACCEPT/REJECT + k*
#   sccb_recommend()       the reported partition: top admissible Stage A configuration at k*
#   sccb_run()             all of the above in one call; optionally writes sccb_*.csv
#
# Extending (add entries by assignment after sourcing; no other code changes needed):
#   SCCB_METHODS$name    <- function(latent, k, seed) -> integer labels (length nrow(latent))
#   SCCB_EMBEDDINGS$name <- function(obj, seed, batch = NULL, ...) -> cells x dims matrix
#                           (obj has already been through sccb_prep(): normalised,
#                            VariableFeatures set, HVGs scaled)
#   SCCB_METHOD_PATHS$name <- function(latent, seed) -> function(k) labels, for a method whose
#                           work can be shared across k (see sccb_louvain_path)
#
# Required: Seurat, igraph, aricode, cluster. Optional (skipped with a message
# or NA metric if absent): mclust, harmony, ROGUE, scran, scDHA+torch, scvi-tools.
#
# Options:  SLURM_CPUS_PER_TASK  worker processes for parallel loops (default N_CORES from config.R, else 1)
#           SCCB_BENCH_WORKERS   cap on workers for sccb_benchmark() metrics (default 6)
#           RETICULATE_PYTHON    python with scvi-tools; required for the scvi embedding
# ============================================================
suppressPackageStartupMessages({ library(Seurat); library(igraph) })

# ---- Setup: load helpers and defaults ----
# Source R/utils.R and R/decision_rule.R relative to this file (found via the source() frame),
# unless the caller already loaded them
local({
  f <- NULL
  for (i in rev(seq_len(sys.nframe()))) { f <- sys.frame(i)$ofile; if (!is.null(f)) break }
  root <- if (is.null(f)) "." else file.path(dirname(f), "..")
  if (!exists("find_resolution_for_k")) sys.source(file.path(root, "R", "utils.R"), envir = globalenv())
  if (!exists("sccb_select")) sys.source(file.path(root, "R", "decision_rule.R"), envir = globalenv())
})
# Defaults normally supplied by config.R, set only if config.R was not sourced (utils.R reads them lazily)
for (.v in list(c("SEED", 1), c("MAX_SIL_CELLS", 10000), c("PUB_BASE_SIZE", 10), c("FIG_DPI", 600)))
  if (!exists(.v[1])) assign(.v[1], as.numeric(.v[2]))
rm(.v)

# ---- Internal helpers ----
# .need(pkg): stops with a clear message if an optional package is missing
.need <- function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) stop(pkg, " not installed", call. = FALSE)
# .ari(a, b), .nmi(a, b): ARI / NMI between two label vectors.
# aricode 1.1.0 errors or returns NA on character labels; factors work in every version.
.ari <- function(a, b) aricode::ARI(factor(a), factor(b))
.nmi <- function(a, b) aricode::NMI(factor(a), factor(b))
# .dummy_counts(ids): empty 2 x n counts matrix, so a Seurat object can carry an embedding alone
# (SeuratObject >= 5.3 rejects a 1-row matrix)
.dummy_counts <- function(ids) Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(2, length(ids)),
                                                    dimnames = list(c("d1", "d2"), ids))
# .dummy_obj(latent): Seurat object holding `latent` (cells x dims) as reduction "lat"
.dummy_obj <- function(latent) {
  colnames(latent) <- paste0("L_", seq_len(ncol(latent)))   # Seurat requires key-prefixed dim names
  o <- suppressWarnings(CreateSeuratObject(.dummy_counts(rownames(latent))))
  o[["lat"]] <- CreateDimReducObject(embeddings = latent, key = "L_", assay = "RNA"); o
}

# ---- Parallel execution ----
# .ncores(): number of worker processes (SLURM_CPUS_PER_TASK, else N_CORES, else 1)
.ncores <- function() max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", if (exists("N_CORES")) N_CORES else 1L)))
# .par(X, FUN, cores): fork-parallel lapply (mclapply); returns the same list as lapply(X, FUN)
# because every task sets its own seed. Stops on the first worker error.
.par <- function(X, FUN, cores = .ncores()) {
  if (isTRUE(getOption("sccb.in_par"))) cores <- 1L      # no nested forking: inner calls run serially inside a worker
  f <- function(x) { options(sccb.in_par = TRUE); FUN(x) }
  out <- if (cores > 1 && length(X) > 1) parallel::mclapply(X, f, mc.cores = min(cores, length(X)), mc.preschedule = FALSE) else lapply(X, FUN)
  bad <- vapply(out, inherits, logical(1), "try-error"); if (any(bad)) stop(as.character(out[[which(bad)[1]]]))
  out
}
# ---- Clustering configurations ----
# sccb_louvain_path(latent, seed, ...): returns a function(k) giving Louvain labels with k clusters (or the
# closest k reached). The Seurat SNN graph is built once and each resolution of a fixed 60-point log grid is
# clustered at most once, so asking for several k reuses the work; each k gets the same labels as a separate
# per-k search (same graph, same grid order, set.seed(seed) before each FindClusters).
# `...` goes to FindNeighbors (e.g. annoy.metric = "cosine", k.param = 30).
sccb_louvain_path <- function(latent, seed = SEED, ...) {
  if (is.null(rownames(latent))) rownames(latent) <- paste0("c", seq_len(nrow(latent)))
  o <- FindNeighbors(.dummy_obj(latent), reduction = "lat", dims = 1:ncol(latent), verbose = FALSE, ...)
  grid <- exp(seq(log(1e-4), log(10), length.out = 60)); labs <- vector("list", length(grid)); nks <- rep(NA_integer_, length(grid))
  at <- function(i) {                                   # labels at grid[i], cached in labs/nks
    if (is.na(nks[i])) { set.seed(seed); lab <- as.integer(Idents(FindClusters(o, resolution = grid[i], algorithm = 1, verbose = FALSE)))
                         labs[[i]] <<- lab; nks[i] <<- length(unique(lab)) }
    labs[[i]]
  }
  # First grid point (lowest resolution) with exactly k clusters, else the closest k over the whole grid
  function(k) { bl <- NULL; bd <- Inf
    for (i in seq_along(grid)) { lab <- at(i); if (abs(nks[i] - k) < bd) { bd <- abs(nks[i] - k); bl <- lab }; if (bd == 0) break }
    bl }
}
# sccb_louvain_to_k(latent, k, seed, ...): Louvain labels for a single k; `...` goes to FindNeighbors
sccb_louvain_to_k <- function(latent, k, seed = SEED, ...) sccb_louvain_path(latent, seed, ...)(k)

# Methods whose labels for several k can share one cached resolution path (see sccb_louvain_path)
SCCB_METHOD_PATHS <- list(louvain = function(latent, seed) sccb_louvain_path(latent, seed))

# Built-in clustering methods: each takes (latent, k, seed) and returns integer labels.
# GMM uses at most the first 20 dimensions and returns a single cluster if Mclust fails.
# Leiden uses the CPM objective on the Seurat SNN graph (find_resolution_for_k in R/utils.R).
SCCB_METHODS <- list(
  kmeans  = function(latent, k, seed) { set.seed(seed); kmeans(latent, k, nstart = 25, iter.max = 100)$cluster },
  gmm     = function(latent, k, seed) {
    .need("mclust"); suppressPackageStartupMessages(library(mclust)); set.seed(seed)   # Mclust needs mclustBIC on the search path
    g <- tryCatch(mclust::Mclust(latent[, 1:min(20, ncol(latent))], G = k, verbose = FALSE), error = function(e) NULL)
    if (is.null(g)) rep(1L, nrow(latent)) else as.integer(g$classification)
  },
  louvain = function(latent, k, seed) sccb_louvain_to_k(latent, k, seed),
  leiden  = function(latent, k, seed) {
    if (is.null(rownames(latent))) rownames(latent) <- paste0("c", seq_len(nrow(latent)))
    colnames(latent) <- paste0("L_", seq_len(ncol(latent)))
    g <- build_seurat_snn_igraph(latent, rownames(latent), .dummy_counts(rownames(latent)),
                                 "lat", "L_", 1:ncol(latent))
    as.integer(find_resolution_for_k(g, k, algorithm = "leiden", res_min = 1e-6, res_max = 5,
                                     max_iter = 30, seed = seed)$labels)
  })

# ---- Preprocessing (done once; embeddings reuse it) ----
# sccb_prep(obj, norm, nfeatures): returns obj normalised (LogNormalize, SCTransform or scran size
# factors), with nfeatures variable genes selected and scaled. RNA layers are joined first.
sccb_prep <- function(obj, norm = c("LogNormalize", "SCT", "scran"), nfeatures = 3000) {
  norm <- match.arg(norm)
  DefaultAssay(obj) <- "RNA"
  if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  obj <- NormalizeData(obj, verbose = FALSE)            # RNA 'data' also used for marker_spec
  if (norm == "SCT") return(SCTransform(obj, variable.features.n = nfeatures, verbose = FALSE))
  if (norm == "scran") {
    .need("scran"); cts <- GetAssayData(obj, layer = "counts")
    sf <- scran::calculateSumFactors(cts, clusters = scran::quickCluster(cts))
    obj <- SetAssayData(obj, layer = "data", new.data = Matrix::Matrix(scuttle::normalizeCounts(cts, size.factors = sf), sparse = TRUE))
  }
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures, verbose = FALSE)
  ScaleData(obj, verbose = FALSE)
}

# ---- Embeddings ----
# .harmonize(lat, obj, batch, seed): Harmony-corrects `lat` on metadata column `batch` (unchanged if
# batch is NULL or harmony is missing). Parameters are the Harmony defaults, as in config.R.
.harmonize <- function(lat, obj, batch, seed) {
  if (is.null(batch)) return(lat)
  if (!requireNamespace("harmony", quietly = TRUE)) { message("  harmony not installed: no batch correction"); return(lat) }
  set.seed(seed)
  harmony::RunHarmony(lat, obj@meta.data[, batch, drop = FALSE], vars_use = batch,
                      theta = 2, lambda = 1, max_iter = 10, verbose = FALSE)
}
# .pca(obj, seed, batch, npcs, features): seeded PCA on `features` (default: the variable genes), then Harmony
.pca <- function(obj, seed, batch, npcs, features = NULL) {
  obj <- RunPCA(obj, features = features, npcs = npcs, seed.use = seed, verbose = FALSE)
  .harmonize(Embeddings(obj, "pca"), obj, batch, seed)
}

# Built-in embeddings: each takes (obj, seed, batch, ...) and returns a cells x dims matrix
SCCB_EMBEDDINGS <- list(
  pca = function(obj, seed, batch = NULL, npcs = 30, ...) .pca(obj, seed, batch, npcs),
  # model-agnostic embedding perturbation: each seed resamples `frac` of the HVGs + seeded randomized PCA
  pca_perturb = function(obj, seed, batch = NULL, npcs = 30, frac = 0.8, ...) {
    set.seed(seed); hv <- VariableFeatures(obj)
    .pca(obj, seed, batch, npcs, features = sample(hv, round(frac * length(hv))))
  },
  # scDHA autoencoder on log2 raw counts (as in step 07); Harmony is applied to its latent space
  scdha = function(obj, seed, batch = NULL, ngenes = 3000, ...) {
    .need("scDHA"); .need("torch")
    x <- log2(t(as.matrix(GetAssayData(obj, assay = "RNA", layer = "counts"))) + 1)
    set.seed(seed)
    sc <- scDHA::scDHA(data = x, k = NULL, method = "scDHA", sparse = FALSE, n = ngenes,
                       ncores = .ncores(), gen_fil = TRUE, do.clus = TRUE, seed = seed)
    lat <- as.matrix(sc$latent); rownames(lat) <- rownames(x)
    .harmonize(lat, obj, batch, seed)
  },
  # scVI via reticulate, on the variable genes; batch is modelled by scVI itself (no Harmony).
  # Opt-in via RETICULATE_PYTHON so reticulate never triggers a python auto-install.
  scvi = function(obj, seed, batch = NULL, n_latent = 10, max_epochs = 200, ...) {
    if (!nzchar(Sys.getenv("RETICULATE_PYTHON"))) stop("scvi: set RETICULATE_PYTHON to a python with scvi-tools", call. = FALSE)
    .need("reticulate")
    if (!reticulate::py_module_available("scvi")) stop("scvi-tools not importable", call. = FALSE)
    scvi <- reticulate::import("scvi"); ad <- reticulate::import("anndata")
    # scVI models integer UMI counts; two GEO samples (P15, GSM8133387/89) were deposited as
    # non-integer (ambient-corrected) matrices, so round to the nearest integer here.
    X <- round(Matrix::t(GetAssayData(obj, assay = "RNA", layer = "counts")[VariableFeatures(obj), ]))
    obs <- if (is.null(batch)) data.frame(row.names = colnames(obj)) else
      data.frame(b = as.character(obj@meta.data[[batch]]), row.names = colnames(obj))
    adata <- ad$AnnData(X = X, obs = obs)
    scvi$settings$seed <- as.integer(seed)
    scvi$model$SCVI$setup_anndata(adata, batch_key = if (is.null(batch)) NULL else "b")
    m <- scvi$model$SCVI(adata, n_latent = as.integer(n_latent)); m$train(max_epochs = as.integer(max_epochs))
    lat <- m$get_latent_representation(); rownames(lat) <- colnames(obj); lat
  })

# sccb_embed(obj, embedding, seed, batch, ...): returns the cells x dims matrix from SCCB_EMBEDDINGS[[embedding]]
# with cell names as rows, or NULL (with a message) if it fails, e.g. when its optional backend is missing
sccb_embed <- function(obj, embedding, seed, batch = NULL, ...) {
  tryCatch({ L <- SCCB_EMBEDDINGS[[embedding]](obj, seed, batch = batch, ...); rownames(L) <- colnames(obj); L },
           error = function(e) { message(sprintf("  [skip] embedding '%s': %s", embedding, conditionMessage(e))); NULL })
}

# ---- Five label-free metrics ----
# balance (1 - Gini of cluster sizes), jacc_min (bootstrap stability), marker_spec (marker specificity),
# rogue (ROGUE purity) and silhouette, all computed without reference labels.
# sccb_rogue_expr(obj): filtered counts matrix for ROGUE, computed once per object (NULL if ROGUE is missing)
sccb_rogue_expr <- function(obj) {
  if (!requireNamespace("ROGUE", quietly = TRUE)) return(NULL)
  ROGUE::matr.filter(as.matrix(GetAssayData(obj, assay = "RNA", layer = "counts")), min.cells = 10, min.genes = 10)
}
# .gini(x): Gini coefficient of cluster sizes (0 = all equal)
.gini <- function(x) { x <- sort(x); n <- length(x); if (n <= 1 || sum(x) == 0) return(0); sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x)) }
# .marker_spec(obj, labels): mean over clusters of the mean (pct.1 - pct.2) * |log2FC| of significant
# positive markers (FindAllMarkers, adjusted p < 0.05); 0 if there are none
.marker_spec <- function(obj, labels) {
  Idents(obj) <- factor(labels)
  m <- tryCatch(FindAllMarkers(obj, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE),
                error = function(e) NULL)
  if (is.null(m) || !nrow(m)) return(0)
  m <- m[m$p_val_adj < 0.05, , drop = FALSE]; if (!nrow(m)) return(0)
  m$score <- (m$pct.1 - m$pct.2) * abs(m$avg_log2FC)
  mean(aggregate(score ~ cluster, data = m, FUN = mean)$score)
}
# .jacc_min(method, latent, ref, k, n_boot, seed, frac): re-clusters n_boot random subsamples of frac of the
# cells with `method`; for each reference cluster takes the best-match Jaccard in every subsample, averages
# over subsamples, and returns the minimum over reference clusters (Hennig 2007 stability)
.jacc_min <- function(method, latent, ref, k, n_boot, seed, frac = 0.8) {
  ref <- as.character(ref)
  # Draw all subsamples up front: methods call set.seed() internally, which would otherwise
  # reset the RNG and make every iteration reuse the same subsample
  set.seed(seed); subs <- replicate(n_boot, sample(nrow(latent), floor(frac * nrow(latent))), simplify = FALSE)
  per <- .par(subs, function(idx) {                        # one bootstrap: best-match Jaccard of each reference cluster
    nl <- as.character(method(latent[idx, , drop = FALSE], k, seed)); rs <- ref[idx]
    vapply(unique(rs), function(cl) { a <- which(rs == cl)
      max(vapply(unique(nl), function(nc) { b <- which(nl == nc); length(intersect(a, b)) / length(union(a, b)) }, numeric(1))) },
      numeric(1)) })
  cls <- unique(ref); min(vapply(cls, function(cl) mean(unlist(lapply(per, function(v) v[names(v) == cl]))), numeric(1)))
}
# .rogue(rogue_expr, labels): mean ROGUE over clusters for the cells in rogue_expr; NA if unavailable/error
.rogue <- function(rogue_expr, labels) if (is.null(rogue_expr)) NA_real_ else tryCatch({
  l <- as.character(labels[colnames(rogue_expr)])
  mean(unlist(ROGUE::rogue(rogue_expr, labels = l, samples = rep(1, length(l)), platform = "UMI")), na.rm = TRUE)
}, error = function(e) NA_real_)
# sccb_metrics(obj, latent, labels, rogue_expr, method, k, n_boot, seed): one-row data.frame of the five
# metrics for one partition. jacc_min is NA when no method is given or n_boot < 1.
sccb_metrics <- function(obj, latent, labels, rogue_expr = NULL, method = NULL, k = length(unique(labels)),
                         n_boot = 20, seed = SEED) {
  names(labels) <- colnames(obj)
  data.frame(balance = 1 - .gini(as.numeric(table(labels))),
             jacc_min = if (is.null(method) || n_boot < 1) NA_real_ else
               tryCatch(.jacc_min(method, latent, labels, k, n_boot, seed), error = function(e) NA_real_),
             marker_spec = .marker_spec(obj, labels), rogue = .rogue(rogue_expr, labels),
             silhouette = compute_silhouette(latent, labels, max_cells = MAX_SIL_CELLS, seed = seed))
}

# ---- Benchmark: embeddings x methods x k ----
# sccb_benchmark(obj, embeddings, methods, k_range, batch, seed, n_boot, rogue, truth, latents, ...):
# returns one row per (embedding_method, target k) with actual_k, the five metrics and cluster sizes
# (plus ARI/NMI against `truth` if given). attr(, "labels") holds every partition, named
# "<embedding>_<method>_k<k>". `latents` (named list) supplies precomputed embeddings; `rogue` may be
# TRUE, FALSE or a precomputed sccb_rogue_expr() matrix; `...` goes to the embedding functions.
sccb_benchmark <- function(obj, embeddings = c("pca"), methods = c("kmeans", "gmm", "louvain", "leiden"),
                           k_range = 2:10, batch = NULL, seed = SEED, n_boot = 20, rogue = TRUE,
                           truth = NULL, latents = NULL, ...) {
  re <- if (is.matrix(rogue)) rogue else if (isTRUE(rogue)) sccb_rogue_expr(obj)   # matrix = precomputed sccb_rogue_expr()
  if (isTRUE(rogue) && is.null(re)) message("  ROGUE not installed: rogue = NA")
  jobs <- list()                                         # 1) all partitions (cheap; Louvain paths shared across k)
  for (e in embeddings) {
    L <- if (!is.null(latents[[e]])) latents[[e]] else sccb_embed(obj, e, seed, batch, ...)
    if (is.null(L)) next
    for (m in methods) {
      pth <- if (!is.null(SCCB_METHOD_PATHS[[m]])) SCCB_METHOD_PATHS[[m]](L, seed) else NULL
      for (k in k_range) {
        lab <- tryCatch(if (!is.null(pth)) pth(k) else SCCB_METHODS[[m]](L, k, seed), error = function(err) {
          message(sprintf("  [skip] %s_%s k=%d: %s", e, m, k, conditionMessage(err))); NULL })
        if (!is.null(lab)) jobs[[length(jobs) + 1]] <- list(e = e, m = m, k = k, L = L, lab = lab)
      }
    }
  }
  rows <- .par(jobs, function(j) {                      # 2) metrics per partition, in parallel (same values as serial)
    r <- data.frame(method = paste(j$e, j$m, sep = "_"), target_k = j$k, actual_k = length(unique(j$lab)),
                    sccb_metrics(obj, j$L, j$lab, re, SCCB_METHODS[[j$m]], j$k, n_boot, seed),
                    sizes_str = paste(table(j$lab), collapse = "/"))
    if (!is.null(truth)) { r$ari <- .ari(j$lab, truth); r$nmi <- .nmi(j$lab, truth) }
    message(sprintf("  %s k=%d (got %d) sil=%.3f jacc=%.3f", r$method, j$k, r$actual_k, r$silhouette, r$jacc_min))
    r }, cores = min(.ncores(), as.integer(Sys.getenv("SCCB_BENCH_WORKERS", "6"))))
  # Partitions are kept so sccb_recommend() can return the labels of the chosen configuration
  labs <- setNames(lapply(jobs, function(j) setNames(j$lab, colnames(obj))), vapply(jobs, function(j) paste0(j$e, "_", j$m, "_k", j$k), ""))
  structure(do.call(rbind, rows), labels = labs)
}

# ---- Across-embedding reproducibility curve R(k) (as in step 07) ----
# sccb_reproducibility(obj, embedding, seeds, k_range, method, batch, latents, rogue, truth,
#                      cluster_seed, ...):
# builds one embedding per seed (or takes the named list `latents`, e.g. cached scDHA runs or a
# hyperparameter family), clusters each with `method` at every k, and returns one row per k with
# mean and min pairwise ARI (R(k)), how many embeddings reached k, the smallest cluster in any
# embedding, and silhouette/ROGUE/balance and sizes of the reference (first) embedding's partition.
# Returns NULL if any embedding cannot be built.
sccb_reproducibility <- function(obj, embedding = "scdha", seeds = c(1, 2, 3, 7, 42), k_range = 2:10,
                                 method = "louvain", batch = NULL, latents = NULL, rogue = FALSE,
                                 truth = NULL, cluster_seed = SEED, ...) {
  if (is.null(latents)) {
    latents <- list()
    for (s in seeds) { L <- sccb_embed(obj, embedding, s, batch, ...); if (is.null(L)) return(NULL)
                       latents[[as.character(s)]] <- L }
  }
  ids <- Reduce(intersect, lapply(latents, rownames))     # compare partitions on the cells common to all embeddings
  if (!is.null(truth)) truth <- as.character(truth)[match(ids, colnames(obj))]
  obj <- obj[, ids]
  re <- if (rogue) sccb_rogue_expr(obj) else NULL
  lab_all <- .par(latents, function(L) {                   # per embedding: labels for every k
    Li <- L[ids, , drop = FALSE]
    f <- if (!is.null(SCCB_METHOD_PATHS[[method]])) SCCB_METHOD_PATHS[[method]](Li, cluster_seed) else function(k) SCCB_METHODS[[method]](Li, k, cluster_seed)
    lapply(k_range, function(k) setNames(f(k), ids)) })
  rows <- lapply(seq_along(k_range), function(ki) { k <- k_range[ki]
    labs <- lapply(lab_all, `[[`, ki)
    pr <- combn(length(labs), 2)
    aris <- apply(pr, 2, function(p) .ari(labs[[p[1]]], labs[[p[2]]]))
    l1 <- labs[[1]]; nk <- vapply(labs, function(l) length(unique(l)), 1L)
    d <- data.frame(k = k, actual_k = nk[1], mean_cross_seed_ari = round(mean(aris), 4),
                    min_cross_seed_ari = round(min(aris), 4),
                    silhouette = round(compute_silhouette(latents[[1]][ids, ], l1, max_cells = MAX_SIL_CELLS, seed = SEED), 4),
                    rogue = round(.rogue(re, l1), 4), balance = round(1 - .gini(as.numeric(table(l1))), 4),
                    n_embeddings_at_k = sum(nk == k),
                    min_cluster_size = min(vapply(labs, function(l) min(table(l)), 1L)),   # flags degenerate splits
                    # reference partition sizes, read by the Stage B size gate in sccb_verdict()
                    seed1_sizes = paste(sort(table(l1), decreasing = TRUE), collapse = "/"))
    if (!is.null(truth)) d$ari_vs_truth <- round(.ari(l1, truth), 4)
    message(sprintf("  R(%d)=%.3f (k got %s)", k, d$mean_cross_seed_ari, paste(nk, collapse = ",")))
    d
  })
  do.call(rbind, rows)
}

# ---- Final answer: which (configuration, k) to report ----
# sccb_recommend(scores, verdict, jacc_min, min_cells): returns list(summary, labels).
# ACCEPT: Stage B fixes the granularity k*; Stage A then ranks the admissible configurations at k*
#         (same filters, mean metric percentile) and the top one is the reported partition.
# REJECT: no finer partition is reported (labels = NULL). `verdict` row 1 (the first repro embedding)
#         decides; if several repro embeddings disagree, that is flagged in `note`.
sccb_recommend <- function(scores, verdict, jacc_min = 0.6, min_cells = 20) {
  v <- verdict[1, ]
  note <- if (length(unique(verdict$verdict)) > 1 || length(unique(na.omit(verdict$k_star))) > 1)
    "repro embeddings disagree: see res$verdict" else ""
  out <- function(method = NA_character_, score = NA_real_, msg = "")
    data.frame(verdict = v$verdict, k_star = v$k_star, method = method, score = score,
               note = trimws(paste(msg, note)))
  if (v$verdict != "ACCEPT") return(list(summary = out(msg = "no reproducible sub-structure; report as one group"), labels = NULL))
  sel <- sccb_select(scores[scores$target_k == v$k_star, ], jacc_min, min_cells)
  if (!any(sel$admissible)) return(list(summary = out(msg = "no admissible partition at k*; see res$selection"), labels = NULL))
  list(summary = out(sel$method[1], sel$score[1]),
       labels = attr(scores, "labels")[[paste0(sel$method[1], "_k", v$k_star)]])
}

# ---- One-call convenience wrapper ----
# sccb_run(obj, ...): sccb_prep, then the benchmark, Stage A, Stage B for each of repro_embeddings, and
# sccb_recommend. Returns list(scores, selection, curve, verdict, recommended, partition); with out_dir,
# each element is also written to <out_dir>/sccb_<name>.csv. `truth` may be a vector of labels or a
# metadata column name (used only for reporting ARI/NMI, never for the decision).
sccb_run <- function(obj, batch = NULL, k_range = 2:10, embeddings = c("pca"),
                     methods = c("kmeans", "gmm", "louvain", "leiden"),
                     repro_embeddings = if (requireNamespace("scDHA", quietly = TRUE)) "scdha" else "pca_perturb",
                     repro_method = "louvain",   # clustering used for R(k); any SCCB_METHODS name
                     seeds = c(1, 2, 3, 7, 42), norm = "LogNormalize", nfeatures = 3000,
                     tau = 0.8, jacc_min = 0.6, min_cells = 20, n_boot = 20, rogue = TRUE,
                     truth = NULL, out_dir = NULL, ...) {
  obj <- sccb_prep(obj, norm, nfeatures)
  if (!is.null(truth) && length(truth) == 1 && truth %in% colnames(obj@meta.data)) truth <- as.character(obj@meta.data[[truth]])
  message("== Stage 0: benchmark"); scores <- sccb_benchmark(obj, embeddings, methods, k_range, batch,
                                                             n_boot = n_boot, rogue = rogue, truth = truth, ...)
  message("== Stage A: select"); sel <- sccb_select(scores, jacc_min = jacc_min, min_cells = min_cells)
  curves <- list(); verdicts <- list()
  for (e in repro_embeddings) {
    message("== Stage B: reproducibility (", e, ")")
    cv <- sccb_reproducibility(obj, e, seeds, k_range, method = repro_method, batch = batch, truth = truth, ...)
    if (is.null(cv)) next
    curves[[e]] <- cbind(embedding = e, cv); verdicts[[e]] <- cbind(embedding = e, sccb_verdict(cv, tau))
  }
  res <- list(scores = scores, selection = sel, curve = do.call(rbind, curves), verdict = do.call(rbind, verdicts))
  if (!is.null(res$verdict)) { rec <- sccb_recommend(scores, res$verdict, jacc_min, min_cells)
    res$recommended <- rec$summary
    if (!is.null(rec$labels)) res$partition <- data.frame(cell = names(rec$labels), cluster = rec$labels) }
  if (!is.null(out_dir)) { dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    for (n in names(res)) if (!is.null(res[[n]])) write.csv(res[[n]], file.path(out_dir, paste0("sccb_", n, ".csv")), row.names = FALSE) }
  print(res$recommended)
  res
}
