# ============================================================
# pipeline/04_phase2_rod_benchmark.R
#
# Phase 2, Stage A screen on the rods. Removes cone-contaminated cells from the rod subset, recomputes a
# rod-level PCA and scDHA embedding, corrects both for sample with Harmony, and clusters each embedding
# with every method in METHODS at every k in K_RANGE. Each candidate partition is scored on silhouette,
# ROGUE, cluster-size balance, marker specificity and minimum bootstrap Jaccard. The rod-level embeddings
# are independent of the whole-retina Harmony embedding from step 01.
#
# Paper:    Phase 2 rod Stage A screen (rod panels of the Stage A figures drawn by
#           pipeline/06_publication_figures.R)
# Inputs:   output/rds/02_rods.rds
#           output/tables/03_scdha_latent.csv + output/tables/_scdha_latent.lock (optional: if both
#           exist, the saved scDHA latent is reused instead of retraining)
# Outputs:  output/rds/03_benchmarking.rds                  (scores, assignments and embeddings bundle)
#           output/tables/03_scores_matrix.csv              (silhouette and ROGUE per partition)
#           output/tables/03_scores_matrix_4metric.csv      (balance, jacc_min, marker_spec, ROGUE, sizes)
#           output/tables/03_all_cluster_assignments.csv    (one column per method x k)
#           output/tables/03_pca_latent.csv, 03_harmony_pca_latent.csv, 03_harmony_scdha_latent.csv
#           output/tables/03_scdha_latent.csv + _scdha_latent.lock (written only when scDHA is trained)
# Run:      Rscript run_pipeline.R 04
# Options:  BENCH_JACC  bootstrap resamples for jacc_min (default BENCH_JACCARD_BOOTSTRAPS in config.R;
#                       0 skips it). Step 12c_rods recomputes jacc_min with independent resamples and
#                       rewrites 03_scores_matrix_4metric.csv, so that column is the one used downstream.
# ============================================================

source("config.R")
source("R/utils.R"); source("R/scClustBench.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(Matrix)
  library(cluster)
  library(mclust)
  library(igraph)
  library(uwot)
  library(dplyr)
  library(tidyr)
  library(ROGUE)
  library(aricode)
  library(RhpcBLASctl)
})

set.seed(SEED)
# Half the cores for BLAS/OpenMP, leaving room for the forked workers used in the clustering sweep
blas_set_num_threads(max(1L, floor(N_CORES / 2)))
omp_set_num_threads(max(1L, floor(N_CORES / 2)))

T0 <- Sys.time()
# Takes a message string; prints it prefixed with the elapsed time since T0
log_step <- function(msg) {
  message(sprintf("[%s] %s",
                  format(Sys.time() - T0, digits = 3), msg))
}

# ============================================================
# Part 1: load rods and remove cone-contaminated cells
# A cell is removed if its cone module score is in the top CONE_QUANTILE_CUT tail and exceeds its rod
# score, or if cone minus rod score exceeds CONE_DELTA_CUT (same filter as load_population("Rods")).
# ============================================================
log_step("=== Part 1: Loading rods + cone filter ===")

rods <- readRDS(file.path(RDS_DIR, "02_rods.rds"))
DefaultAssay(rods) <- "RNA"
rods[["RNA"]] <- JoinLayers(rods[["RNA"]])
rods <- NormalizeData(rods, verbose = FALSE)
log_step(sprintf("  Loaded %d rod cells", ncol(rods)))

cone_g <- intersect(CONE_PANEL, rownames(rods))
rod_g  <- intersect(ROD_PANEL,  rownames(rods))
rods <- AddModuleScore(rods, features = list(cone_g), name = "cone_score",
                        seed = SEED)
rods <- AddModuleScore(rods, features = list(rod_g),  name = "rod_score",
                        seed = SEED)
rods$cone_score <- rods$cone_score1
rods$rod_score  <- rods$rod_score1

cone_q <- quantile(rods$cone_score, CONE_QUANTILE_CUT, na.rm = TRUE)
contam <- (rods$cone_score >= cone_q & rods$cone_score > rods$rod_score) |
          ((rods$cone_score - rods$rod_score) > CONE_DELTA_CUT)
n_contam <- sum(contam)
clean <- subset(rods, cells = colnames(rods)[!contam])
log_step(sprintf("  Removed %d cone-contaminated cells (%.2f%%); clean: %d",
                 n_contam, 100 * n_contam / ncol(rods), ncol(clean)))
rm(rods); gc()

# ============================================================
# Part 2: rod-level PCA and scDHA embeddings
# ============================================================
log_step("=== Part 2: PCA + scDHA on clean rods ===")

clean <- NormalizeData(clean, verbose = FALSE) |>
         FindVariableFeatures(nfeatures = 3000, verbose = FALSE) |>
         ScaleData(verbose = FALSE) |>
         RunPCA(npcs = N_PCS_CLUSTERING, verbose = FALSE)
pca_clean <- Embeddings(clean, "pca")
log_step(sprintf("  PCA: %d cells x %d PCs",
                 nrow(pca_clean), ncol(pca_clean)))

# Uncorrected rod PCA, saved for reproducibility
write.csv(data.frame(cell_id = rownames(pca_clean), pca_clean),
          file.path(TABLE_DIR, "03_pca_latent.csv"), row.names = FALSE)

# ---- scDHA latent: reuse the saved latent when both it and its lock file exist, otherwise train ----
# The lock file records the seed and dimensions of the saved latent; delete both files to retrain.
scdha_path <- file.path(TABLE_DIR, "03_scdha_latent.csv")
scdha_lock <- file.path(TABLE_DIR, "_scdha_latent.lock")
if (file.exists(scdha_path) && file.exists(scdha_lock)) {
  log_step("  Loading frozen scDHA latent from sidecar ...")
  df <- read.csv(scdha_path, check.names = FALSE,
                  stringsAsFactors = FALSE)
  ids <- as.character(df[[1]])
  scdha_clean <- as.matrix(df[, -1, drop = FALSE])
  rownames(scdha_clean) <- ids
  stopifnot(setequal(rownames(scdha_clean), colnames(clean)))
  scdha_clean <- scdha_clean[colnames(clean), , drop = FALSE]
  log_step(paste0("  Frozen lock: ",
                   paste(readLines(scdha_lock), collapse = " | ")))
} else {
  log_step("  Training scDHA from scratch ...")
  data_scdha <- t(as.matrix(GetAssayData(clean, assay = "RNA",
                                           layer = "counts")))
  data_scdha <- log2(data_scdha + 1)   # scDHA takes a cells x genes log2(count + 1) matrix

  suppressPackageStartupMessages({
    library(torch); library(scDHA)
  })
  if (!torch_is_installed()) {
    log_step("  Installing torch backend ...")
    install_torch()
  }
  set.seed(SEED)
  sc_res <- scDHA(data    = data_scdha,
                   k       = NULL,
                   method  = "scDHA",
                   sparse  = FALSE,
                   n       = SCDHA_N_GENES,
                   ncores  = N_CORES,
                   gen_fil = TRUE,
                   do.clus = TRUE,
                   seed    = SEED)
  scdha_clean <- as.matrix(sc_res$latent)
  rownames(scdha_clean) <- rownames(data_scdha)
  rm(data_scdha); gc()

  write.csv(data.frame(cell_id = rownames(scdha_clean), scdha_clean),
            scdha_path, row.names = FALSE)
  writeLines(c(paste0("trained_at=", Sys.time()),
                paste0("seed=", SEED),
                paste0("n_cells=", nrow(scdha_clean)),
                paste0("dims=", ncol(scdha_clean))),
              scdha_lock)
  log_step(sprintf("  scDHA latent: %d cells x %d dims (saved + locked)",
                    nrow(scdha_clean), ncol(scdha_clean)))
}

# ============================================================
# Part 3: Harmony correction (by sample) of the rod-level PCA and scDHA embeddings
# ============================================================
log_step("=== Part 3: Harmony on rod-level PCA + scDHA ===")

har_meta <- data.frame(sample = clean$sample, row.names = colnames(clean))
set.seed(SEED)
har_pca   <- harmony::RunHarmony(pca_clean,   har_meta, vars_use = "sample",
                                  theta = HARMONY_THETA,
                                  lambda = HARMONY_LAMBDA,
                                  max_iter = HARMONY_MAX_ITER,
                                  verbose = FALSE)
set.seed(SEED)
har_scdha <- harmony::RunHarmony(scdha_clean, har_meta, vars_use = "sample",
                                  theta = HARMONY_THETA,
                                  lambda = HARMONY_LAMBDA,
                                  max_iter = HARMONY_MAX_ITER,
                                  verbose = FALSE)
log_step("  Harmony complete on both rod-level latents")

write.csv(data.frame(cell_id = rownames(har_pca),   har_pca),
          file.path(TABLE_DIR, "03_harmony_pca_latent.csv"),
          row.names = FALSE)
write.csv(data.frame(cell_id = rownames(har_scdha), har_scdha),
          file.path(TABLE_DIR, "03_harmony_scdha_latent.csv"),
          row.names = FALSE)

cell_ids   <- colnames(clean)
counts_raw <- GetAssayData(clean, assay = "RNA", layer = "counts")
sample_vec <- as.character(clean$sample)

# ============================================================
# Part 4: ROGUE gene filter (the filtered matrix is shared by every partition)
# ============================================================
log_step("=== Part 4: ROGUE prefilter ===")

counts_dense <- as.matrix(counts_raw)
rogue_expr   <- ROGUE::matr.filter(counts_dense, min.cells = 10,
                                    min.genes = 10)
ROGUE::SE_fun(rogue_expr)
log_step(sprintf("  ROGUE-filtered: %d genes x %d cells",
                 nrow(rogue_expr), ncol(rogue_expr)))
rm(counts_dense); gc()

# ROGUE of the unsplit rod population (single label), a reference for the per-partition values
rogue_baseline <- tryCatch({
  one_lab <- rep("1", ncol(rogue_expr))
  r <- ROGUE::rogue(rogue_expr, labels = one_lab,
                    samples = rep(1, length(one_lab)),
                    platform = "UMI")
  mean(unlist(r), na.rm = TRUE)
}, error = function(e) NA_real_)
log_step(sprintf("  ROGUE baseline (unsplit rods): %.4f", rogue_baseline))

# ============================================================
# Part 5: SNN graphs on the Harmony-corrected embeddings (used by the Leiden methods)
# ============================================================
log_step("=== Part 5: Building SNN graphs (for Leiden) ===")

g_har_pca   <- build_seurat_snn_igraph(har_pca,   cell_ids,
                                        counts_for_dummy = counts_raw,
                                        reduction_name = "harpca",
                                        key_prefix = "HARPC_",
                                        dims_use = 1:N_DIMS_NEIGHBORS)
log_step(sprintf("  PCA-Harmony SNN: %d edges", igraph::ecount(g_har_pca)))

g_har_scdha <- build_seurat_snn_igraph(har_scdha, cell_ids,
                                        counts_for_dummy = counts_raw,
                                        reduction_name = "harscdha",
                                        key_prefix = "HARSC_",
                                        dims_use = 1:ncol(har_scdha))
log_step(sprintf("  scDHA-Harmony SNN: %d edges",
                 igraph::ecount(g_har_scdha)))

# ============================================================
# Helper functions for the sweep
# ============================================================
# Takes a numeric vector; returns its Gini coefficient (0 = all equal)
gini <- function(x) {
  x <- sort(x); n <- length(x)
  if (n <= 1L || sum(x) == 0) return(0)
  sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
}
# Takes cluster labels; returns 1 - Gini of the cluster sizes (1 = perfectly balanced)
balance_of <- function(l) 1 - gini(as.numeric(table(l)))

# Takes a Seurat object and labels; returns the mean over clusters of the mean (pct.1 - pct.2) * |log2FC|
# of significant positive markers (0 if none)
marker_specificity <- function(seurat_obj, labels) {
  Idents(seurat_obj) <- factor(labels)
  m <- tryCatch(
    FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25,
                   logfc.threshold = 0.5, verbose = FALSE),
    error = function(e) NULL)
  if (is.null(m) || nrow(m) == 0) return(0)
  m <- m[m$p_val_adj < 0.05, , drop = FALSE]
  if (nrow(m) == 0) return(0)
  m$score <- (m$pct.1 - m$pct.2) * abs(m$avg_log2FC)
  per_cl <- aggregate(score ~ cluster, data = m, FUN = mean)
  mean(per_cl$score)
}

# Takes a method name, a subsampled latent matrix and a target k; returns integer labels from re-running
# that method on the subsample (used by the bootstrap Jaccard)
recluster_subsample <- function(method, latent_sub, target_k) {
  if (target_k == 1) return(rep(1L, nrow(latent_sub)))
  if (method %in% c("PCA_KMeans", "scDHA_KMeans")) {
    set.seed(SEED)
    return(kmeans(latent_sub, centers = target_k,
                  nstart = 10, iter.max = 100)$cluster)
  }
  if (method == "PCA_GMM") {
    g <- tryCatch(
      Mclust(latent_sub[, 1:min(20, ncol(latent_sub))],
              G = target_k, verbose = FALSE),
      error = function(e) NULL)
    if (is.null(g)) return(rep(1L, nrow(latent_sub)))
    return(as.integer(g$classification))
  }
  if (method %in% c("PCA_LL", "scDHA_LL")) {
    obj <- CreateSeuratObject(matrix(0, 2, nrow(latent_sub),
                              dimnames = list(c("d1", "d2"), rownames(latent_sub))))
    obj[["lat"]] <- CreateDimReducObject(embeddings = latent_sub,
                                          key = "L_", assay = "RNA")
    obj <- FindNeighbors(obj, reduction = "lat",
                          dims = 1:ncol(latent_sub), verbose = FALSE)
    best_lab <- NULL; best_diff <- Inf
    # Log-spaced resolution grid; keep the labels whose cluster count is closest to target_k
    for (r in exp(seq(log(1e-3), log(2), length.out = 30))) {
      set.seed(SEED)
      obj <- FindClusters(obj, resolution = r, algorithm = 1,
                           verbose = FALSE)
      lab <- as.integer(Idents(obj)); k_now <- length(unique(lab))
      if (abs(k_now - target_k) < best_diff) {
        best_diff <- abs(k_now - target_k); best_lab <- lab
      }
      if (best_diff == 0) break
    }
    return(best_lab)
  }
  if (method %in% c("PCA_Leiden", "scDHA_Leiden")) {
    obj <- CreateSeuratObject(matrix(0, 2, nrow(latent_sub),
                              dimnames = list(c("d1", "d2"), rownames(latent_sub))))
    obj[["lat"]] <- CreateDimReducObject(embeddings = latent_sub,
                                          key = "L_", assay = "RNA")
    obj <- FindNeighbors(obj, reduction = "lat",
                          dims = 1:ncol(latent_sub), verbose = FALSE)
    # Convert the Seurat SNN matrix to a weighted igraph (upper triangle only, one edge per pair)
    snn_mat <- methods::as(obj@graphs[["RNA_snn"]], "TsparseMatrix")
    keep    <- snn_mat@i < snn_mat@j
    edges   <- cbind(snn_mat@i[keep] + 1L, snn_mat@j[keep] + 1L)
    g_sub   <- igraph::make_graph(t(edges), n = ncol(snn_mat),
                                    directed = FALSE)
    igraph::E(g_sub)$weight <- snn_mat@x[keep]
    igraph::V(g_sub)$name   <- rownames(snn_mat)
    res <- find_resolution_for_k(g_sub, target_k = target_k,
                                  algorithm = "leiden",
                                  res_min = 1e-6, res_max = 5,
                                  max_iter = 20, seed = SEED)
    return(as.integer(res$labels))
  }
  rep(1L, nrow(latent_sub))
}

# Takes a method, latent, reference labels and target k; returns the minimum over reference clusters
# of the mean best-match Jaccard across n_iter subsamples of frac of the cells. All subsamples are
# drawn before any reclustering, because recluster_subsample() reseeds the RNG. Step 12c computes the
# same quantity with .jacc_min() in R/scClustBench.R.
per_cluster_jaccard <- function(method, latent, ref_labels, target_k,
                                 n_iter = JACCARD_BOOTSTRAPS, frac = 0.8) {
  n <- nrow(latent)
  ref <- as.character(ref_labels)
  per_iter_max <- vector("list", length(unique(ref)))
  names(per_iter_max) <- unique(ref)
  for (cl in names(per_iter_max)) per_iter_max[[cl]] <- numeric(0)
  set.seed(SEED)
  subsamples <- lapply(seq_len(n_iter), function(i) sample(seq_len(n), floor(frac * n)))
  for (it in seq_len(n_iter)) {
    samp_idx <- subsamples[[it]]
    sub_lat  <- latent[samp_idx, , drop = FALSE]
    new_lab  <- recluster_subsample(method, sub_lat, target_k)
    new_lab  <- as.character(new_lab)
    ref_sub  <- ref[samp_idx]
    for (cl in unique(ref_sub)) {
      ref_set <- which(ref_sub == cl)
      jac <- vapply(unique(new_lab), function(nc) {
        new_set <- which(new_lab == nc)
        length(intersect(ref_set, new_set)) /
          length(union(ref_set, new_set))
      }, numeric(1))
      per_iter_max[[cl]] <- c(per_iter_max[[cl]], max(jac))
    }
  }
  per_cluster_mean <- sapply(per_iter_max, mean)
  if (length(per_cluster_mean) == 0) return(NA_real_)
  min(per_cluster_mean)
}

# ============================================================
# Part 6: cluster every (method, k) and score each partition
# ============================================================
log_step("=== Part 6: Clustering sweep + metrics ===")

sil_rogue_rows  <- list()
fourmet_rows    <- list()
all_assignments <- data.frame(cell_id = cell_ids)

# Takes a method name; returns the Harmony-corrected latent it clusters (scDHA or PCA)
method_to_latent <- function(m) {
  if (grepl("scDHA", m)) har_scdha else har_pca
}
# Takes a method name; returns the matching SNN graph (Leiden methods)
method_to_graph <- function(m) {
  if (grepl("scDHA", m)) g_har_scdha else g_har_pca
}

# Writes the score tables and assignments accumulated so far; returns nothing
checkpoint <- function() {
  if (length(sil_rogue_rows) > 0) {
    write.csv(do.call(rbind, sil_rogue_rows),
              file.path(TABLE_DIR, "03_scores_matrix.csv"),
              row.names = FALSE)
  }
  if (length(fourmet_rows) > 0) {
    write.csv(do.call(rbind, fourmet_rows),
              file.path(TABLE_DIR, "03_scores_matrix_4metric.csv"),
              row.names = FALSE)
  }
  write.csv(all_assignments,
            file.path(TABLE_DIR, "03_all_cluster_assignments.csv"),
            row.names = FALSE)
}

# Methods run in parallel (one worker each). Each Louvain method builds its SNN graph once and reuses the
# resolution path across k (same labels as a separate search per k). BENCH_JACC=0 skips the bootstrap
# Jaccard here; step 12c recomputes that column with independent resamples.
BENCH_JACC <- as.integer(Sys.getenv("BENCH_JACC", BENCH_JACCARD_BOOTSTRAPS))
# Takes a method name; returns list(sr = silhouette/ROGUE rows, fm = 4-metric rows, asg = labels per k)
run_method <- function(mname) {
  sil_rogue_rows <- list(); fourmet_rows <- list(); asg <- list(); lpath <- NULL
  lat <- method_to_latent(mname)
  for (k in K_RANGE) {
    log_step(sprintf("--- %s @ k=%d ---", mname, k))
    labels <- NULL

    tryCatch({
      if (mname == "PCA_KMeans") {
        set.seed(SEED)
        labels <- kmeans(lat, centers = k, nstart = 25,
                         iter.max = 100)$cluster
      } else if (mname == "PCA_GMM") {
        set.seed(SEED)
        gmm <- tryCatch(
          Mclust(lat[, 1:min(20, ncol(lat))], G = k, verbose = FALSE),
          error = function(e) NULL)
        labels <- if (is.null(gmm)) rep(1L, nrow(lat)) else
                  as.integer(gmm$classification)
      } else if (mname == "scDHA_KMeans") {
        set.seed(SEED)
        labels <- kmeans(lat, centers = k, nstart = 25,
                         iter.max = 100)$cluster
      } else if (mname %in% c("PCA_LL", "scDHA_LL")) {
        if (is.null(lpath)) lpath <- sccb_louvain_path(lat, SEED)
        labels <- lpath(k)
        log_step(sprintf("  Louvain actual_k=%d", length(unique(labels))))
      } else if (mname %in% c("PCA_Leiden", "scDHA_Leiden")) {
        g <- method_to_graph(mname)
        res <- find_resolution_for_k(g, target_k = k, algorithm = "leiden",
                                      res_min = 1e-6, res_max = 5,
                                      max_iter = 30, seed = SEED)
        labels <- as.integer(res$labels[cell_ids])
        log_step(sprintf("  Leiden res=%.4g, actual_k=%d",
                          res$resolution, res$actual_k))
      }
      names(labels) <- cell_ids
      actual_k <- length(unique(labels))
      col <- paste0(mname, "_k", k)
      asg[[col]] <- labels[cell_ids]

      # Internal metrics for this partition
      sil <- compute_silhouette(lat, labels,
                                 max_cells = MAX_SIL_CELLS, seed = SEED)
      rg <- tryCatch({
        rogue_lab <- labels[colnames(rogue_expr)]
        r <- ROGUE::rogue(rogue_expr,
                          labels   = as.character(rogue_lab),
                          samples  = rep(1, length(rogue_lab)),
                          platform = "UMI")
        mean(unlist(r), na.rm = TRUE)
      }, error = function(e) {
        log_step(sprintf("  ROGUE error: %s", e$message)); NA_real_ })
      bal <- balance_of(labels)
      ms  <- marker_specificity(clean, labels)
      jac <- tryCatch(
        if (BENCH_JACC < 1) NA_real_ else per_cluster_jaccard(mname, lat, labels, target_k = k,
                             n_iter = BENCH_JACC),
        error = function(e) {
          log_step(sprintf("  Jaccard error: %s", e$message)); NA_real_ })

      log_step(sprintf("  sil=%.4f rogue=%.4f bal=%.4f jacc=%.4f ms=%.4f",
                        sil, rg, bal, jac, ms))

      sil_rogue_rows[[length(sil_rogue_rows) + 1L]] <- data.frame(
        method = mname, target_k = k, actual_k = actual_k,
        silhouette = round(sil, 4), rogue = round(rg, 4),
        stringsAsFactors = FALSE)
      fourmet_rows[[length(fourmet_rows) + 1L]] <- data.frame(
        method = mname, target_k = k, actual_k = actual_k,
        balance = round(bal, 4), jacc_min = round(jac, 4),
        marker_spec = round(ms, 4), rogue = round(rg, 4),
        sizes_str = paste(table(labels), collapse = "/"),
        stringsAsFactors = FALSE)

    }, error = function(e) {
      # A failed partition is kept as an NA row so every (method, k) appears in the tables
      log_step(sprintf("  ERROR for %s k=%d: %s", mname, k, e$message))
      sil_rogue_rows[[length(sil_rogue_rows) + 1L]] <<- data.frame(
        method = mname, target_k = k, actual_k = NA,
        silhouette = NA, rogue = NA, stringsAsFactors = FALSE)
      fourmet_rows[[length(fourmet_rows) + 1L]] <<- data.frame(
        method = mname, target_k = k, actual_k = NA,
        balance = NA, jacc_min = NA, marker_spec = NA, rogue = NA,
        sizes_str = NA_character_, stringsAsFactors = FALSE)
    })
  }
  list(sr = sil_rogue_rows, fm = fourmet_rows, asg = asg)
}
res_m <- .par(METHODS, run_method)
# Merge the per-method results in METHODS order
for (r in res_m) { sil_rogue_rows <- c(sil_rogue_rows, r$sr); fourmet_rows <- c(fourmet_rows, r$fm)
                   for (nm in names(r$asg)) all_assignments[[nm]] <- r$asg[[nm]] }
checkpoint()

# ============================================================
# Part 7: write final tables and the RDS bundle used by step 05
# ============================================================
log_step("=== Part 7: Saving final outputs ===")

scores_silrogue <- do.call(rbind, sil_rogue_rows)
scores_4metric  <- do.call(rbind, fourmet_rows)

write.csv(scores_silrogue,
          file.path(TABLE_DIR, "03_scores_matrix.csv"), row.names = FALSE)
write.csv(scores_4metric,
          file.path(TABLE_DIR, "03_scores_matrix_4metric.csv"),
          row.names = FALSE)
write.csv(all_assignments,
          file.path(TABLE_DIR, "03_all_cluster_assignments.csv"),
          row.names = FALSE)

bundle <- list(
  scores_silrogue     = scores_silrogue,
  scores_4metric      = scores_4metric,
  cluster_assignments = all_assignments,
  har_pca             = har_pca,
  har_scdha           = har_scdha,
  scdha_clean         = scdha_clean,
  pca_clean           = pca_clean,
  cell_ids            = cell_ids,
  rogue_baseline      = rogue_baseline,
  n_cells             = length(cell_ids),
  n_contam_removed    = n_contam
)
saveRDS(bundle, file.path(RDS_DIR, "03_benchmarking.rds"))

log_step("=== Step 04 done ===")
log_step(sprintf("Total wall-time: %s",
                  format(Sys.time() - T0, digits = 3)))
