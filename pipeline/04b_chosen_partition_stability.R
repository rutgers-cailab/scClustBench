# ============================================================
# pipeline/04b_chosen_partition_stability.R
#
# Bootstrap stability of scDHA + Louvain partitions within one fixed embedding. For each k, the full
# embedding is clustered to a reference partition; then JACCARD_BOOTSTRAPS subsamples of 80% of the cells
# (without replacement) are re-clustered to the same k. Reports the per-cluster mean best-match Jaccard
# (ordered by cluster size), its minimum, and the mean and SD of the subsample-vs-reference ARI. This is
# the within-embedding test; step 04c tests reproducibility across scDHA training seeds.
#
# Paper:    Supplementary Figure S6 (rod table, drawn by pipeline/04d_stability_figures.R); bipolar
#           bootstrap stability (POP=BC)
# Inputs:   POP=Rods: output/tables/03_harmony_scdha_latent.csv (step 04)
#           POP=BC:   output/rds/07_bc_positive_control.rds (step 07_bc; seed-1 scDHA + Harmony embedding)
# Outputs:  output/tables/04b_chosen_partition_stability.csv   (POP=Rods)
#           output/tables/04b_bc_partition_stability.csv       (POP=BC)
# Run:      POP=Rods Rscript pipeline/04b_chosen_partition_stability.R   (Rscript run_pipeline.R 04b)
#           POP=BC Rscript pipeline/04b_chosen_partition_stability.R     (Rscript run_pipeline.R 04b_bc)
# Options:  POP     Rods or BC (default Rods)
#           K_LIST  comma-separated k values for POP=BC (default 9,16); rods use STABILITY_K_RANGE
# ============================================================

source("config.R")
source("R/utils.R"); source("R/scClustBench.R")

suppressPackageStartupMessages({
  library(Seurat); library(igraph); library(aricode)
})

T0 <- Sys.time()
# Takes a message; prints it with the elapsed time
log_step <- function(m) message(sprintf("[%s] %s", format(Sys.time() - T0, digits = 3), m))

# ---- Load the embedding ----
# POP=BC: the bipolar seed-1 scDHA + Harmony embedding from step 07, at the k values in K_LIST.
# Otherwise: the rod Harmony-corrected scDHA latent from step 04, at STABILITY_K_RANGE.
POP <- Sys.getenv("POP", "Rods")
if (POP == "BC") {
  lat <- readRDS(file.path(RDS_DIR, "07_bc_positive_control.rds"))$seed_lat[["1"]]
  STABILITY_K_RANGE <- as.integer(strsplit(Sys.getenv("K_LIST", "9,16"), ",")[[1]])
} else {
  latent_path <- file.path(TABLE_DIR, "03_harmony_scdha_latent.csv")
  stopifnot(file.exists(latent_path))
  df  <- read.csv(latent_path, check.names = FALSE)
  lat <- as.matrix(df[, -1]); rownames(lat) <- as.character(df[[1]])
}
log_step(sprintf("Loaded harmony-scDHA latent: %d x %d", nrow(lat), ncol(lat)))

# Takes a latent matrix and target k; returns named integer Louvain labels from the resolution (searched on
# a log grid of grid_n points in [lo, hi]) whose cluster count is closest to target_k, as in step 04
louvain_to_k <- function(latent, target_k, grid_n = 60, lo = 1e-4, hi = 10) {
  obj <- CreateSeuratObject(matrix(0, 2, nrow(latent), dimnames = list(c("d1", "d2"), rownames(latent))))
  obj[["lat"]] <- CreateDimReducObject(embeddings = latent, key = "L_", assay = "RNA")
  obj <- FindNeighbors(obj, reduction = "lat", dims = 1:ncol(latent), verbose = FALSE)
  bl <- NULL; bd <- Inf
  for (r in exp(seq(log(lo), log(hi), length.out = grid_n))) {
    set.seed(SEED); obj <- FindClusters(obj, resolution = r, algorithm = 1, verbose = FALSE)
    lab <- as.integer(Idents(obj)); nk <- length(unique(lab))
    if (abs(nk - target_k) < bd) { bd <- abs(nk - target_k); bl <- lab }
    if (bd == 0) break
  }
  names(bl) <- rownames(latent); bl
}

# Takes a latent and target k; returns a one-row data.frame of bootstrap stability results
stability_at_k <- function(latent, target_k, n_iter = JACCARD_BOOTSTRAPS, frac = 0.8) {
  ref  <- as.character(louvain_to_k(latent, target_k))
  names(ref) <- rownames(latent)
  refsz <- sort(table(ref), decreasing = TRUE)
  n <- nrow(latent)
  per_cl <- setNames(vector("list", length(unique(ref))), unique(ref))
  aris <- numeric(n_iter)
  # Draw all subsamples before re-clustering: louvain_to_k() calls set.seed(SEED) internally, so
  # drawing inside the loop would reuse one subsample for every iteration.
  # Subsamples are re-clustered on a coarser 30-point resolution grid to keep the bootstrap affordable.
  set.seed(SEED)
  subs <- replicate(n_iter, sample(seq_len(n), floor(frac * n)), simplify = FALSE)
  boot <- .par(subs, function(idx) {                      # bootstraps in parallel; subsamples fixed above
    nl <- as.character(louvain_to_k(latent[idx, , drop = FALSE], target_k, grid_n = 30, lo = 1e-3, hi = 2))
    rs <- ref[idx]
    # per reference cluster: best Jaccard overlap with any subsample cluster
    list(ari = ARI(rs, nl), jac = vapply(unique(rs), function(cl) { rset <- which(rs == cl)
      max(vapply(unique(nl), function(nc) { ns <- which(nl == nc); length(intersect(rset, ns)) / length(union(rset, ns)) }, numeric(1))) },
      numeric(1)))
  })
  aris <- vapply(boot, `[[`, numeric(1), "ari")
  for (b in boot) for (cl in names(b$jac)) per_cl[[cl]] <- c(per_cl[[cl]], b$jac[[cl]])
  pcm <- sapply(per_cl, mean)
  # order per-cluster jaccard by descending cluster size
  cl_by_size <- names(sort(table(ref), decreasing = TRUE))
  data.frame(
    method = CHOSEN_METHOD, target_k = target_k,
    sizes = paste(as.integer(refsz), collapse = " / "),
    jacc_per_cluster = paste(round(pcm[cl_by_size], 4), collapse = " / "),
    jacc_min = round(min(pcm), 4),
    mean_ari = round(mean(aris), 4), sd_ari = round(sd(aris), 4),
    n_iter = n_iter, stringsAsFactors = FALSE)
}

# ---- Run each k and write the table ----
rows <- list()
for (k in STABILITY_K_RANGE) {
  log_step(sprintf("=== bootstrap stability @ k=%d (n_iter=%d) ===", k, JACCARD_BOOTSTRAPS))
  rows[[length(rows) + 1]] <- stability_at_k(lat, k)
}
out <- do.call(rbind, rows)
write.csv(out, file.path(TABLE_DIR, if (POP == "BC") "04b_bc_partition_stability.csv" else "04b_chosen_partition_stability.csv"), row.names = FALSE)
log_step("=== STABILITY TABLE ==="); print(out)
log_step("DONE")
