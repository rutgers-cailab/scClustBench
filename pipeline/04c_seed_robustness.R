# ============================================================
# pipeline/04c_seed_robustness.R
#
# Tests whether the illustrative rod partition (scDHA + Louvain at CHOSEN_K) is a property of the data or
# of one scDHA training run. The cone-cleaned rods are embedded with scDHA under each seed in SEED_LIST,
# corrected with Harmony, and clustered with Louvain at CHOSEN_K; only the scDHA seed varies. Reports the
# cluster sizes per seed and the pairwise ARI and NMI between seeds. When step 07_rods has run, its five
# rod scDHA + Harmony embeddings (same cells, seeds and settings) are reused instead of retraining.
#
# Paper:    Supplementary Figure S6 (drawn by pipeline/04d_stability_figures.R)
# Inputs:   output/rds/02_rods.rds
#           output/rds/07_rods_positive_control.rds (optional; step 07_rods embeddings)
# Outputs:  output/tables/04c_seed_sizes.csv               (cluster sizes per seed at CHOSEN_K)
#           output/tables/04c_seed_ari_matrix.csv          (pairwise ARI between seeds)
#           output/tables/04c_seed_robustness_summary.csv  (mean/min/max pairwise ARI, mean NMI)
#           output/rds/04c_seed_robustness.rds
# Run:      Rscript run_pipeline.R 04c
# ============================================================

source("config.R")
source("R/utils.R")

suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(Matrix)
  library(igraph); library(aricode)
})

SEED_LIST <- c(1, 2, 3, 7, 42)   # scDHA training seeds to compare (same as SCDHA_SEED_LIST in config.R)

T0 <- Sys.time()
# Takes a message; prints it with the elapsed time
log_step <- function(m) message(sprintf("[%s] %s", format(Sys.time() - T0, digits = 3), m))

# ---- Load rods and apply the cone-contamination filter (same as step 04) ----
log_step("Loading rods + cone filter")
rods <- readRDS(file.path(RDS_DIR, "02_rods.rds"))
DefaultAssay(rods) <- "RNA"
rods[["RNA"]] <- JoinLayers(rods[["RNA"]])
rods <- NormalizeData(rods, verbose = FALSE)

cone_g <- intersect(CONE_PANEL, rownames(rods))
rod_g  <- intersect(ROD_PANEL,  rownames(rods))
rods <- AddModuleScore(rods, features = list(cone_g), name = "cone_score", seed = SEED)
rods <- AddModuleScore(rods, features = list(rod_g),  name = "rod_score",  seed = SEED)
rods$cone_score <- rods$cone_score1; rods$rod_score <- rods$rod_score1
cone_q <- quantile(rods$cone_score, CONE_QUANTILE_CUT, na.rm = TRUE)
contam <- (rods$cone_score >= cone_q & rods$cone_score > rods$rod_score) |
          ((rods$cone_score - rods$rod_score) > CONE_DELTA_CUT)
clean <- subset(rods, cells = colnames(rods)[!contam])
log_step(sprintf("Clean rod cells: %d", ncol(clean)))
rm(rods); gc(verbose = FALSE)

har_meta <- data.frame(sample = clean$sample, row.names = colnames(clean))
# scDHA input: cells x genes log2(count + 1)
counts_log2 <- log2(t(as.matrix(GetAssayData(clean, assay = "RNA", layer = "counts"))) + 1)

# Takes a latent matrix and target k; returns named integer Louvain labels whose cluster count is closest
# to target_k over a 60-point log resolution grid (same search as step 04)
louvain_to_k <- function(latent, target_k) {
  obj <- CreateSeuratObject(matrix(0, 2, nrow(latent), dimnames = list(c("d1", "d2"), rownames(latent))))
  obj[["lat"]] <- CreateDimReducObject(embeddings = latent, key = "L_", assay = "RNA")
  obj <- FindNeighbors(obj, reduction = "lat", dims = 1:ncol(latent), verbose = FALSE)
  bl <- NULL; bd <- Inf
  for (r in exp(seq(log(1e-4), log(10), length.out = 60))) {
    set.seed(SEED); obj <- FindClusters(obj, resolution = r, algorithm = 1, verbose = FALSE)
    lab <- as.integer(Idents(obj)); nk <- length(unique(lab))
    if (abs(nk - target_k) < bd) { bd <- abs(nk - target_k); bl <- lab }
    if (bd == 0) break
  }
  names(bl) <- rownames(latent); bl
}

suppressPackageStartupMessages({ library(torch); library(scDHA) })
if (!torch_is_installed()) install_torch()

# ---- One embedding per seed: reuse step 07_rods embeddings if present, else train scDHA + Harmony ----
pc_f <- file.path(RDS_DIR, "07_rods_positive_control.rds")
pre <- if (file.exists(pc_f)) readRDS(pc_f)$seed_lat else NULL
if (!is.null(pre)) log_step("Using the five stage-07 rod embeddings (no retraining)")
parts <- list(); size_rows <- list()
for (s in SEED_LIST) {
  log_step(sprintf("=== scDHA seed %d ===", s))
  if (!is.null(pre)) {        # step 07_rods embeddings: same cells, seeds and settings
    har <- pre[[as.character(s)]]; stopifnot(identical(rownames(har), rownames(counts_log2)))
  } else {
    set.seed(s)
    sc <- scDHA(data = counts_log2, k = NULL, method = "scDHA", sparse = FALSE,
                n = SCDHA_N_GENES, ncores = N_CORES, gen_fil = TRUE,
                do.clus = TRUE, seed = s)
    lat <- as.matrix(sc$latent); rownames(lat) <- rownames(counts_log2)
    set.seed(s)
    har <- harmony::RunHarmony(lat, har_meta, vars_use = "sample",
                               theta = HARMONY_THETA, lambda = HARMONY_LAMBDA,
                               max_iter = HARMONY_MAX_ITER, verbose = FALSE)
  }
  lab <- louvain_to_k(har, CHOSEN_K)
  parts[[as.character(s)]] <- lab
  sz <- sort(as.integer(table(lab)), decreasing = TRUE)
  size_rows[[length(size_rows)+1]] <- data.frame(
    seed = s, k = length(unique(lab)),
    sizes = paste(sz, collapse = " / "),
    smallest = min(sz), smallest_pct = round(100*min(sz)/length(lab), 2))
  log_step(sprintf("  seed %d k=%d sizes: %s", s, CHOSEN_K, paste(sz, collapse = " / ")))
}

size_df <- do.call(rbind, size_rows)
write.csv(size_df, file.path(TABLE_DIR, "04c_seed_sizes.csv"), row.names = FALSE)

# ---- Pairwise ARI and NMI between seeds, on the cells common to all (the same cells for every seed) ----
ids <- Reduce(intersect, lapply(parts, names))
sl <- names(parts)
ariM <- matrix(NA, length(sl), length(sl), dimnames = list(sl, sl))
nmiM <- ariM
for (i in seq_along(sl)) for (j in seq_along(sl)) {
  ariM[i,j] <- ARI(parts[[i]][ids], parts[[j]][ids])
  nmiM[i,j] <- NMI(parts[[i]][ids], parts[[j]][ids])
}
write.csv(data.frame(seed = sl, round(ariM, 4)),
          file.path(TABLE_DIR, "04c_seed_ari_matrix.csv"), row.names = FALSE)

# Summary over distinct seed pairs (upper triangle, diagonal excluded)
offdiag <- ariM[upper.tri(ariM)]
summary_df <- data.frame(
  metric = c("n_seeds","mean_pairwise_ARI","min_pairwise_ARI","max_pairwise_ARI","mean_pairwise_NMI"),
  value  = c(length(sl), round(mean(offdiag),4), round(min(offdiag),4),
             round(max(offdiag),4), round(mean(nmiM[upper.tri(nmiM)]),4)))
write.csv(summary_df, file.path(TABLE_DIR, "04c_seed_robustness_summary.csv"), row.names = FALSE)

log_step("=== SEED ROBUSTNESS SUMMARY ===")
print(size_df); cat("\nPairwise ARI:\n"); print(round(ariM,3)); cat("\n"); print(summary_df)
saveRDS(list(parts=parts, sizes=size_df, ari=ariM, nmi=nmiM, summary=summary_df),
        file.path(RDS_DIR, "04c_seed_robustness.rds"))
log_step("DONE")
