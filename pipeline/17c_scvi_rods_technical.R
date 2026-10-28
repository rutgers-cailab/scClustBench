# ============================================================
# pipeline/17c_scvi_rods_technical.R
#
# Tests whether the scVI rod k = 3 split (13s/13t; R = 0.84, but one of the five
# embeddings does not reach k = 3, so it fails the validity gates) is technical. The
# split separates a minority of rods with higher sequencing depth and
# ribosomal fraction, plus a small third cluster. The scVI R(k) curve and
# Stage B verdict are recomputed under four conditions:
#   base        as in 13s (reference)
#   thinned     counts binomially thinned to one common depth (20th percentile
#               of UMI per cell; shallower cells dropped)
#   covariates  scVI with log(UMI) and ribosomal fraction as continuous
#               covariates (sample stays the batch key)
#   no_tech     ribosomal, mitochondrial, sex-linked, Malat1 and Gm lncRNA genes
#               removed from the HVGs, and the small third cluster (13t cluster 3) removed
# If R(k) falls below tau, or the minority group is no longer recovered
# (column ari_vs_13t_axis), once depth and composition are controlled, the
# split is technical; if it persists, it is reported as an axis of
# unexplained origin.
#
# Paper:    scVI rod split, technical-covariate test
# Inputs:   output/rds/07_rods_positive_control.rds
#           output/tables/13t_scvi_rods_assignments.csv
# Outputs:  output/tables/17c_scvi_rods_technical_curves.csv
#           output/tables/17c_scvi_rods_technical_verdicts.csv
#           output/tables/17c_parts/<condition>_{curve,verdict}.csv and <condition>_seed<s>.rds
# Run:      Rscript run_pipeline.R 17c
#           or  RETICULATE_PYTHON=~/micromamba/envs/scvi/bin/python Rscript pipeline/17c_scvi_rods_technical.R
# Options:  RETICULATE_PYTHON    python with scvi-tools (required unless all seeds are cached)
#           COND                 run only this condition (base, thinned, covariates, no_tech); default: all
#           SEED_ONLY            with COND: train and save only this seed, then exit (default: unset)
#           COMBINE              merge the 17c_parts tables into the final tables, no training (default: unset)
#           SLURM_CPUS_PER_TASK  torch threads in SEED_ONLY mode (default 8)
# Runtime:  about 4 h on CPU (4 conditions x 5 scVI trainings, 200 epochs, plus the
#           Louvain search); about 1 h on a GPU. RAM about 32 GB.
# ============================================================
source("config.R"); source("R/scClustBench.R")
suppressPackageStartupMessages({ library(aricode); library(Matrix) })
KS <- 2:6

# ---- Rod data and per-cell technical covariates ----
pc  <- readRDS(file.path(RDS_DIR, "07_rods_positive_control.rds"))
raw <- pc$bc
if (inherits(raw[["RNA"]], "Assay5")) raw[["RNA"]] <- JoinLayers(raw[["RNA"]])
# Rounded because two samples were deposited as non-integer (ambient-corrected) counts
cnt <- round(GetAssayData(raw, assay = "RNA", layer = "counts"))
raw$log_umi  <- log1p(colSums(cnt))
raw$pct_ribo <- colSums(cnt[grepl("^Rp[sl]", rownames(cnt)), ]) / colSums(cnt)

# ---- Reference scVI split from 13t (seed-1 partition) ----
# Column 2 is the seed-1 partition. Louvain label numbers are arbitrary, so the clusters are
# identified by size: the second largest is the minority group, the smallest the small third group
ref <- read.csv(file.path(TABLE_DIR, "13t_scvi_rods_assignments.csv"), check.names = FALSE)
by_size <- names(sort(table(ref[[2]]), decreasing = TRUE))
ref_axis <- setNames(ifelse(ref[[2]] == by_size[2], "minor", "major"), ref$cell_id)
ref_k3   <- setNames(ref[[2]], ref$cell_id)

# thin(x, target, seed): sparse count matrix -> binomially thinned copy with about `target`
# UMIs per cell; cells below target are dropped
thin <- function(x, target, seed) {
  x <- x[, colSums(x) >= target]; set.seed(seed)
  # j: column index of every stored non-zero, so each count is thinned with its cell's keep probability
  p <- target / colSums(x); j <- rep(seq_len(ncol(x)), diff(x@p))
  x@x <- as.numeric(rbinom(length(x@x), x@x, p[j])); drop0(x)
}
# Technical genes removed from the HVGs in no_tech: ribosomal, mitochondrial, sex-linked
# (Xist, Tsix and Y genes), Malat1 and Gm-prefixed lncRNAs
TECH <- "^Rp[sl]|^mt-|^Xist$|^Tsix$|^Ddx3y$|^Eif2s3y$|^Uty$|^Kdm5d$|^Malat1$|^Gm[0-9]"

# scvi_cov(obj, seed, batch, ...): latent matrix (cells x n_latent) from scVI with continuous
# covariates; same as the registry 'scvi' embedding plus continuous_covariate_keys
SCCB_EMBEDDINGS$scvi_cov <- function(obj, seed, batch = NULL, n_latent = 10, max_epochs = 200,
                                     covs = c("log_umi", "pct_ribo"), ...) {
  scvi <- reticulate::import("scvi"); ad <- reticulate::import("anndata")
  X <- round(Matrix::t(GetAssayData(obj, assay = "RNA", layer = "counts")[VariableFeatures(obj), ]))
  obs <- data.frame(b = as.character(obj@meta.data[[batch]]), obj@meta.data[, covs], row.names = colnames(obj))
  adata <- ad$AnnData(X = X, obs = obs); scvi$settings$seed <- as.integer(seed)
  scvi$model$SCVI$setup_anndata(adata, batch_key = "b", continuous_covariate_keys = as.list(covs))
  m <- scvi$model$SCVI(adata, n_latent = as.integer(n_latent)); m$train(max_epochs = as.integer(max_epochs))
  lat <- m$get_latent_representation(); rownames(lat) <- colnames(obj); lat
}

# ---- The four conditions: a prepared Seurat object and the embedding to train on it ----
cond <- list(
  base       = list(obj = sccb_prep(raw, "LogNormalize", 3000), emb = "scvi"),
  thinned    = local({ tgt <- floor(quantile(colSums(cnt), 0.2)); message("thinning target ", tgt, " UMI")
                       x <- thin(cnt, tgt, SEED)
                       list(obj = sccb_prep(CreateSeuratObject(x, meta.data = raw@meta.data[colnames(x), ]), "LogNormalize", 3000),
                            emb = "scvi") }),
  covariates = list(obj = sccb_prep(raw, "LogNormalize", 3000), emb = "scvi_cov"),
  # Keeps only cells present in the 13t table that are not in its smallest cluster
  no_tech    = local({ o <- sccb_prep(raw[, colnames(raw) %in% names(ref_k3)[ref_k3 != by_size[3]]], "LogNormalize", 3000)
                       VariableFeatures(o) <- grep(TECH, VariableFeatures(o), value = TRUE, invert = TRUE)
                       list(obj = o, emb = "scvi") }))

# ---- COMBINE mode: merge per-condition parts and exit ----
# Parallel use: COND=<name> runs one condition and writes 17c_parts/<name>_{curve,verdict}.csv;
# COMBINE=1 merges the parts into the final tables (no training).
PARTS <- file.path(TABLE_DIR, "17c_parts"); dir.create(PARTS, showWarnings = FALSE)
if (nzchar(Sys.getenv("COMBINE"))) {
  # rd(p): row-bind every file in PARTS whose name matches p
  rd <- function(p) do.call(rbind, lapply(list.files(PARTS, p, full.names = TRUE), read.csv))
  write.csv(rd("_curve\\.csv$"), file.path(TABLE_DIR, "17c_scvi_rods_technical_curves.csv"), row.names = FALSE)
  v <- rd("_verdict\\.csv$"); write.csv(v, file.path(TABLE_DIR, "17c_scvi_rods_technical_verdicts.csv"), row.names = FALSE)
  print(v); quit(save = "no")
}

# ---- Train (or load) five scVI embeddings per condition, then R(k) and the verdict ----
if (nzchar(Sys.getenv("COND"))) cond <- cond[Sys.getenv("COND")]
curves <- list(); verd <- list()
for (nm in names(cond)) {
  message("== ", nm, " (", ncol(cond[[nm]]$obj), " cells)")
  # SEED_ONLY=<s> (with COND): train and save that one embedding (seeds can run as parallel jobs);
  # otherwise cached seeds are loaded and missing ones trained here
  sf <- function(s) file.path(PARTS, sprintf("%s_seed%s.rds", nm, s)); so <- Sys.getenv("SEED_ONLY")
  if (nzchar(so)) { reticulate::import("torch")$set_num_threads(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
    saveRDS(sccb_embed(cond[[nm]]$obj, cond[[nm]]$emb, as.integer(so), batch = "sample"), sf(so)); quit(save = "no") }
  lat <- lapply(setNames(SCDHA_SEED_LIST, SCDHA_SEED_LIST), function(s)
    if (file.exists(sf(s))) readRDS(sf(s)) else sccb_embed(cond[[nm]]$obj, cond[[nm]]$emb, s, batch = "sample"))
  if (any(vapply(lat, is.null, TRUE))) stop("scVI unavailable: set RETICULATE_PYTHON")
  cv <- sccb_reproducibility(cond[[nm]]$obj, k_range = KS, latents = lat)
  # does this condition's seed-1 k=3 partition still recover the 13t minority axis?
  l3 <- setNames(sccb_louvain_to_k(lat[[1]], 3, SEED), rownames(lat[[1]])); sh <- intersect(names(l3), names(ref_axis))
  cv$ari_vs_13t_axis <- round(ARI(as.integer(factor(l3[sh])), as.integer(factor(ref_axis[sh]))), 3)
  cv$ari_vs_13t_k3   <- round(ARI(as.integer(factor(l3[sh])), as.integer(factor(ref_k3[sh]))), 3)
  curves[[nm]] <- cbind(condition = nm, n_cells = ncol(cond[[nm]]$obj), cv)
  verd[[nm]]   <- cbind(condition = nm, sccb_verdict(cv, 0.8), ari_vs_13t_axis = cv$ari_vs_13t_axis[1])
  # Per-condition parts are always written, so COMBINE can merge runs made in separate jobs
  write.csv(curves[[nm]], file.path(PARTS, paste0(nm, "_curve.csv")), row.names = FALSE)
  write.csv(verd[[nm]], file.path(PARTS, paste0(nm, "_verdict.csv")), row.names = FALSE)
  if (!nzchar(Sys.getenv("COND"))) write.csv(do.call(rbind, curves), file.path(TABLE_DIR, "17c_scvi_rods_technical_curves.csv"), row.names = FALSE)
}
if (nzchar(Sys.getenv("COND"))) quit(save = "no")
v <- do.call(rbind, verd)
write.csv(v, file.path(TABLE_DIR, "17c_scvi_rods_technical_verdicts.csv"), row.names = FALSE); print(v)
