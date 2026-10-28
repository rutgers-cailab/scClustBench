# ============================================================
# pipeline/14_preprocessing_sensitivity.R
#
# Tests whether the benchmark and the decision rule depend on preprocessing.
# Factorial design: normalisation {LogNormalize, SCTransform, scran pooling} x
# kNN graph distance {euclidean, cosine (FindNeighbors annoy.metric)}.
#   (a) Whole retina (Phase 1): Harmony-PCA + Louvain at k = 2 to 12, the five
#       label-free metrics, Stage A selection, and post hoc ARI/NMI against the
#       step 02 cell_type annotation.
#   (b) Bipolar cells and rods: R(k) from the PCA gene-resampling family (five
#       seeds) with the same Louvain/distance variant, and the Stage B verdict.
#   (c) Combines the tagged tables written by parallel (a) and (b) jobs.
#
# Paper:    Preprocessing sensitivity analysis (supplementary tables)
# Inputs:   output/rds/02_annotated.rds, output/rds/02_rods.rds (via load_population())
#           PART=c: output/tables/14_preprocessing_{phase1_scores,curves}_<TAG>.csv
# Outputs:  output/tables/14_preprocessing_phase1_scores.csv    (all partitions)
#           output/tables/14_preprocessing_phase1_selection.csv (Stage A top pick per variant)
#           output/tables/14_preprocessing_curves.csv           (R(k), bipolar cells and rods)
#           output/tables/14_preprocessing_verdicts.csv
#           output/figures/publication/14_preprocessing_sensitivity.{pdf,png}
#           With TAG set, (a) and (b) write <table>_<TAG>.csv instead.
# Run:      Rscript run_pipeline.R 14 14_tables
#           or  PART=b CELLTYPES=Rods NORMS=SCT TAG=b_Rods_SCT Rscript pipeline/14_preprocessing_sensitivity.R
# Options:  PART       which parts to run, any of a, b, c (default ab)
#           NORMS      comma list of normalisations (default LogNormalize,SCT,scran)
#           DISTS      comma list of kNN distances (default euclidean,cosine)
#           CELLTYPES  comma list of populations for (b) (default BC,Rods)
#           TAG        suffix for the output tables of (a) and (b) (default: none)
#           SCCB_BENCH_WORKERS  parallel workers in the benchmark grid (default 6)
# Runtime:  about 10 to 14 h in one process, dominated by (a): 6 variants x 11 k on
#           about 56k cells with ROGUE, FindAllMarkers and 3 Jaccard bootstraps each;
#           (b) about 2 h. RAM about 64 GB (dense counts for ROGUE, SCTransform).
# ============================================================
source("config.R"); source("R/scClustBench.R")
suppressPackageStartupMessages(library(ggplot2))

# ---- Settings from the environment ----
# env_list(v, d): environment variable v (default d) split on commas into a character vector
env_list <- function(v, d) strsplit(Sys.getenv(v, d), ",")[[1]]
NORMS <- env_list("NORMS", "LogNormalize,SCT,scran"); DISTS <- env_list("DISTS", "euclidean,cosine")
CELLTYPES <- env_list("CELLTYPES", "BC,Rods"); TAG <- Sys.getenv("TAG")
# out(f): path in TABLE_DIR for table f, with _<TAG> inserted before .csv when TAG is set
out <- function(f) file.path(TABLE_DIR, if (nzchar(TAG)) sub("\\.csv$", paste0("_", TAG, ".csv"), f) else f)
PART  <- Sys.getenv("PART", "ab")
# The cosine variant is registered as another clustering method, so the benchmark and
# R(k) code run it unchanged
SCCB_METHODS$louvain_cosine <- function(latent, k, seed) sccb_louvain_to_k(latent, k, seed, annoy.metric = "cosine")
SCCB_METHOD_PATHS$louvain_cosine <- function(latent, seed) sccb_louvain_path(latent, seed, annoy.metric = "cosine")
meth <- c(euclidean = "louvain", cosine = "louvain_cosine")

# ---- (a) Phase 1, whole retina ----
if (grepl("a", PART)) {
  all <- load_population("ALL"); truth <- as.character(all$cell_type); re <- sccb_rogue_expr(all)
  sc <- list()
  for (nm in NORMS) {
    # One Harmony-PCA embedding per normalisation, shared by both distances
    o <- sccb_prep(all, nm, 3000)
    L <- sccb_embed(o, "pca", SEED, "sample", npcs = N_PCS_CLUSTERING)
    for (d in DISTS) {
      message("=== Phase 1: ", nm, " x ", d, " ===")
      s <- sccb_benchmark(o, "pca", meth[[d]], POC_K_RANGE, n_boot = BENCH_JACCARD_BOOTSTRAPS,
                          truth = truth, latents = list(pca = L), rogue = re)
      sc[[paste(nm, d)]] <- cbind(norm = nm, distance = d, s)
      # Written after every variant so a stopped job keeps its finished variants
      write.csv(do.call(rbind, sc), out("14_preprocessing_phase1_scores.csv"), row.names = FALSE)
    }
  }
  sel <- do.call(rbind, lapply(sc, function(s) sccb_select(s)[1, ]))
  write.csv(sel, out("14_preprocessing_phase1_selection.csv"), row.names = FALSE)
  print(sel[, c("norm", "distance", "target_k", "score", "ari", "nmi")])
}

# ---- (b) Stage B verdicts, bipolar cells and rods ----
if (grepl("b", PART)) {
  curves <- list()
  for (ct in CELLTYPES) {
    pop <- load_population(ct); ks <- if (ct == "BC") 2:18 else 2:10
    for (nm in NORMS) {
      o <- sccb_prep(pop, nm, 3000)
      # Five gene-resampled PCA embeddings, reused for both distances
      lat <- lapply(setNames(SCDHA_SEED_LIST, SCDHA_SEED_LIST), function(s)
        sccb_embed(o, "pca_perturb", s, "sample", npcs = N_PCS_CLUSTERING))
      for (d in DISTS) {
        message("=== ", ct, ": ", nm, " x ", d, " ===")
        curves[[paste(ct, nm, d)]] <- cbind(celltype = ct, norm = nm, distance = d,
                                            sccb_reproducibility(o, latents = lat, k_range = ks, method = meth[[d]]))
      }
    }
  }
  cur <- do.call(rbind, curves); rownames(cur) <- NULL
  ver <- do.call(rbind, lapply(curves, function(d) cbind(d[1, c("celltype", "norm", "distance")], sccb_verdict(d))))
  write.csv(cur, out("14_preprocessing_curves.csv"), row.names = FALSE)
  write.csv(ver, out("14_preprocessing_verdicts.csv"), row.names = FALSE)
  print(ver)
}

# ---- (c) Combine the tagged outputs of parallel jobs into the untagged tables ----
if (grepl("c", PART)) {
  # rb(pat): row-bind every table in TABLE_DIR whose name matches pat
  rb <- function(pat) do.call(rbind, lapply(list.files(TABLE_DIR, pat, full.names = TRUE), read.csv))
  sc <- rb("^14_preprocessing_phase1_scores_.+\\.csv$"); sc <- sc[!duplicated(sc[, c("norm", "distance", "method", "target_k")]), ]
  write.csv(sc, file.path(TABLE_DIR, "14_preprocessing_phase1_scores.csv"), row.names = FALSE)
  # Stage A selection is rerun per variant on the combined scores
  sel <- do.call(rbind, lapply(split(sc, paste(sc$norm, sc$distance)), function(s) sccb_select(s)[1, ]))
  write.csv(sel, file.path(TABLE_DIR, "14_preprocessing_phase1_selection.csv"), row.names = FALSE)
  # Verdicts are recomputed from the combined curves
  cur <- rb("^14_preprocessing_curves_.+\\.csv$")
  ver <- do.call(rbind, lapply(split(cur, paste(cur$celltype, cur$norm, cur$distance)), function(d)
    cbind(d[1, c("celltype", "norm", "distance")], sccb_verdict(d))))
  write.csv(cur, file.path(TABLE_DIR, "14_preprocessing_curves.csv"), row.names = FALSE)
  write.csv(ver, file.path(TABLE_DIR, "14_preprocessing_verdicts.csv"), row.names = FALSE)
  print(sel[, intersect(c("norm", "distance", "target_k", "score", "ari", "nmi"), names(sel))]); print(ver)
}

# ---- Working figure of R(k) per variant (whenever curves are available) ----
if (grepl("b|c", PART) && exists("cur")) {
  p <- ggplot(cur[cur$actual_k == cur$k, ], aes(k, mean_cross_seed_ari, colour = norm, linetype = distance)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") + geom_line() + geom_point(size = 1.2) +
    facet_wrap(~celltype, scales = "free_x") +
    labs(x = "k", y = "R(k) (pca_perturb)", colour = "Normalisation", linetype = "kNN distance") + pub_theme()
  save_fig(p, "14_preprocessing_sensitivity", FIG_PUB, width = 8, height = 3.8)
}
