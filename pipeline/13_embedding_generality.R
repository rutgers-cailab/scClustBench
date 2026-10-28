# ============================================================
# pipeline/13_embedding_generality.R
#
# Computes the Stage B curve R(k) for bipolar cells and rods under several embedding families,
# so the verdict can be checked for dependence on the embedding. Each family is a set of
# embeddings, each Harmony-corrected on sample (scVI uses sample as its batch key):
#   scdha        the five scDHA embeddings saved by step 07 (retrained here if missing and scDHA
#                is installed; otherwise the family is skipped)
#   pca_perturb  five seeds, each a PCA on a random 80% of the HVGs with a seeded randomized PCA
#   pca_hyper    nine PCAs over HVGs {2000, 3000, 5000} x PCs {20, 30, 50}
#   scvi         five scVI seeds (only without POP, and only if RETICULATE_PYTHON points to a
#                Python with scvi-tools; with POP the scVI family is step 13s)
# R(k) is the mean pairwise ARI between Louvain partitions of the embeddings in a family.
#
# Paper:    Stage B embedding-family verdicts (via pipeline/13c_combine.R and 12_decision_rule.R)
# Inputs:   output/rds/07_bc_positive_control.rds, output/rds/07_rods_positive_control.rds (step 07);
#           if missing, the population is loaded from output/rds/02_*.rds;
#           without POP also output/tables/16_null_calibration_tau.csv (step 16; optional)
# Outputs:  with POP:    output/tables/13_embedding_generality_curves_<pop>.csv
#           without POP: output/tables/13_embedding_generality_curves.csv
#                        output/tables/13_embedding_generality_verdicts.csv
#                        output/figures/publication/13_embedding_generality.{pdf,png}
# Run:      Rscript run_pipeline.R 13_bc 13_rods   (POP=BC and POP=Rods; then 13s and 13c)
#           or  POP=BC Rscript pipeline/13_embedding_generality.R
#           Without POP (Rscript pipeline/13_embedding_generality.R) both populations run in one
#           process, scVI included, and the combined tables are written directly.
# Options:  POP   "BC" or "Rods": run one population, without the scVI family (default: both)
# Runtime:  several hours (the Louvain resolution search dominates), plus about 1 h if scDHA is
#           retrained and, without POP, 1 h (GPU) to 3 h (CPU) for scVI. About 32 GB RAM for rods.
# ============================================================
source("config.R"); source("R/scClustBench.R")
suppressPackageStartupMessages(library(ggplot2))
POPS <- c(BC = 18, Rods = 10)          # largest k per population (BC as in step 07b)
# POP=<BC|Rods>: run one population and write 13_embedding_generality_curves_<pop>.csv (populations run as
# parallel jobs; 13c_combine.R merges the parts with the scVI parts of 13s).
# Without POP: both populations, combined outputs.
if (nzchar(Sys.getenv("POP"))) POPS <- POPS[Sys.getenv("POP")]

curves <- list()
for (ct in names(POPS)) {
  ks <- 2:POPS[[ct]]; tag <- tolower(ct)
  message("=== ", ct, " ===")
  # Reuse the population and scDHA embeddings saved by step 07 when available
  rds <- file.path(RDS_DIR, sprintf("07_%s_positive_control.rds", tag))
  pc  <- if (file.exists(rds)) readRDS(rds) else NULL
  obj <- sccb_prep(if (is.null(pc)) load_population(ct) else pc$bc, "LogNormalize", 3000)
  fam <- list()
  fam$scdha <- if (!is.null(pc)) sccb_reproducibility(obj, latents = pc$seed_lat, k_range = ks) else
    sccb_reproducibility(obj, "scdha", SCDHA_SEED_LIST, ks, batch = "sample")
  fam$pca_perturb <- sccb_reproducibility(obj, "pca_perturb", SCDHA_SEED_LIST, ks, batch = "sample", npcs = N_PCS_CLUSTERING)
  # pca_hyper: one PCA per (number of HVGs, number of PCs), all with the same seed
  hyper <- list()
  for (nf in c(2000, 3000, 5000)) { o <- sccb_prep(obj, "LogNormalize", nf)
    for (np in c(20, 30, 50)) hyper[[sprintf("hvg%d_pc%d", nf, np)]] <- sccb_embed(o, "pca", SEED, "sample", npcs = np) }
  fam$pca_hyper <- sccb_reproducibility(obj, latents = hyper, k_range = ks)
  # With POP set, the scVI family is left to step 13s
  if (!nzchar(Sys.getenv("POP"))) fam$scvi <- sccb_reproducibility(obj, "scvi", SCDHA_SEED_LIST, ks, batch = "sample")
  # Families whose backend is unavailable return NULL and are dropped
  for (f in names(fam)) if (!is.null(fam[[f]])) curves[[paste(ct, f)]] <- cbind(celltype = ct, embedding_family = f, fam[[f]])
}
cur <- do.call(rbind, curves); rownames(cur) <- NULL
if (nzchar(Sys.getenv("POP"))) {
  write.csv(cur, file.path(TABLE_DIR, sprintf("13_embedding_generality_curves_%s.csv", tolower(Sys.getenv("POP")))), row.names = FALSE)
  quit(save = "no")
}
# ---- Verdict per (population, family), tables and the summary figure ----
# tau: 0.8, raised to the step 16 null calibration of the family when that is higher (as in step 12)
calib_f <- file.path(TABLE_DIR, "16_null_calibration_tau.csv"); calib <- if (file.exists(calib_f)) read.csv(calib_f) else NULL
ver <- do.call(rbind, lapply(split(cur, list(cur$celltype, cur$embedding_family), drop = TRUE), function(d)
  cbind(celltype = d$celltype[1], embedding_family = d$embedding_family[1], sccb_verdict(d, tau = sccb_tau(d$embedding_family[1], calib)))))
write.csv(cur, file.path(TABLE_DIR, "13_embedding_generality_curves.csv"), row.names = FALSE)
write.csv(ver, file.path(TABLE_DIR, "13_embedding_generality_verdicts.csv"), row.names = FALSE)
print(ver)

# Plot only k values that the reference partition actually reached
p <- ggplot(cur[cur$actual_k == cur$k, ], aes(k, mean_cross_seed_ari, colour = embedding_family)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  geom_line() + geom_point(size = 1.5) + facet_wrap(~celltype, scales = "free_x") +
  scale_y_continuous(limits = c(NA, 1)) +
  labs(x = "k", y = "R(k): mean pairwise ARI across embeddings", colour = "Embedding family",
       title = "Across-embedding reproducibility is not scDHA-specific") + pub_theme()
save_fig(p, "13_embedding_generality", FIG_PUB, width = 8, height = 3.8)
