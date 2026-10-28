# ============================================================
# pipeline/13s_scvi_family.R
#
# Builds the scVI embedding family for Stage B: five scVI models trained with
# different seeds (sample as the batch key), each clustered at every k, and
# R(k) computed as the mean pairwise ARI between the five partitions (bipolar
# cells k = 2 to 18, rods k = 2 to 10). With POP set (the mode used by the
# pipeline) it handles one population and leaves the merge to 13c_combine.R.
# Without POP it runs both populations in one process and merges the scVI
# curves into the stage 13 tables and figure itself.
#
# Paper:    Stage B scVI family (Supplementary Figure S4, via 13c_combine.R and 13f)
# Inputs:   output/rds/07_{bc,rods}_positive_control.rds (from 07_positive_control_bc.R)
#           output/tables/13_embedding_generality_curves.csv (standalone mode only, if present)
# Outputs:  POP mode: output/rds/13s_scvi_<pop>_seed<s>.rds (one per seed),
#           output/rds/13s_scvi_<pop>_latents.rds (all five, reused by 13t),
#           output/tables/13s_scvi_curves_<pop>.csv
#           Standalone mode: output/tables/13s_scvi_curves.csv, plus updated
#           output/tables/13_embedding_generality_{curves,verdicts}.csv and
#           output/figures/publication/13_embedding_generality.{pdf,png}
# Run:      Rscript run_pipeline.R 13s_bc 13s_rods
#           or  POP=Rods RETICULATE_PYTHON=~/micromamba/envs/scvi/bin/python Rscript pipeline/13s_scvi_family.R
# Options:  POP                  BC or Rods; one population (default: unset, standalone mode)
#           SEED_ONLY            train and save only this seed, then exit (POP mode; default: unset)
#           FORCE                standalone mode: retrain even if 13s_scvi_curves.csv exists (default: unset)
#           RETICULATE_PYTHON    python with scvi-tools (required)
#           SLURM_CPUS_PER_TASK  torch threads in POP mode (default 8)
# ============================================================
source("config.R"); source("R/scClustBench.R")
suppressPackageStartupMessages(library(ggplot2))
# Largest k swept for each population
POPS <- c(BC = 18, Rods = 10)

# ---- POP mode: one population, seeds trained separately and cached ----
# POP=<BC|Rods>: one population; its five scVI embeddings are saved (output/rds/13s_scvi_<pop>_latents.rds,
# reused by 13t) and its curve is written to 13s_scvi_curves_<pop>.csv for 13c_combine.R.
if (nzchar(Sys.getenv("POP"))) {
  ct <- Sys.getenv("POP"); tag <- tolower(ct)
  pc  <- readRDS(file.path(RDS_DIR, sprintf("07_%s_positive_control.rds", tag)))
  obj <- sccb_prep(pc$bc, "LogNormalize", 3000)
  # SEED_ONLY=<s>: train that seed only (seeds can run as parallel jobs) and save it; otherwise
  # train any seed not yet saved, then load all of them.
  # On CPU, R's BLAS must run single-threaded (OPENBLAS_NUM_THREADS=1, set by run_scvi_cpu.sbatch)
  # or it competes with torch for the cores and training becomes many times slower.
  reticulate::import("torch")$set_num_threads(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
  # sf(s): path of the cached embedding for seed s
  sf <- function(s) file.path(RDS_DIR, sprintf("13s_scvi_%s_seed%s.rds", tag, s))
  so <- Sys.getenv("SEED_ONLY")
  for (s in SCDHA_SEED_LIST) if ((so == "" && !file.exists(sf(s))) || so == as.character(s)) {
    L <- sccb_embed(obj, "scvi", s, batch = "sample"); if (is.null(L)) stop("scVI embedding unavailable"); saveRDS(L, sf(s)) }
  if (nzchar(so)) quit(save = "no")
  lats <- lapply(setNames(SCDHA_SEED_LIST, SCDHA_SEED_LIST), function(s) readRDS(sf(s)))
  saveRDS(lats, file.path(RDS_DIR, sprintf("13s_scvi_%s_latents.rds", tag)))
  r <- sccb_reproducibility(obj, latents = lats, k_range = 2:POPS[[ct]])
  write.csv(cbind(celltype = ct, embedding_family = "scvi", r), file.path(TABLE_DIR, sprintf("13s_scvi_curves_%s.csv", tag)), row.names = FALSE)
  quit(save = "no")
}

# ---- Standalone mode: both populations in one process ----
# Reuses 13s_scvi_curves.csv if it exists (set FORCE to retrain)
cf <- file.path(TABLE_DIR, "13s_scvi_curves.csv")
if (!file.exists(cf) || nzchar(Sys.getenv("FORCE"))) {
  curves <- list()
  for (ct in names(POPS)) {
    message("=== ", ct, " (scvi) ===")
    # Use the population saved by step 07 if available, otherwise rebuild it from the annotated object
    rds <- file.path(RDS_DIR, sprintf("07_%s_positive_control.rds", tolower(ct)))
    pc  <- if (file.exists(rds)) readRDS(rds) else NULL
    obj <- sccb_prep(if (is.null(pc)) load_population(ct) else pc$bc, "LogNormalize", 3000)
    r <- sccb_reproducibility(obj, "scvi", SCDHA_SEED_LIST, 2:POPS[[ct]], batch = "sample")
    if (is.null(r)) stop("scVI embedding unavailable")
    curves[[ct]] <- cbind(celltype = ct, embedding_family = "scvi", r)
  }
  write.csv(do.call(rbind, curves), cf, row.names = FALSE)
}

# ---- Merge the scVI curves into the stage 13 tables and redraw the stage 13 figure ----
sc <- read.csv(cf)
f13 <- file.path(TABLE_DIR, "13_embedding_generality_curves.csv")
# tau: 0.8, raised to the step 16 null calibration of the family when that is higher (as in step 12)
calib_f <- file.path(TABLE_DIR, "16_null_calibration_tau.csv"); calib <- if (file.exists(calib_f)) read.csv(calib_f) else NULL
if (file.exists(f13)) {
  # Replace any earlier scVI rows, then recompute the verdict for every population x family
  cur <- read.csv(f13); cur <- rbind(cur[cur$embedding_family != "scvi", ], sc[, names(cur)])
  ver <- do.call(rbind, lapply(split(cur, list(cur$celltype, cur$embedding_family), drop = TRUE), function(d)
    cbind(celltype = d$celltype[1], embedding_family = d$embedding_family[1], sccb_verdict(d, tau = sccb_tau(d$embedding_family[1], calib)))))
  write.csv(cur, f13, row.names = FALSE)
  write.csv(ver, file.path(TABLE_DIR, "13_embedding_generality_verdicts.csv"), row.names = FALSE)
  print(ver)
  # Working figure (not a paper figure); only k values where the reference embedding reached k are drawn
  p <- ggplot(cur[cur$actual_k == cur$k, ], aes(k, mean_cross_seed_ari, colour = embedding_family)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
    geom_line() + geom_point(size = 1.5) + facet_wrap(~celltype, scales = "free_x") +
    scale_y_continuous(limits = c(NA, 1)) +
    labs(x = "k", y = "R(k): mean pairwise ARI across embeddings", colour = "Embedding family",
         title = "Across-embedding reproducibility is not scDHA-specific") + pub_theme()
  save_fig(p, "13_embedding_generality", FIG_PUB, width = 8, height = 3.8)
} else message("13_embedding_generality_curves.csv not found yet; rerun this script after stage 13 to merge.")
